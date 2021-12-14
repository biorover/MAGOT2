#!/usr/bin/env python
import sys, pysam, os
import numpy as np
from intervaltree import IntervalTree, Interval
from collections import OrderedDict
from typing import Union
from pathlib import Path

def mode(a,prec = 2) -> float:
    vals,counts =  np.unique(np.around(a,prec),return_counts=True)
    return vals[np.argmax(counts)]

def translate(x: str,padXs: bool = False,frames: Union[list,int] = 1) -> Union[list,str]:
    """
    Translates DNA/RNA sequence into amino acid characters

    :params x: string dna sequence
    """
    ttabs = []
    codons = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    if 'u' in x or 'U' in x:
        x = x.replace('U','T').replace('u','t')
    rstring = False
    if type(frames) == int:
        frames = [frames]
        rstring = True
    elif frames in ['all']:
        frames = [1,2,3,-1,-2,-3]
    for frame in frames:
        ttab = []
        if frame >= 0:
            seq = x
            start = frame - 1
        else:
            seq = revcomp(x)
            start = -1 * frame - 1
        for i in range(start,len(seq),3):
            codon = x[i:i+3]
            if len(codon) == 3 or padXs:
                if codon.upper() in codons:
                    ttab.append(codons[codon])
                else:
                    ttab.append('X')
        ttabs.append("".join(ttab))
    if rstring:
        return ttabs[0]
    else:
        return ttabs

def revcomp(x: str) -> str:
    comp = str.maketrans('ATCGMKRYVBHDatcgmkryvbhd','TAGCKMYRBVDHtagckmyrbvdh')
    seq = x[::-1].translate(comp)
    return seq

def orfs(x: str,best: bool = False) -> Union[list,str]: #todo: #1 add info table
    orfs = []
    for frames in translate(x,frames = 'all'):
        orfs += frames.split('*')
    if best:
        lens = [len(k) for k in orfs]
        return orfs[np.argmax(lens)]
    else:
        return orfs

def read_gff(gfffile: Path,version: Union[str,int] = 'auto') -> OrderedDict:
    annotdict = OrderedDict()
    t2g = {}
    with open(gfffile) as gff:
        for i,line in enumerate(gff):
            fields = line.split('\t')
            if version == 'auto' and i==0:
                if "#gff-version" in line:
                    version = int(line.split()[1])
                elif "ID=" in line:
                    version = 3
                elif "gene_id " in line:
                    version = 2
                else:
                    sys.stderr.write('Cannot determine gff format version of file from first line\n')
                    return None
            if line[0] != "#":
                feature = fields[2]
                if version == 3:
                    attrs = {k:v for (k,v) in [i.split('=') for i in fields[8].split(';')] }
                elif version == 2:
                    attrs = {k:v.replace('"','') for (k,v) in [(i.split(' ')[0],' '.join(i.split(' ')[1:])) for i in fields[8].replace('; ',';').split(';')] }
                if version == 3:
                    if feature == 'gene':
                        gene_id = attrs['ID']
                        transcript_id = 'none' # this is so I can access the coords of the whole gene under annotdict[gene_id]["none"]['gene']
                    elif feature in ['mRNA','transcript']:
                        transcript_id = attrs['ID']
                        gene_id = attrs['Parent']
                        t2g[transcript_id] = gene_id
                    else:
                        transcript_id = attrs['Parent']
                        gene_id = t2g[transcript_id]
                elif version == 2:
                    gene_id = attrs['gene_id']
                    transcript_id = attrs['transcript_id']
                if not gene_id in annotdict:
                    annotdict[gene_id] = OrderedDict()
                if not transcript_id in annotdict[gene_id]:
                    annotdict[gene_id][transcript_id] = dict()
                if not feature in annotdict[gene_id][transcript_id]:
                    annotdict[gene_id][transcript_id][feature] = IntervalTree()
                start = int(fields[3]) - 1
                stop = int(fields[4])
                attrs['seqid'] = fields[0]
                attrs['score'] = fields[5]
                attrs['strand'] = fields[6]
                attrs['phase'] = fields[7]
                annotdict[gene_id][transcript_id][feature][start:stop] = attrs
    return annotdict

def annot2seqs(annotdict: dict, fasta_file: Path, which_transcript: str = 'all', 
            seq_from: str = 'CDS', seq_level: str = 'transcript', seq_type: str = "nucl") -> OrderedDict:
    """
    takes an annotdict (nested objects coordinate_IntervalTree -> feature_dictionary -> \
transcript_dictionary -> gene_dictionary) ant returns a dictionary of all sequences

    :param annotdict: dict. All annotations
    :param fasta_file: Path. Path to fasta file
    :param which_transcript: str default "all". Options are "all", "longest", "first", and "best_scoring"    
    :param seq_from: str default "CDS". Name of feature to extract sequence from (usually "CDS" or "exon")
    :param seq_level: str default "transcript". What to return the sequences of - usually "transcript" or "gene" \
(gene only valid if which_transcript = "first", "longest", or "best_scoring"). \
Could also be "CDS" or "exon" or other sub-feature (must then match "seq_from" argument)
    :param seq_type: str default "nucl". Whether to output "nucl" (nucleotide), "aa" (translated amino acid) \
or "lorfaa" (longest orf amino acid) sequence
    """
    if not os.path.exists(str(fasta_file) + '.fai'):
        pysam.faidx(fasta_file)
    if not seq_level in ['transcript','gene']:
        if seq_level != seq_from:
            sys.stderr.write('error: invalid choice for "seq_level": ' + str(name_from) + 
                            '! Currently supported: "transcript", "gene", or else must match "seq_from" argument\n')
            return None
    fasta = pysam.FastaFile(fasta_file)
    outseqs = OrderedDict()
    for gene_id in annotdict:
        tseqs = []
        tscores = []
        for transcript_id in annotdict[gene_id]:
            if seq_from in annotdict[gene_id][transcript_id]:
                seq = []
                scores = []
                for i,interval in enumerate(sorted(annotdict[gene_id][transcript_id][seq_from])):
                    if i == 0:
                        seqid = interval[2]['seqid']
                        strand = interval[2]['strand']
                    elif seqid != interval[2]['seqid']:
                        sys.stderr.write('error: annotations contain transcript with features on different sequence! Not currently supported\n')
                        return None
                    elif strand != interval[2]['strand']:
                        sys.stderr.write('error: annotations contain transcript with features on different strands! Not currently supported, and doesn\'t really make sense!\n')
                        return None
                    iseq = fasta.fetch(seqid,interval[0],interval[1])
                    if not seq_level in ['transcript','gene']:
                        outseqs[transcript_id + ":" + seq_from + str(i)] = iseq
                    else:
                        seq.append(iseq)
                        if which_transcript == 'best_scoring':
                            scores.append(interval[2]['score'])
                if seq_level in ['transcript','gene']:
                    seq = "".join(seq)
                    if strand == '-':
                        seq = revcomp(seq)
                    if which_transcript == 'all':
                        outseqs[transcript_id] = seq
                    else:
                        tseqs.append(seq[:])
                    if which_transcript == 'best_scoring':
                        tscores.append(np.mean(scores))
        if which_transcript == 'first' and seq_level in ['transcript','gene']:
            outseqs[locals()[seq_level + '_id'] ] = tseqs[0]
        elif which_transcript == 'longest' and seq_level in ['transcript','gene']:
            seqlens = [len(k) for k in tseqs]
            outseqs[locals()[seq_level + '_id'] ] = tseqs[np.argmax(seqlens)]
        elif which_transcript == 'best_scoring' and seq_level in ['transcript','gene']:
            outseqs[locals()[seq_level + '_id'] ] = tseqs[np.argmax(tscores)]
    if seq_type == 'aa':
        outseqs = OrderedDict((k,translate(v)) for k,v in outseqs.items())
    elif seq_type == 'lorfaa':
        outseqs = OrderedDict((k,orfs(v,best = True)) for k,v in outseqs.items())
    return outseqs
