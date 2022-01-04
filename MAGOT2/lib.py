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
    """
    Reads a reasonably sane gff3 or gtf file and fills and returns a nested dictionary containing annotation information
    (format feature_type_interval_tree -> dictionary_of_feature_types -> dictionary_of_transcripts -> dictionary_of_genes)
    attributes of features are stored as a dictionary within each interval

    :param gfffile: Path. Path to gff3 or gtf
    :param version: str or int. Valid options are 2, 3, or "auto" (in which case version will be determined from first line)
    """
    annotdict = OrderedDict()
    t2g = {}
    with open(gfffile) as gff:
        for i,line in enumerate(gff):
            fields = line.strip().split('\t')
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
                attrs['source'] = fields[1]
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
        pysam.faidx(str(fasta_file))
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

def annot2gtf(annotdict: dict, features2write: Union[list,str] = 'CDS') -> str:
    """
    takes a MAGOT2 format annotation dictionary and returns a string for specific annotations in gtf format

    :params annotdict: dict. MAGOT2 format annotation dictionary (nested objects coordinate_IntervalTree -> feature_dictionary -> \
transcript_dictionary -> gene_dictionary)
    :params features2write: list or str default "CDS". What features to write out (usually "CDS", sometimes also "exon")
    """
    if type(features2write) == str:
        f2w = [features2write]
    else:
        f2w = features2write
    linelist = []
    for gene_id in annotdict:
        for transcript_id in annotdict[gene_id]:
            for feature in f2w:
                if feature in annotdict[gene_id][transcript_id]:
                    for ivl in annotdict[gene_id][transcript_id][feature]:
                        attrs = ivl[2]
                        linelist.append('\t'.join([attrs['seqid'],attrs['source'],feature,str(ivl[0] + 1), str(ivl[1]), 
                                                attrs['score'],attrs['strand'],attrs['phase'],
                                                'gene_id ' + gene_id + ';transcript_id ' + transcript_id]))
    return "\n".join(linelist)

def annot2gff3(annotdict: dict) -> str:
    """
    takes a MAGOT2 format annotation dictionary and returns a string for gene, mRNA/transcript, exon, and CDS features in gff3 format

    :params annotdict: dict. MAGOT2 format annotation dictionary (nested objects coordinate_IntervalTree -> feature_dictionary -> \
transcript_dictionary -> gene_dictionary)
    """
    linelist = []
    for gene_id in annotdict:
        tlines = []
        seqid,gstart,gend,strand,gscore,source = '.',np.inf,0,'.',0,'.'
        for transcript_id in annotdict[gene_id]:
            flines = []
            tstart,tend,tscore = np.inf,0,0
            for feature in ['exon','CDS']:
                if feature in annotdict[gene_id][transcript_id]:
                    for i,ivl in enumerate(annotdict[gene_id][transcript_id][feature]):
                        attrs = ivl[2]
                        if ivl[0] - 1 < tstart: tstart = ivl[0] -1
                        if ivl[1] > tend: tend = ivl[1]
                        if attrs['score'] > tscore: tscore = attrs['score']
                        if seqid == '.': seqid = attrs['seqid']
                        if strand == '.': strand = attrs['strand']
                        if source == '.': source = attrs['source']
                        flines.append('\t'.join([attrs['seqid'],attrs['source'],feature,str(ivl[0] + 1), str(ivl[1]), 
                                                attrs['score'],attrs['strand'],attrs['phase'],
                                                'ID=' + transcript_id + ':' + feature + str(i) + ';Parent=' + transcript_id]))
            if 'CDS' in  annotdict[gene_id][transcript_id]:
                ttype = 'mRNA'
            elif 'exon' in annotdict[gene_id][transcript_id]:
                ttype = 'transcript'
            else:
                continue
            tlines.append(
                '\t'.join([seqid,source,ttype,str(tstart),str(tend), str(tscore),strand,'.',
                        'ID=' + transcript_id + ';Parent=' + gene_id])
            )
            tlines.extend(flines)
            if tstart < gstart: gstart = tstart
            if tend > gend: gend = tend
            if tscore > gscore: gscore = tscore
        if len(tlines) > 0:
            linelist.append(
                '\t'.join([seqid,source,'gene',str(gstart),str(gend), str(gscore),strand,'.',
                        'ID=' + gene_id])
            )
    return "\n".join(linelist)

def dipvcf2tripvcf(vcfin: Path,vcfout: Path,expand_on: str = "SR"):
    """
    reads a diploid vcf file and writes a triploid vcf file, chosing which allele to duplicate based on the "expand_on" param.

    :param vcfin: Path. Path to diploid input vcf file
    :param vcfout: Path. Path to output triploid vcf file
    :param expand_on: str defualt "SR". What information to use to chose which allele to duplicate. Currently only valid option \
is "SR", which will use the "SR" attribute from the info collumn to add the allele with the highest depth to the alleles field
    """
    with open(vcfin) as dip:
        with open(vcfout,'w') as trip:
            for line in dip:
                if line[0] == "#":
                    trip.write(line)
                    continue
                elif ';SR=' in line:
                    srs = line.split(';SR=')[1].split(';')[0].split('\t')[0]
                elif '\tSR=' in line:
                    srs = line.split('\tSR=')[1].split(';')[0].split('\t')[0]
                else:
                    continue
                srs = [int(k) for k in srs.split(',')]
                ads = [srs[i] + srs[i+1] for i in range(0,len(srs),2)]
                allele = np.argmax(ads)
                fields = line.strip().split('\t')
                subfields = fields[9].split(':')
                gtindex = fields[8].split(':').index('GT')
                alleles = [int(k) for k in subfields[gtindex].split('/')] + [allele]
                alleles.sort()
                subfields[gtindex] = '/'.join([str(k) for k in alleles])
                fields[9] = ":".join(subfields)
                trip.write("\t".join(fields) + '\n')

