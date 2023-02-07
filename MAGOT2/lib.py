#!/usr/bin/env python
import sys, pysam, os, re, glob
import numpy as np
from intervaltree import IntervalTree, Interval
from collections import OrderedDict
from typing import Union
from pathlib import Path
import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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
                    ttab.append(codons[codon.upper()])
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

def merge_tables(table_path: str,*args,**kwargs):
    """
    Takes a path w/ wildcards and merges tables with new columns for variable folders / filenames. Positional args are used for column names and, \
    if iterable, regexes to match specific values. Keyword args are passed to pandas read_csv function.

    :params table_path: path to tables with wildcards for variables to expand in new columns
    """
    ppl = table_path.split('*')
    globpaths = glob.glob(table_path)
    relist,colnames = [],[]
    for col in args:
        if type(col) == str:
            relist.append(re.compile('.*'))
            colnames.append(col)
        else:
            relist.append(re.compile(col[1]))
            colnames.append(col[0])
    #print(relist)
    #print(globpaths)
    for k,tp in enumerate(globpaths):
        wcs = []
        regexfail = False
        tpstump = tp
        for i,pp in enumerate(ppl[:-1]):
            tpstump = tpstump[len(pp):]
            wc = tpstump.split(ppl[i+1])[0]
            wcs.append(wc)
            tpstump = tpstump[len(wc):]
            if len(relist) > i:
                if not relist[i].search(wcs[-1]):
                    regexfail = True
                    continue
        if regexfail:
            continue
        tpdf = pd.read_csv(tp,**kwargs)
        for i,col in enumerate(colnames):
            tpdf[col] = wcs[i]
        for i,wc in enumerate(wcs[len(colnames):]):
            tpdf['newcol' + str(i)] = wc
        if k == 0:
            df = tpdf.copy()
        else:
            df = df.append(tpdf)
    return df

def plot_chrms(infile, ref_fai, infmt = 'auto', colors= ['royalblue','grey'],order = 'auto',
               centromere_bed = None, min_colorchange_width = 2000000, 
               legend_names = ['phaseblocks','unphased','assembly gap'],
               min_tig_len = 100000,figsize = (10,5),gap_bed = None,ax = None, feature_bed = None,
               feature_name = None, legend = True, title = "",swath_bed = None,swath_name = None,
               show_alignment_gap = True,chrname_exclude_chars = ['_','EBV'], gapcolor = 'plum',
               refgapcolor = 'white',algngapcolor = 'yellow',swath_and_feature_colors = ['red',
               'gold','grey','purple', 'turqoise','skyblue','darkblue'],
               chrstrip = False, same_scaf_suffix = None, same_scaf_colors = ['skyblue','lightgrey']
               ):
    """
    plots intervals (from paf, bed, or gtf files)
    """
    
    if infmt == 'auto':
        if type(infile) == str:
            if infile.strip('.gz')[-3:] in ['bed','gtf','paf']:
                infmt = infile.strip('.gz')[-3:]
            else:
                infmt = 'paf'
        else:
            infmt = 'paf'
            
    if infile == None:
            intdf = pd.DataFrame(columns = ['sid','sstart','send','qid'])
    elif infmt == 'paf':
            intdf = pd.read_csv(infile,sep="\t",usecols=list(range(12)),header = None)
            intdf.columns = ['qid','qlen','qstart','qend','orient','sid','slen','sstart','send','matches','alen','ascore']
    elif infmt == 'bed':
            intdf = pd.read_csv(infile, sep = '\t', usecols= [0,1,2])
            intdf.columns = ['sid','sstart','send']
            intdf['qid'] = list(intdf.index)
    elif infmt == 'gtf':
            intdf = pd.read_csv(infile, sep = '\t', usecols= [0,1,2,3,4,5,6,7,8],header = None)
            intdf.columns = ['sid','src','kind','sstart','send','score','orient','phase','attributes']
            intdf['qid'] = intdf['attributes'].str.split('gene_id ').str[1].str.split(';'
                ).str[0].str.replace('"','').str.replace(' ','')
    #print(intdf.head())
    intdf['alen'] = intdf['send'] - intdf['sstart']
    intdf = intdf[intdf['alen'] > min_tig_len].sort_values('alen',ascending=False)
    #print(intdf.query('sid == "chr16_MATERNAL"').head())
    if centromere_bed:
        cendf = pd.read_csv(centromere_bed,header=None,sep="\t",index_col = 0)
    if gap_bed:
        gapdf = pd.read_csv(gap_bed,header=None,sep="\t").iloc[:,:3]
        gapdf.columns = ['seq','start','stop']
        gapdf['length'] = gapdf['stop'] - gapdf['start']
    chrdict = {}
    chrlens = {}
    for i,row in intdf.iterrows():
        coords = [row['sstart'],row['send']]
        loc = row['sid']
        if not loc in chrdict:
            chrdict[loc] = IntervalTree()
        ov = False
        for overlap in chrdict[loc].overlap(coords[0],coords[1]):
            if overlap[0] <= coords[0] and overlap[1] >= coords[1]:
                ov = True
        if not ov:
            chrdict[loc][coords[0]:coords[1]] = row['qid']
    excludechars = set(chrname_exclude_chars)
    for line in open(ref_fai):
        fields = line.split('\t')
        if not len(set(fields[0]) & excludechars) > 0:
            chrlens[fields[0]] = int(fields[1])
    if not ax:
        fig,ax = plt.subplots(1,1,figsize=figsize)
        axpassed = False
    else:
        axpassed = True
    xlabels = []
    if order == 'auto':
        enumerator = enumerate(sorted([ [chrlens[k],k] for k in chrlens],reverse=True))
    elif order == 'human':
        enumerator = enumerate([[0, 'chr' + str(k)] for k in list(range(1,23)) + ['M','X','Y'] ])
    else:
        enumerator = enumerate([[0,k] for k in order])
    
    if swath_bed:
        if type(swath_bed) != list:
            swath_bed = [swath_bed]
        swath_list = []
        for bed in swath_bed:
            swath_list.append({})
            for i in open(bed):
                fields = i.split('\t')
                if len(fields) > 1:
                    if not fields[0] in swath_list[-1]:
                        swath_list[-1][fields[0]] = []
                    swath_list[-1][fields[0]].append([int(fields[1]),int(fields[2])])
    else:
        swath_list = []
    
    if feature_bed:
        if type(feature_bed) != list:
            feature_bed = [feature_bed]
        feature_list = []
        for bed in feature_bed:
            feature_list.append({})
            for i in open(bed):
                fields = i.split('\t')
                if len(fields) > 1:
                    if not fields[0] in feature_list[-1]:
                        feature_list[-1][fields[0]] = []
                    feature_list[-1][fields[0]].append(int(fields[1]))
    
    for i,locp in enumerator:
        ax.bar(x=i+1,height = chrlens[locp[1]] / 1000000,color =gapcolor,width = 0.75,edgecolor=colors[0],label = 'gap')
        xlabels.append(locp[1])
        wqid = ""
        color = colors[0]
        astart = 0
        astop = 0
        switch = [True,0]
        if locp[1] in chrdict:
            for k,coords in enumerate(sorted(chrdict[locp[1]])):
                if wqid != coords[2]:
                    pcoords = coords
                    astart = coords[0]
                    switch[0] = True
                else:
                    pcoords = [astart,coords[1]]
                    switch[1] == coords[1]
                if switch[0] and pcoords[0] - switch[1] > min_colorchange_width: # this controls minimum band length I think
                    if same_scaf_suffix: # this expands tig colors out to four, with closer colors for tigs from the same scaffold
                        if wqid.split(same_scaf_suffix)[0] == coords[2].split(same_scaf_suffix)[0]:
                            if color == colors[0]:
                                color = same_scaf_colors[0]
                            elif color == same_scaf_colors[0]:
                                color = colors[0]
                            elif color == colors[1]:
                                color = same_scaf_colors[1]
                            elif color == same_scaf_colors[1]:
                                color = colors[1]
                        else:
                            if color in (color[0],same_scaf_colors[0]):
                                color = colors[1]
                            else:
                                color = colors[0]
                    if color == colors[0]:
                        color = colors[1]
                    elif color == colors[1]:
                        color = colors[0]
                    switch = [False,coords[0]]
                wqid = coords[2]
                astop = coords[1]
                ax.bar(x=i+1,y=pcoords[0]/ 1000000 ,height=  (pcoords[1]-pcoords[0]) / 1000000,color=color,
                       width=0.75,linewidth = 0)

            laststop = 0
            lastid = ""
            for coords in sorted(chrdict[locp[1]]):
                if coords[0] - laststop > 10000 and lastid == coords[2] and show_alignment_gap:
                    ax.bar(x=i+1,y=laststop/ 1000000,height= (coords[0]- laststop) / 1000000,color=(0,0,0,0),
                           edgecolor=algngapcolor,width=0.75,hatch = '//////',linewidth=0)
                elif coords[0] - laststop > 300000 and lastid != coords[2]:
                    ax.bar(x=i+1,y=laststop/ 1000000,height= (coords[0]- laststop) / 1000000,color=(0,0,0,0),
                           edgecolor=colors[0],width=0.75,linewidth=1)
                lastid = coords[2]
                laststop = coords[1]
            ax.bar(x=i+1,height = chrlens[locp[1]] / 1000000,color =(0, 0, 0, 0),
                   width = 0.75,edgecolor=colors[0],linewidth = 1)
        if gap_bed:
            if locp[1] in set(gapdf['seq']):
                chrgapdf = gapdf.query('seq == "' + locp[1] + '"' )
                for index,row in chrgapdf.iterrows():
                    if row['length'] > 3000000:
                        ax.bar(x=i+1,y = row['start']/ 1000000,height = row['length'] / 1000000,color = 'white',
                               width = 0.75,edgecolor=colors[0],label = 'gap')
        if centromere_bed:
            if locp[1] in cendf.index:
                censtart = cendf.at[locp[1],1]
                cenend = cendf.at[locp[1],2]
                ax.bar(x=i+1 - 0.3, y = censtart / 1000000, height = (cenend - censtart) / 1000000,color =(1, 1, 1, 1),
                       width = 0.15,edgecolor=colors[0])
                ax.bar(x=i+1 - 0.35, y = censtart / 1000000, height = (cenend - censtart) / 1000000,color =(1, 1, 1, 1),
                       width = 0.15)
                ax.bar(x=i+1 + 0.3, y = censtart / 1000000, height = (cenend - censtart) / 1000000,color =(1, 1, 1, 1),
                       width = 0.15,edgecolor=colors[0])
                ax.bar(x=i+1 + 0.35, y = censtart / 1000000, height = (cenend - censtart) / 1000000,color =(1, 1, 1, 1),
                       width = 0.15)
        
        if swath_bed:
            for swath_index,swath_dict in enumerate(swath_list):
                if locp[1] in swath_dict:
                    
                    for swath in swath_dict[locp[1]]:
                        ax.bar(x=i+1,y= swath[0]/ 1000000,height= (swath[1] - swath[0]) / 1000000,
                                color= swath_and_feature_colors[swath_index],
                           edgecolor='black',width=0.75,linewidth=0)
        
        if feature_bed:
            for feature_index,feature_dict in enumerate(feature_list):
                if locp[1] in feature_dict:
                    elipses = {}
                    for feature in feature_dict[locp[1]]:
                        elipses[feature] = mpl.patches.Ellipse(
                            (i+1,feature / 1000000),
                            width= 0.5,
                            height = 5,
                            facecolor = swath_and_feature_colors[len(swath_list) + feature_index],
                        )

                        ax.add_patch(elipses[feature])
    
    if legend == True:
        if infile != None:
            pa1 = Patch(facecolor=colors[0], edgecolor=colors[0])
            pa2 = Patch(facecolor=colors[1], edgecolor=colors[0])
            pb1 = Patch(facecolor=refgapcolor, edgecolor=colors[0])
            pb2 = Patch(facecolor=refgapcolor, edgecolor=colors[0])
            pc1= Patch(facecolor=gapcolor,edgecolor=colors[0])
            pc2= Patch(facecolor=gapcolor,edgecolor=colors[0])
            handles_c1 = [pa1, pc1, pb1]
            handles_c2 = [pa2, pc2, pb2]
            labels = legend_names
            if same_scaf_suffix != None:
                handles_c1.append(mpl.patches.Ellipse((-1,-1),4,2,facecolor = colors[0] ))
                handles_c2.append(mpl.patches.Ellipse((-1,-1),4,2,facecolor = same_scaf_colors[0] ))
                handles_c1.append(mpl.patches.Ellipse((-1,-1),4,2,facecolor = colors[1] ))
                handles_c2.append(mpl.patches.Ellipse((-1,-1),4,2,facecolor = same_scaf_colors[1] ))
                labels.append('contigs on same')
                labels.append('    scaffold')
        else:
            handles_c1 = []
            handles_c2 = []
            labels = []
        
        if swath_name:
            if type(swath_name) != list:
                swath_name = [swath_name]
            for findex,fname in enumerate(swath_name):
                handles_c1.append(mpl.patches.Ellipse((-1,-1),4,2,facecolor = swath_and_feature_colors[ findex]))
                handles_c2.append(mpl.patches.Ellipse((-1,-1),4,2,facecolor = swath_and_feature_colors[ findex]))
                labels.append(fname)

        if feature_name:
            if type(feature_name) != list:
                feature_name = [feature_name]
            else:
                for findex,fname in enumerate(feature_name):
                    handles_c1.append(mpl.patches.Ellipse((-1,-1),4,2,
                        facecolor = swath_and_feature_colors[len(swath_list) + findex]))
                    handles_c2.append(mpl.patches.Ellipse((-1,-1),4,2,
                        facecolor = swath_and_feature_colors[len(swath_list) + findex]))
                    labels.append(fname)
        if len(labels) > 0:
            labels = ['' for k in labels] + labels
            handles = handles_c1 + handles_c2
            ax.legend(handles= handles,
                    labels=labels,
                    ncol=2, handletextpad=0.5, handlelength=1.0, columnspacing=-0.5,
                    loc=1, fontsize=10,frameon=False)

    #ax.bar(x=5,y=220,height = 1, color = colors[0])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    ax.set_ylim([-1,np.around(max(chrlens.values())/50000000, decimals=0)*50])
    ax.xaxis.set_ticks(np.arange(1, len(xlabels) + 1))
    ax.tick_params(axis=u'x', which=u'both',length=0)
    ax.axes.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    if chrstrip:
        ax.set_xticklabels([k.replace('chr','') for k in xlabels],rotation = 0)
    else:
        ax.set_xticklabels(xlabels)
        plt.xticks(rotation=90)
    ax.set_title(title,ha='center')

    plt.tight_layout()
    if not axpassed:
        return fig,ax