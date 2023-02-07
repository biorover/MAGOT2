#!/usr/bin/env python
import numpy as np
import pandas as pd
from . import lib
import sys
from pathlib import Path
import ete3
from typing import Iterable, List, Optional, Union

def plot_chrs(*, paf: Optional[List[Path]] = None, ref_fai: Optional[Path] = None,
                centromere_bed: Optional[Path] = None,outprefix: Path = "genome_plots",
                chr_order: Union[List[str]] = 'auto'):
    """
    Generates plots for denovo genome assembly quality evaluation.
    
    :param paf: paf file(s) of genome aligned to a reference (generates dot plots and, if ref_fai \
specified, chromosome plots)
    :param ref_fai: samtools faidx index file(s) for reference (for generating chromosome plots)
    :param centromere_bed: bedfile of centromere locations in reference genome (for making \
chromosome plots pretty)
    :param outprefix: prefix for plot output files
    :param chr_order: order from chromosomes in chromosome plots. Options are "auto" (length sorted), "human" (chr1-22+XY),\
 or list of chromosomes
    """
    if type(chr_order) == list and len(chr_order) == 1:
        chr_order = chr_order[0]
    if paf != False and ref_fai != False:
        for i,paffile in enumerate(paf):
            fig,ax = lib.plot_chrms(paffile, ref_fai, infmt = 'paf',centromere_bed = centromere_bed,
                                    legend_names = ['scaffolds','gap','reference gap'],
                                    same_scaf_suffix = "_contig",order = chr_order)
            fig.savefig(f'{outprefix}.chrmPlot{i}.pdf')

def split_scaffolds(fasta: Path):
    """
    Splits scaffolds in fasta file into contigs (splits at 10+ Ns)

    :param fasta: fasta file input
    """

    ws = []
    for line in open(fasta):
        if line[0] == '>':
            if ws != []:
                seq = "".join(ws)
                seqlines = seq.split('NNNNNNNNNN')
                sl = len(seqlines)
                for i,seqline in enumerate(seqlines):
                if seqline != '':
                    if sl > 1:
                        print(">" + id + "_contig" + str(i))
                    else:
                        print('>' + id)
                    print(seqline.strip('N'))
            id = line[1:].split()[0]
            ws = []
        else:
            ws.append(line.strip())

def sumstats(table:str, *,column: int = 0, cname: str = None, N: bool = False, 
            deciles: bool = False, delim: str = "\t", mode_precision: int = 3):
    """
    Calulate summary stats for a column in an input file

    :param table:  input table with column to summarize
    :param column: int default 0. Column to summarize (zero-based index)
    :param cname: str default None. If specified, selects column number based on field name in header
    :param N: bool default False. If set, returns N0-N100 (in increments of ten) for input data
    :param deciles: bool default False. If set, returns deciles for input data
    :param delim: str default "\\t". Table column delimiter character
    """
    if table == '-':
        table = '/dev/stdin'
    with open(table) as f:
        nlist = []
        for i,line in enumerate(f):
            fields = line.strip().split(delim)
            if i==0 and cname:
                column = fields.index(cname)
            try:
                nlist.append(float(fields[column]))
            except:
                pass
        myarray = np.sort(nlist)[::-1]
        mymax = myarray.max()
        mysum = myarray.sum()
        mymin = myarray.min()
        prec = int(mode_precision - np.log10(mymax))
        sys.stdout.write('Sum: ' + str(mysum) + '\nCount: ' + str(len(myarray)) + '\nMean: ' + str(myarray.mean()) + '\nMedian: ' + 
            str(myarray[int(len(myarray) / 2)]) + "\nMode: " + str(lib.mode(myarray,prec)) + 
            '\nStdev: ' + str(myarray.std()) + '\n')
        if N or deciles:
            toprint = []
            running_total = 0
            tot_len = len(myarray)
            Ns = [100,90,80,70,60,50,40,30,20,10,0]
            Ps = [100,90,80,70,60,50,40,30,20,10,0]
            for i,value in enumerate(myarray):
                running_total += value
                if len(Ns) > 0:
                    if 100 * running_total / mysum >= Ns[-1] and N:
                        sys.stdout.write("N" + str(Ns[-1]) + ": " + str(value) + "\n")
                        Ns.pop()
                if len(Ps) > 0:
                    if 100 * i / tot_len >= Ps[-1] and deciles:
                        toprint.append(str(Ps[-1]) + 'th percentile: ' + str(value) + '\n')
                        Ps.pop()
                if i==len(myarray) - 1:
                    if N and len(Ns) > 0:
                        sys.stdout.write("N" + str(Ns[-1]) + ": " + str(value) + "\n")
                    if deciles and len(Ps) > 0:
                        toprint.append(str(Ps[-1]) + 'th percentile: ' + str(value) + '\n')
            if deciles:
                sys.stdout.write(''.join(toprint))

def gff2fasta(gff: Path, fasta: Path, *, which_transcript: str = 'all', seq_level: str = 'transcript', 
            seq_from: str = 'CDS', seq_type: str = 'nucl', output: Path = '/dev/stdout'):
    """
    Fetches sequences of transcripts from gff (should work with sane gff3 and gtf files) from a fasta

    :param gff: Path. Path to gff file
    :param fasta: Path. Path to fasta file
    :param which_transcript: str default "all". Which transcript to select for each gene. Options are \
"all", "first", "longest", or "best_scoring"
    :param seq_level: str default "transcript". What to return the sequences of - usually "transcript" or "gene" \
(gene only valid if which_transcript = "first", "longest", or "best_scoring"). \
Could also be "CDS" or "exon" or other sub-feature (must then match "seq_from" argument)
    :param seq_from: str default "CDS". Feature type from which to extract sequence (usually "CDS" or "exon")
    :param seq_type: str default "nucl". Whether to output "nucl" (nucleotide), "aa" (translated peptide/amino acid) \
or "lorfaa" (longest orf amino acid) sequence
    """
    annots = lib.read_gff(gff)
    seqs = lib.annot2seqs(annots,fasta,which_transcript = which_transcript, seq_level = seq_level,
                        seq_from = seq_from, seq_type = seq_type)
    with open(output,'w') as f:
        for k,v in seqs.items():
            f.write('>' + k + '\n' + v + '\n')

def fai2tiginfo(fai: Path):
    """
    Annotates tigs with cummulative sequence content and returns info sorted by length

    :params fai: Path. Path to fai file (or any tab-delimeted table with contig name in first column and length in second)
    """
    df = pd.read_csv(fai,header = None, sep = "\t",comment = '#').iloc[:,:2]
    df.columns = ['tig','tiglen']
    df = df.sort_values('tiglen',ascending = False)
    lensum,runsum,tigi = df['tiglen'].sum(),0,0
    for i,row in df.iterrows():
        runsum += row['tiglen']
        tigi += 1
        sys.stdout.write("\t".join([str(k) for k in [row['tig'],row['tiglen'], tigi, int(100 * runsum / lensum), runsum] ]) + '\n')

def gff2gtf(gff: Path, *, output: Path = '/dev/stdout'):
    """
    Takes input gff (sane gff3s or gtfs) and writes out clean CDS gtf for O-G MAGOT / HAPpy-ABCENTH compatibility

    :params gff: Path. Path to input gff file
    """
    gfflines = lib.annot2gtf(lib.read_gff(gff))
    with open(output,'w') as f:
        f.write(gfflines + '\n')

def dipvcf2tripvcf(vcfin: Path,vcfout: Path,*,expand_on: str = "SR"):
    """
    reads a diploid vcf file and writes a triploid vcf file, chosing which allele to duplicate based on the "expand_on" param.

    :param vcfin: Path. Path to diploid input vcf file
    :param vcfout: Path. Path to output triploid vcf file
    :param expand_on: str defualt "SR". What information to use to chose which allele to duplicate. Currently only valid option \
is "SR", which will use the "SR" attribute from the info collumn to add the allele with the highest depth to the alleles field
    """
    lib.dipvcf2tripvcf(vcfin,vcfout,expand_on = expand_on)

def compress_homopolymers(infile: Path,outfile: Path = '/dev/stdout'):
    """
    reads a fasta file and writes a homopolymer compressed fasta file

    :param infile: Path. Input fasta file
    :param ofile: Path. Output fasta file
    """
    lastchar = None
    isdef = False
    with open(infile) as f:
        with open(outfile,'w') as o:
            while True:
                c = f.read(1).upper()
                if not c:
                    break
                if c == '>':
                    isdef = True
                elif c == "\n":
                    isdef = False
                if isdef or not c == lastchar:
                    o.write(c)
                if c != '\n':
                    lastchar = c

def build_clusters(tree_file: Path,seq_file: Path, out_dir: Path, cluster_threshold: float):
    """Reads an ultrametric newick tree and breaks tree into clusters of genes \
    within a specified branch distance of each other, writing each cluster to a separate \
    fasta file in "out_dir".
    
    :param tree_file: Path. Input tree file in newick format (expected to be ultrametric)
    :param out_dir: Path. Directory to which to write output fasta files
    :param seq_file: Path. Fasta file of sequences to split
    :param cluster_threshold: float. Cutoff at which to split clusters
    """
    prot_seq_dict = {}
    for line in open(seq_file):
        if line[0] == ">":
            active_seq = line[1:].strip()
            prot_seq_dict[active_seq] = ""
        else:
            prot_seq_dict[active_seq] += line.strip()
    prottree = ete3.Tree(open(tree_file).read() + ';')
    clusters = []
    dont_append = False
    while prottree.get_distance(prottree.get_leaves()[0]) > cluster_threshold and len(prottree.get_leaves()) > 1:
        children = prottree.get_children()
        child_one_dist = children[0].get_distance(children[0].get_leaves()[0])
        child_two_dist = children[1].get_distance(children[1].get_leaves()[0])
        while child_one_dist > cluster_threshold and child_two_dist > cluster_threshold:
            children = children[0].get_children()
            child_one_dist = children[0].get_distance(children[0].get_leaves()[0])
            child_two_dist = children[1].get_distance(children[1].get_leaves()[0])
        if child_one_dist < cluster_threshold:
            clusters.append(children[0].get_leaf_names())
            try:
                prottree.prune(list(set(prottree.get_leaves()) - set(children[0].get_leaves())), preserve_branch_length = True)
            except:
                if len(prottree.get_leaves()) < 2:
                    dont_append = True
                    break
        if child_two_dist < cluster_threshold:
            clusters.append(children[1].get_leaf_names())
            try:
                prottree.prune(list(set(prottree.get_leaves()) - set(children[1].get_leaves())), preserve_branch_length = True)
            except:
                if len(prottree.get_leaves()) < 2:
                    dont_append = True
                    break
    if not dont_append:
        clusters.append(prottree.get_leaf_names())
    for i,cluster in enumerate(clusters):
        with open(f'{out_dir}/cluster_{i}.fasta','w') as o:
            for gene in cluster:
                seq = prot_seq_dict[gene]
                o.write(f'>{gene}\n{seq}\n')

def simple_cluster(infile: Path,outfile: Path = '/dev/stdout'):
    """
    takes sorted, whitespace seperated pairs of names and outputs \
    clusters of names with any degree of connection.

    :param infile: Path. Input pair file
    :param ofile: Path. Output cluster file   
    """
    clustdict = {}
    with open(infile) as f:
        for line in f:
            p = line.split()
            if p[0] in clustdict:
                if p[1] in clustdict:
                    pass
                else:
                    clustdict[p[1]] = clustdict[p[0]]
            elif p[1] in clustdict:
                clustdict[p[0]] == clustdict[p[1]]
            else:
                clustdict[p[0]] == p[0]
                clustdict[p[1]] == p[0] 
    with open(outfile,'w') as f:
        for node in clusterdict:
            print("function not yet ready")