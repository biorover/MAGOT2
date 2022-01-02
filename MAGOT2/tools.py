#!/usr/bin/env python
import numpy as np
import pandas as pd
from . import lib
import sys
from pathlib import Path

def sumstats(table:str, *,column: int = 0, cname: str = None, N: bool = False, 
            deciles: bool = False, delim: str = "\t"):
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
            fields = line.split(delim)
            if i==0 and cname:
                column = fields.index(cname)
            try:
                nlist.append(float(fields[column]))
            except:
                pass
        myarray = np.sort(nlist)[::-1]
        mymax = myarray.max()
        mysum = myarray.sum()
        prec = int(3 - np.log10(mymax))
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