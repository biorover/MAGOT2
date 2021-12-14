#!/usr/bin/env python
import numpy as np
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
                        toprint.append(str(Ps[-2]) + 'th percentile: ' + str(value) + '\n')
                        Ps.pop()
            if deciles:
                sys.stdout.write(''.join(toprint))

def gff2fasta(gff: Path, fasta: Path, *, name_from: str = 'transcript', which_transcript: str = 'all', 
            seq_from: str = 'CDS', seq_type: str = 'nucl', out_file: Path = '/dev/stdout'):
    """
    Fetches sequences of transcripts from gff (should work with sane gff3 and gtf files) from a fasta

    :param gff: Path. Path to gff file
    :param fasta: Path. Path to fasta file
    :param name_from: str default "transcript". Whether to name sequences after "transcript" or "gene" identifier \
(gene identifier only valid if which_transcript = "first", "longest", or "best_scoring")
    :param which_transcript: str default "all". Which transcript to select for each gene. Options are \
"all", "first", "longest", or "best_scoring"
    :param seq_from: str default "CDS". Feature type from which to extract sequence (usually "CDS" or "exon")
    :param seq_type: str default "nucl". Whether to output "nucl" (nucleotide), "aa" (translated peptide/amino acid) \
or "lorfaa" (longest orf amino acid) sequence
    """
    annots = lib.read_gff(gff)
    seqs = lib.annot2seqs(annots,fasta,which_transcript = which_transcript, name_from = name_from,
                        seq_from = seq_from, seq_type = seq_type)
    out = open(out_file,'w')
    for k,v in seqs.items():
        out.write('>' + k + '\n' + v + '\n')