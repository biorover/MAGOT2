#!/usr/bin/env python
import sys
import numpy as np
from intervaltree import IntervalTree,Interval

def mode(a,prec = 2):
    vals,counts =  np.unique(np.around(a,prec),return_counts=True)
    return vals[np.argmax(counts)]

def translate(x,padNs = False,frames = 1):
    ttabs = []
    codons = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
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
            if len(codon) == 3 or padNs:
                if codon.upper() in codons:
                    ttab.append(codons[codon])
                else:
                    ttab.append('N')
        ttabs.append("".join(ttab))
    if rstring:
        return ttabs[0]
    else:
        return ttabs

def revcomp(x):
    comp = {'A':'T','T':'A','C':'G','G':'C',
            'a':'t','t':'a','c':'g','g':'c'}
    rctab = []
    for char in x[::-1]:
        rtab.append(comp[char])
    return ''.join(rctab)

def orfs(x,best = False): #todo: #1 add info table
    orfs = []
    for frames in translate(x,frames = 'all'):
        orfs += frames.split('*')
    if best:
        lens = [len(k) for k in orfs]
        return orfs[np.argmax(lens)]
    else:
        return orfs



