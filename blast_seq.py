#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
from glob import glob
from Bio.Blast import NCBIWWW

def parse_fasta(filename):
#    f = open(filename, 'r')
    sequences = {}
    for line in filename:
        if line.startswith(">"):
            name = line.rstrip("\n")
            sequences[name] = ""
        else:
            sequences[name] = sequences[name] + line.rstrip("\n")

    return [(i[0]+"\n"+i[1]) for i in sequences.items()]

seq_file = glob("*.seq")

sav_fil = open("result.xml",'w')

for i in seq_file:
    with open(i, "r") as in_put:
        seqL = parse_fasta(in_put.readlines())
        sort_seqs = sorted(seqL, key=lambda x:len(x.split("\n")[1]), reverse=True)
        a = sort_seqs[0]
        result_handle = NCBIWWW.qblast("blastn", "nt",a, entrez_query='all [filter] NOT(environmental samples[organism] OR metagenomes[organism])', format_type="Text", hitlist_size = 25)
        sav_fil.write(result_handle.read())
        #for i in sorted(seqL, key=lambda x:len(x.split("\n")[1]), reverse=True):
        #    print i

sav_fil.close()
