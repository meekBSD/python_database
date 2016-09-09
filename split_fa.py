#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO

records = SeqIO.parse("ExtractFa.txt", "fasta")

for i in records:
    f_name = str(i.id) + ".fasta"
    f = open(f_name, "w")
    f.write(">"+i.id+"\n")
    f.write(str(i.seq)+"\n")
    f.close()


