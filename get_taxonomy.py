#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
import argparse

Help = """USAGE: python dereplicate_fasta.py New_phylum.faa NR_pre.faa NR_file.faa

       You shall provide two fasta files and one output file name."""


parser = argparse.ArgumentParser(description = Help)

parser.add_argument("in_newphyla")
parser.add_argument("NR_pro")
parser.add_argument("out_file")

args = parser.parse_args()

records_new = SeqIO.to_dict(SeqIO.parse(args.in_newphyla, "fasta"))

Pretreated_fa = SeqIO.parse(args.NR_pro, "fasta" )

result = open(args.out_file, "w")

for rec in Pretreated_fa:
    flag = True
    for k, i in records_new.iteritems():
        if str(i.seq) == str(rec.seq):
            flag =False
            break
    if flag :
        SeqIO.write(rec, result, "fasta")

result.close()



d = {}

file_ta = open("/home/qiaozy/Projects/protokary_tax.txt", "r")

for i in file_ta:
    co = i.rstrip().split("\t")
    if co[1].find("Candidatus") != -1:
        taxNam = "#".join(co[1].replace("Candidatus ","").split())
    elif co[1].find("candidate division") != -1:
        taxNam = "#".join(co[1].replace("candidate division ","").split())
    else:
        taxNam = "#".join(co[1].split())
    d[co[0]]= taxNam

file_ta.close()

result_fa = open("Output_prokaryotic_NR.faa", "w")


for Seq in SeqIO.parse(args.out_file, "fasta"):
    a = Seq.id.split("|")[3]
    if a in d:
        Seq.id = d[a]+"_" + Seq.id
        SeqIO.write(Seq, result_fa, "fasta")
    else:
        Seq.id = "Unclassified_"+Seq.id
        SeqIO.write(Seq, result_fa, "fasta")

result_fa.close()

