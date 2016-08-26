#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from glob import glob
from Bio import SeqIO
from shutil import copy
import os
import sys

cur_path = os.getcwd()

fa_file = glob("*.fasta")

d_stats = {}

for i in fa_file:
    scaff_num = 0
    records = SeqIO.parse(i, "fasta")
    for rec in records:
        scaff_num += 1
    d_stats[i.rstrip(".fasta")] = scaff_num


res_path = os.path.join(cur_path , "essResult")
if not os.path.isdir(res_path):
    os.makedirs(res_path)
ass_hmm_res = os.popen('find ./ -name "assembly_*hmm.orfs.txt"').readlines()

for i in ass_hmm_res:
    copy(i.rstrip(), res_path)


os.chdir(res_path)

print len(os.listdir("."))

rev_genes = ['PF00252', 'PF00318', 'PF00327', 'PF00333', 'PF00380',
             'PF00416', 'PF00466', 'PF00572', 'PF00935', 'PF01000', 
             'PF01157', 'PF01172', 'PF01200', 'PF01280', 'PF01655',
             'PF01920', 'PF01984', 'PF03439', 'PF03719', 'PF04406', 
             'PF09239', 'PF09377', 'PF00164', 'PF00177', 'PF00189',
             'PF00203', 'PF00237', 'PF00281', 'PF00344', 'PF00347', 
             'PF00410', 'PF00411', 'PF00562', 'PF00623', 'PF00673',
             'PF00679', 'PF00687', 'PF00749', 'PF00750', 'PF00867', 
             'PF01090', 'PF01139', 'PF01191', 'PF01194', 'PF01351',
             'PF02834', 'PF02996', 'PF03484', 'PF03764', 'PF03876', 
             'PF03947', 'PF04560', 'PF04561', 'PF04563', 'PF04565', 
             'PF04566', 'PF04567', 'PF04919', 'PF04983', 'PF04997', 
             'PF05000', 'PF05746', 'PF07541', 'PF09249', 'PF00238', 
             'PF00297', 'PF00298', 'PF00366', 'PF00573', 'PF00752', 
             'PF00827', 'PF01015', 'PF01192', 'PF01201', 'PF01253', 
             'PF01780', 'PF01866', 'PF01951', 'PF01981', 'PF02005', 
             'PF03946', 'PF03950', 'PF06093', 'PF09173', 'PF11987']

d_rev = {}

original_cutoff = open("/home/huangjm/database/HMM_essential/162_archaeal_cutoff", "r")
for i in original_cutoff:
    o = i.rstrip().split()
    if o[0] in rev_genes:
        d_rev[o[0]] = float(o[1])

original_cutoff.close()

#for i in d_rev:
#    print i, d_rev.get(i, "None")

hmm_orfs = glob("assembly_*hmm.orfs.txt")


for f in hmm_orfs:
    with open(f, "r") as hmmSearchRes, open("Cutoff_"+f, "w") as out:
        for i in hmmSearchRes:
            if i.startswith("#"):
                continue
            c = i.rstrip().split()
            pfam_ID = c[3].split(".")[0]
            bitScore = float(c[5])
            if pfam_ID in rev_genes and bitScore >= d_rev.get(pfam_ID):
            #if pfam_ID in common_gene:
                # print c[0]+"\t"+c[2] +"\t"+c[3] + "\t" + c[5]
                #d_Each_PF[c[3]].append(c[0]+":"+c[3]+":"+c[5])
                out.write(c[0]+"\t"+c[2] +"\t"+c[3]+"\t"+c[5]+"\n")

cutoff_Res = glob("Cutoff_*")
essen_res = open("essential_result.txt", "w")

for file_cut in cutoff_Res:
    with open(file_cut, "r") as fc:
        essen_res.write(file_cut+"\t")
        for f_line in fc:
            split_line = f_line.rstrip().split("\t")
            essen_res.write("\t"+split_line[2])
        essen_res.write("\n")

essen_res.close()

input_file = open("essential_result.txt", "r")

all_ess = []

pfam_numStat = open("essential_stat.txt",'w')

for num, i in enumerate(input_file.readlines()):
    l = i.rstrip().split("\t")
    pfam_numStat.write(l[0] + "\t" + str(len(l[2:])) +"\n")
    tot_ess = [ x.split(".")[0] for x in l[2:]]
    singleton_num = 0
    uniq_num = len(set(tot_ess))
    for j in rev_genes:
        if tot_ess.count(j) == 1:
            singleton_num += 1
    all_ess.append(l[0].rstrip(".hmm.orfs.txt")[16:]+"\t"+str(len(tot_ess))+"\t"+str(uniq_num)+"\t"+str(singleton_num))

    if num == 0:
        pfams = set(l[1:])
        #print pfams
        #print len(pfams)
    elif num > 0:
        pfams = set(l[1:]) & pfams
        #print len(pfams)

pfam_numStat.close()

print len(pfams)
print len(all_ess)

all_pfam = sorted(all_ess, key = lambda x: int(x.split("\t")[1]), reverse=True)
specie_pfam_num = open("pfams_stat.txt", "w")

specie_pfam_num.write("species\ttot_ess\tuniq_ess\tsingleton\tscaff_num\n")

for i in  all_pfam:
    each_spe = i.split("\t")[0]
    specie_pfam_num.write(i+"\t"+str(d_stats.get(each_spe, "_"))+"\n")

specie_pfam_num.close()
#for i in pfams:
#    print i




filein = open("essential_result.txt", 'r')
outputf = open("res2.txt", "w")

outputf.write("genes_PFAM\t")
for num, pfam in enumerate(rev_genes):
    outputf.write("%s\t" % pfam)
outputf.write("\n")

for i in filein:
    c = [ n.split(".")[0] for n in i.rstrip().split("\t") ]
    outputf.write(c[0].rstrip(".orfs.hmm_stat.txt")+"\t")
    for j in rev_genes:
        gene_count = c.count(j)
        outputf.write(str(gene_count)+"\t")
        #if j in c:
        #    outputf.write("5\t")
        #else:
        #    outputf.write("1\t")
    outputf.write("\n")
filein.close()
outputf.close()


