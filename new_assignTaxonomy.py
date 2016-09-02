#!/usr/bin/env python
# -*- coding: UTF-8 -*-

################################################################################################
##      Input File : assembly.orfs.hmm.NewPhylum.blast.txt                               ##
##      Three Output Files: test_tmp.txt                                                      ##
##                          test_blast_tab.txt                                                ##
##                          assembly.tax.consensus_"+ num + ".txt                             ##
##                                                                                            ##
################################################################################################

from __future__ import with_statement
from collections import defaultdict
from collections import Counter
import os
import sys
import argparse

Help = """USAGE: python % prog.py % blast_result_file . You shall provide blast result file name."""


parser = argparse.ArgumentParser(description = Help)

parser.add_argument("blastfile")
args = parser.parse_args()


d=defaultdict(list)

NewPhyla = set(["Firmicutes", "Euryarchaeota", "Proteobacteria",
             "delta/epsilon#subdivisions", "Chloroflexi" , "Bacteroidetes", 
             "Planctomycetes" , "Actinobacteria" , "Spirochaetes" , 
             "Thaumarchaeota" ,  "Parcubacteria" , "Crenarchaeota" , 
             "Cyanobacteria" , "Archaea" ,"Ascomycota",
             "unclassified#Bacteria" , "NC10", "OP1", 
             "SR1", "Zixibacteria",  "WWE1",     
             "TM7", "Atribacteria", "Kryptonia",
             "TA06","Omnitrophica","Aerophobetes",
             "Spirochaetes", "Poribacteria","WOR-1",
             "Hyd24-12","Latescibacteria","KD3-62",
             "Cloacamonas","BRC1", "WOR-3",
             "OP11","OD1","WWE3",
             "Peregrinibacteria", "Saccharibacteria","OP3",
             "Microgenomates", "Korarchaeota","Geoarchaeota",
             "Bathyarchaeota","Diapherotrites","Woesearchaeota",
             "Parvarchaeota","Aenigmarchaeota","Micrarchaeota",
             "Nanoarchaeota","Nanohaloarchaeota", "Unclassified"
    ])

## CD12 is "Aerophobetes",  WOR-2 is "Omnitrophica"

bestHit_E = {}

with open(args.blastfile, 'r') as tan:
    for i in tan.readlines():
        allcol = i.rstrip().split("\t")
        
        Evalue = float(allcol[4])                # get the e-value of each hit
        if d[allcol[0]] == []:
            bestHit_E[allcol[0]] = Evalue
            d[allcol[0]].append(allcol[2].split("_")[0])
        elif Evalue <= 10 * bestHit_E[allcol[0]]:
            d[allcol[0]].append(allcol[2].split("_")[0])

        #d[allcol[0]].append(allcol[2].split("_")[0])

print os.getcwd()
cur_path = os.getcwd()

res_path = os.path.join(cur_path , "assign_result")
if not os.path.isdir(res_path):
    os.makedirs(res_path)

os.chdir(res_path)
print ("Your work will be stored in the %s" % res_path)


out_tmp = open("test_tmp.txt",'w')
for k,v in d.iteritems():
    kk =k.split("_")
    target_n = float(len(v))
    cnt = Counter(v).most_common(2)
    if len(cnt) == 2 and cnt[0][0] in NewPhyla:
        out_tmp.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(kk[0], kk[1],v[0], cnt[0][0], cnt[0][1]/target_n, cnt[1][0], cnt[1][1]/target_n, target_n))
        #print("{0}\t{1}\t{2}\t{3}\t{4}".format(k, cnt[0][0],cnt[0][1], cnt[1][0], cnt[1][1]))
    elif len(cnt) == 1 and cnt[0][0] in NewPhyla:
        out_tmp.write("{0}\t{1}\t{2}\t{3}\t{4:.2f}\t{5}\t{6}\t{7}\n".format(kk[0], kk[1],v[0], cnt[0][0],cnt[0][1]/target_n, "-", "-", target_n))
        #print("{0}\t{1}\t{2}\t{3}\t{4}".format(k, cnt[0][0],cnt[0][1], "-", "-"))

out_tmp.close()

blast_tab = open("test_blast_tab.txt",'w')
d_cons = defaultdict(list)

with open("test_tmp.txt",'r') as tre:
    for i in tre.readlines():
        line = i.rstrip().split("\t")
        scaffold = line[0]
        orf_order = line[1]
        d_cons[scaffold].append(line[2])
        blast_tab.write("{0}\t{1}\t{2}\n".format(scaffold, orf_order, line[2]))

assmbly_results = [i for i in os.listdir(".") if i.startswith("assembly.tax.consensus")]

outfile = "assembly.tax.consensus_"+str(len(assmbly_results)+1) + ".txt"

scaff_tax = {}

most_tax = []
for k,v in d_cons.iteritems():    
    cnt = Counter(v)
    mostcnt = cnt.most_common(1)
    phylum = mostcnt[0][0]
    most_tax.append(mostcnt[0][0])
    all_phylum = ";".join([i[0] for i in cnt.items()])
    scaff_tax[k] = "\t".join([phylum, all_phylum])

d_taxonomy = {}
tax_assign = Counter(most_tax)

tax_color_file = open("temp_color_assign.txt", "w")
tax_lab = 0
tax_not_in = 53
for ki, vi in tax_assign.most_common():
    if ki in NewPhyla:
        tax_lab += 1
        d_taxonomy[ki] = tax_lab
    else:
        tax_not_in += 1
        d_taxonomy[ki] = tax_not_in

tax_color = sorted(d_taxonomy.iteritems(), key = lambda x: int(x[1]))

tax_color_file = open("temp_color_assign.txt", "w")

for k,v in tax_color:
    tax_color_file.write("{0}\t{1}\n".format(k, v))


tax_color_file.close()

out_consensus = open(outfile, 'w')

out_consensus.write("scaffold\tphylum\ttaxcolor\tall.assignments\n")

for kf, vf in scaff_tax.iteritems():
    scaff_x = vf.split("\t")
    phylum = scaff_x[0]
    out_consensus.write("{0}\t{1}\t{2}\t{3}\n".format(kf, phylum, d_taxonomy.get(phylum, "118"), scaff_x[1]))

out_consensus.close()
#print num


