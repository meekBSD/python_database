#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import with_statement
import sys
import re

nodescon = open("nodes.dmp","r")

classifi_t = {}

#phylist = []
#pattern = re.compile(r"\t(.*)phylum\t")

for i in nodescon:
    c = i.rstrip().split("\t")
    classifi_t[c[0]] = c[2]

    #if pattern.search(i):
    #    phylist.append(i.split("\t")[0])

del classifi_t["1"]

nodescon.close()

tax_dict = {}
input_taname = open("names.dmp", "r")

for each_line in input_taname:
    ac = each_line.rstrip().split("\t")
    if each_line.find("scientific name") != -1:
        tax_dict[ac[0]] = ac[2]

input_taname.close()


taxonomyID = open("taxid", 'r')

all_tax = set(classifi_t.keys())


count = 0
test_out = open("test_result_name.txt" , "w")

for i in taxonomyID:

    ta = i.rstrip()
    test_out.write(ta + "\t")
    test_out.write(tax_dict.get(ta, "-"))
    while True:
        if ta in all_tax:
            ta = classifi_t[ta] 
            #test_out.write(ta + "\t")
            test_out.write(";"+ tax_dict.get(ta,"-"))
        else:
            break
    test_out.write("\n")    
 
taxonomyID.close()
test_out.close()


phylumSet = set()


with open("nodes.dmp",'r') as n:
    for i in n.readlines():
        idp = i.rstrip().split("\t")[0]
        ThirdColumn = i.rstrip().split("\t")[4]
        if ThirdColumn.find("phylum") != -1:
            phylumSet.add(tax_dict.get(idp, "-"))


outfile = open("result_taname_2.txt",'w')

all_ta = open("test_result_name.txt", 'r')

for i in all_ta:
    tax_id = i.rstrip().split("\t")[0]
    name = i.rstrip().split("\t")[1]
    Notfound = True
    for i in phylumSet:
        if name.find(i) != -1:
            outfile.write(tax_id + "\t" + i + "\t" + name + "\n")
            Notfound = False
            break
    if Notfound:
        if name.find("virus") != -1 or name.find("Virus") != -1:
            outfile.write(tax_id + "\t" + "Virus" + "\t" + name +"\n")
        elif name.find("Eukaryota") != -1:
            outfile.write(tax_id + "\t" + "Eukaryota" +"\t" + name +"\n")
        elif name.find("division WS6") != -1:
            outfile.write(tax_id + "\t" + "WS6" +"\t"+ name +"\n")
        elif name.find("Parcubacteria group") != -1:
            outfile.write(tax_id + "\t" + "Parcubacteria" + "\t" + name +"\n")
        elif name.find("division Kazan-3B-28") != -1:
            outfile.write(tax_id + "\t" + "Kazan" +"\t" + name +"\n")
        elif name.find("Bacteria (miscellaneous)") != -1:
            outfile.write(tax_id + "\t" + "miscellaneous unclassified" +"\t" + name +"\n")
        elif name.find("Candidatus Dependentiae")!= -1:
            outfile.write(tax_id + "\t" + "Dependentiae" + "\t" + name +"\n")
        elif name.find("division BRC1") != -1:
            outfile.write(tax_id + "\t" + "BRC1" + "\t" + name + "\n")
        elif name.find("Microgenomates group") != -1:
            outfile.write(tax_id + "\t" + "Microgenomates" + "\t" + name + "\n")
        elif name.find("Haloplasmatales") != -1:
            outfile.write(tax_id + "\t" + "Haloplasmatales" + "\t" + name + "\n")
        elif name.find("division TA06") != -1:
            outfile.write(tax_id + "\t" + "TA06" + "\t" + name + "\n")
        elif name.find("division WOR-1") != -1:
            outfile.write(tax_id + "\t" + "WOR-1" + "\t" + name + "\n")
        elif name.find("division SR1") != -1:
            outfile.write(tax_id + "\t" + "SR1" + "\t" + name + "\n")
        elif name.find("Thermobaculum") != -1:
            outfile.write(tax_id + "\t" + "Thermobaculum" + "\t" + name + "\n")
        elif name.find("samples;Bacteria") != -1 or name.find("Bacteria;cellular") != -1:
            outfile.write(tax_id + "\t" + "unclassified Bacteria" + "\t" + name +"\n")	
        elif name.find("prokaryotic environmental") != -1:
            outfile.write(tax_id + "\t" + "Prokaryotic" + "\t" + name +"\n")	
        elif name.find("Archaea") != -1:
            outfile.write(tax_id + "\t" + "Archaea" + "\t" + name +"\n")	
        else:
            outfile.write(tax_id + "\tUnclassified" + "\t" + name + "\n")

all_ta.close()
outfile.close()

        
d={}

a_in = open("result_taname_2.txt", "r")

d_pro = {}

d_Euk = {}
for line in a_in:
    csplit = line.rstrip().split("\t")
    if csplit[2].find("Eukaryota") != -1:
        d_Euk[csplit[0]] = csplit[1]
    else:
        d_pro[csplit[0]]=csplit[1]

a_in.close()

pro = open("protokary_tax.txt", "w")
euk = open("eukaryota_tax.txt", "w")

proteinID_ta = open("protID2taxID.txt", "r")

for line in proteinID_ta:
    s = line.rstrip().split("\t")
    if s[1] in d_pro:
        pro.write(s[0]+"\t"+d_pro[s[1]]+"\n")
       # print s[0]+"\t"+d[s[1]]
    elif s[1] in d_Euk:
        euk.write(s[0]+"\t"+d_Euk[s[1]]+"\n")

    else:
        pro.write(s[0]+"\tUnclassified\n")
        #print s[0]+"\t"+"Unclassified"

pro.close()
euk.close()

