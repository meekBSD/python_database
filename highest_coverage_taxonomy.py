#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
from collections import Counter
from glob import glob
import StringIO
import csv

cov_files = glob("*coverage.csv")

d_cons = defaultdict(list)

with open("test_tmp.txt",'r') as tre:
    for i in tre.readlines():
        line = i.rstrip().split("\t")
        scaffold = line[0]
        orf_order = line[1]
        d_cons[scaffold].append(line[2])

top_Num = len(d_cons)/100

d_ta ={}

for cfile in cov_files:
    with open(cfile, "r") as cov_in:
        data_i = StringIO.StringIO(cov_in.read())
    try:
        reader_cov = csv.reader(data_i, delimiter=",")
        for num, j in enumerate(reader_cov):
            if j[0] in d_cons:
                if j[0] in d_ta:
                    d_ta[j[0]] += float(j[1])
                else:
                    d_ta[j[0]] = float(j[1])
    except Exception, e:
        print str(e)

w_color = sorted(d_ta.iteritems(), key = lambda x: x[1], reverse=True)

scaff_col = {}
count = 0
for s,cov in w_color:
    cnt = Counter(d_cons[s])
    mostcnt = cnt.most_common(1)
    phylum = mostcnt[0][0]
    scaff_col[s] = phylum
    
    count += 1
    if count == top_Num:
        break

Cnt_col = Counter(scaff_col.values())

col_num = 1

for k,v in Cnt_col.most_common():    
    print k+"\t"+str(col_num)
    col_num += 1
    
