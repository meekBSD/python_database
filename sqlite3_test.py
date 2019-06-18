#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sqlite3
import csv

conn = sqlite3.connect('test_gene.db')
cursor = conn.cursor()
cursor.execute('create table IF NOT EXISTS hg18_gene(geneName varchar(20) primary key, chr varchar(20), seq_start integer, seq_end integer)')

filename = 'AceView.ncbi_37.gene2chromosome2coordinates.txt'
# 打开文件
with open(filename) as f:
    # 创建cvs文件读取器
    reader = csv.reader(f, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    # 读取第一行，这行是表头数据
    header_row = next(reader, None)
    # print(header_row)
    # 读取第二行，这行是真正的数据
    #first_row = next(reader)
    for row in reader:
        if row == []:
            break
        cursor.execute('insert into hg18_gene (geneName, chr, seq_start, seq_end) values (\'{0}\', \'{1}\', {2}, {3})'.format(row[0], row[1], row[2], row[3]))

#print(cursor.rowcount)
cursor.close()

conn.commit()
conn.close()

conn = sqlite3.connect('test_gene.db')
# 创建一个cursor:
cursor = conn.cursor()
# 建立索引
cursor.execute('CREATE INDEX index_genes on hg18_gene (chr, seq_start, seq_end);')
#执行查询语句:
#cursor.execute('select *from hg18_gene where chr=? AND seq_start<? AND seq_end>?', ('1',78025000, 78025000,))
#cursor.execute('select *from hg18_gene where chr=? AND seq_start>? AND seq_end<?', ('1',78029000, 78029000,))
cursor.execute('select *from hg18_gene where chr=? AND (seq_start-?)*(seq_end-?)<0', ('1',78029000,78029000,))

#使用featchall获得结果集(list)
values = cursor.fetchall()
print(values) 

# 关闭cursor
cursor.close()
# 关闭conn
conn.close()

