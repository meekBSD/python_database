import sys, os

inputdb = sys.argv[1]     # db file name
queryseq= sys.argv[2]     # sequences file name, fasta format

evalue = sys.argv[3]

list1 = []
list2 = []
######### initialize or format a blast database #
cmd1 = 'makeblastdb -in %s -dbtype prot -parse_seqids'%(inputdb)  # makeblast database

print(cmd1)
os.system(cmd1)
cmd2 = 'blastn -query %s -db %s -out blastResult.txt -outfmt 7'%(queryseq, inputdb)  # query similar sequences with db

print(cmd2)
os.system(cmd2)

############
inputfile = open('./blastResult.txt', 'r')
outputfile = open('./out.txt', 'w')

resultlist = inputfile.readlines()

for f1 in resultlist:
    if f1.startswith("A"):    # seq name of sample Name
        f1.rstrip("\n")
        a = f1.split("\t")
        list2.append(a[1])
        if float(a[-2]) < float('%s'%(evalue)):
            list1.append(a[1])
inputfile.close()

print(len(list1))

for i in list1:
    outputfile.writelines(i + "\n")

outputfile.close()
############
cmd3 = 'blastdbcmd -db %s -entry_batch out.txt > out1.fas'%(inputdb)  ## extract sequences

print(cmd3) 
os.system(cmd3)






