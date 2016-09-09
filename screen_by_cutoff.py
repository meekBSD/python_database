#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
from collections import defaultdict
from Bio import SeqIO

USAGE = "python screen_by_cutoff.py *hmm.orfs.txt outputfile"

#USAGE = ''' python ExtractFastaAccordingID.py g3.txt assembly.fa '''


origi_ess_list = ['PF00113', 'PF00121', 'PF00137', 'PF00162', 'PF00164', 
                  'PF00177', 'PF00181', 'PF00189', 'PF00203', 'PF00237',
                  'PF00238', 'PF00252', 'PF00253', 'PF00276', 'PF00281',
                  'PF00297', 'PF00298', 'PF00312', 'PF00318', 'PF00327',
                  'PF00333', 'PF00334', 'PF00344', 'PF00347', 'PF00366',
                  'PF00368', 'PF00380', 'PF00410', 'PF00411', 'PF00416', 
                  'PF00466', 'PF00562', 'PF00572', 'PF00573', 'PF00623', 
                  'PF00670', 'PF00673', 'PF00679', 'PF00687', 'PF00709', 
                  'PF00736', 'PF00749', 'PF00750', 'PF00752', 'PF00827', 
                  'PF00831', 'PF00832', 'PF00833', 'PF00861', 'PF00867', 
                  'PF00900', 'PF00919', 'PF00935', 'PF00958', 'PF01000',
                  'PF01015', 'PF01035', 'PF01090', 'PF01092', 'PF01139',
                  'PF01157', 'PF01172', 'PF01174', 'PF01180', 'PF01191', 
                  'PF01192', 'PF01194', 'PF01198', 'PF01200', 'PF01201', 
                  'PF01246', 'PF01253', 'PF01269', 'PF01280', 'PF01282',
                  'PF01287', 'PF01351', 'PF01496', 'PF01509', 'PF01599',
                  'PF01655', 'PF01667', 'PF01725', 'PF01780', 'PF01798', 
                  'PF01813', 'PF01849', 'PF01864', 'PF01866', 'PF01868',
                  'PF01884', 'PF01896', 'PF01907', 'PF01912', 'PF01920',
                  'PF01922', 'PF01941', 'PF01948', 'PF01949', 'PF01951',
                  'PF01967', 'PF01973', 'PF01981', 'PF01982', 'PF01984', 
                  'PF01990', 'PF01991', 'PF02005', 'PF02006', 'PF02748', 
                  'PF02834', 'PF02938', 'PF02978', 'PF02996', 'PF03439', 
                  'PF03484', 'PF03719', 'PF03764', 'PF03853', 'PF03874', 
                  'PF03876', 'PF03911', 'PF03946', 'PF03947', 'PF03950', 
                  'PF04010', 'PF04019', 'PF04104', 'PF04127', 'PF04135',
                  'PF04406', 'PF04560', 'PF04561', 'PF04563', 'PF04565', 
                  'PF04566', 'PF04567', 'PF04919', 'PF04981', 'PF04983', 
                  'PF04997', 'PF05000', 'PF05221', 'PF05670', 'PF05746',
                  'PF06026', 'PF06093', 'PF06418', 'PF06463', 'PF07541',
                  'PF08068', 'PF08069', 'PF08071', 'PF08608', 'PF08704', 
                  'PF09173', 'PF09239', 'PF09249', 'PF09377', 'PF10559',
                  'PF11987', 'PF04414']
print len(origi_ess_list)

#revis_genes = [i.rstrip().split(".")[0] for i in open("genes.txt", 'r').readlines()]

#print revis_genes

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
print len(rev_genes)

d_cutoff = {}

#original_cutoff = open("162_archaeal_cutoff", "r")
original_cutoff = open("/home/huangjm/database/HMM_essential/162_archaeal_cutoff", "r")
for i in original_cutoff:
    o = i.rstrip().split()
    if o[0] != "PFAM":
        d_cutoff[o[0]] = float(o[1])

original_cutoff.close()

Gene_111 = ['PF00162', 'PF00276', 'PF00281', 'PF00297', 'PF00347', 'PF00366',
            'PF00380', 'PF00410', 'PF00411', 'PF00416', 'PF00466', 'PF00573', 
            'PF00750', 'PF01025', 'PF01795', 'TIGR00001', 'TIGR00002', 'TIGR00009',
            'TIGR00012', 'TIGR00019', 'TIGR00029', 'TIGR00043', 'TIGR00059', 'TIGR00060',
            'TIGR00061', 'TIGR00062', 'TIGR00064', 'TIGR00082', 'TIGR00086', 'TIGR00092', 
            'TIGR00115', 'TIGR00116', 'TIGR00152', 'TIGR00158', 'TIGR00165', 'TIGR00166',
            'TIGR00168', 'TIGR00234', 'TIGR00337', 'TIGR00344', 'TIGR00362', 'TIGR00388', 
            'TIGR00389', 'TIGR00392', 'TIGR00396', 'TIGR00408', 'TIGR00409', 'TIGR00414', 
            'TIGR00418', 'TIGR00420', 'TIGR00422', 'TIGR00435', 'TIGR00436', 'TIGR00442', 
            'TIGR00459', 'TIGR00460', 'TIGR00468', 'TIGR00471', 'TIGR00472', 'TIGR00487',
            'TIGR00496', 'TIGR00575', 'TIGR00631', 'TIGR00663', 'TIGR00775', 'TIGR00810', 
            'TIGR00855', 'TIGR00922', 'TIGR00952', 'TIGR00959', 'TIGR00963', 'TIGR00964', 
            'TIGR00967', 'TIGR00981', 'TIGR01009', 'TIGR01011', 'TIGR01017', 'TIGR01021', 
            'TIGR01024', 'TIGR01029', 'TIGR01030', 'TIGR01031', 'TIGR01032', 'TIGR01044',
            'TIGR01049', 'TIGR01050', 'TIGR01059', 'TIGR01063', 'TIGR01066', 'TIGR01067', 
            'TIGR01071', 'TIGR01079', 'TIGR01164', 'TIGR01169', 'TIGR01171', 'TIGR01391', 
            'TIGR01393', 'TIGR01632', 'TIGR01953', 'TIGR02012', 'TIGR02013', 'TIGR02027',
            'TIGR02191', 'TIGR02350', 'TIGR02386', 'TIGR02387', 'TIGR02397', 'TIGR02432',
            'TIGR02729', 'TIGR03263', 'TIGR03594']


d_rev ={}
for i in Gene_111:
    if i in rev_genes:
        d_rev.setdefault(i, d_cutoff.get(i))

print len(d_rev)

common_gene = d_rev.keys()

for i in d_rev:
    print("{0} => {1}".format(i, d_rev.get(i)))

#print len(d_cutoff)

# check for argument number
if len(sys.argv) != 4:
    sys.exit(USAGE)

in_file = sys.argv[1]
out_file = sys.argv[2]
tmp_contig_ID = sys.argv[3]


filein = open(in_file, 'r')
outputf = open(out_file, "w")

Eleven_gene_fragments = [ "PF00380.16" , "PF00410.16" , "PF00466.17", "PF00281.16", "PF00347.20"  , "PF00573.19" ,"PF00297.19", "PF00366.17" , "PF00411.16", "PF00750.16", "PF00416.19" ]

d_Each_PF = defaultdict(list)


for i in filein:
    if i.startswith("#"):
        continue
    c = i.rstrip().split()
    pfam_ID = c[3].split(".")[0]
    bitScore = float(c[5])
    # if pfam_ID in common_gene and bitScore >= d_rev.get(pfam_ID):
    if pfam_ID in common_gene:
        # print c[0]+"\t"+c[2] +"\t"+c[3] + "\t" + c[5]
        d_Each_PF[c[3]].append(c[0]+":"+c[3]+":"+c[5])
        outputf.write(c[0]+"\t"+c[2] +"\t"+c[3]+"\t"+c[5]+"\n")

filein.close()
outputf.close()

out_contig = open(tmp_contig_ID, "w")

for i in Eleven_gene_fragments:
    for k,v in d_Each_PF.iteritems():
        
        if k == i and len(v) > 1:
            v1 = sorted(v, key=lambda x: float(x.split(":")[2]), reverse=True)
            #print v1[0].split(":")[0] + "\t" + v1[0].split(":")[1]
            out_contig.write("\t".join(v1[0].split(":")[:2]) + "\n")
            break
        elif k == i:
            #print v[0].split(":")[0] + "\t" + v[0].split(":")[1]
            out_contig.write("\t".join(v[0].split(":")[:2]) + "\n" )
            break

out_contig.close()

ID_list = ["##".join(i.rstrip().split("\t")) for i in open(sys.argv[3],'r').readlines()]

fa_name = ".".join(sys.argv[1].split(".")[:-3])

result_fasta = open(fa_name+"_result.fa", "w")

specie_contig_stat = open("Contig_" + fa_name.lstrip("assembly_") + "_stat.txt", "w")

contig_number = 0

fa_dict = SeqIO.to_dict(SeqIO.parse(fa_name+".orfs.faa", "fasta" ) )

for contig_pfam in ID_list:
    seq_ID = contig_pfam.split("##")[0]
    for kc , vc in fa_dict.iteritems():
        if kc == seq_ID:
            contig_number += 1
            if not str(vc.seq).endswith("*"):
                result_fasta.write(">"+contig_pfam+"\n")
                result_fasta.write(str(vc.seq)+"X\n")
                break
            else:

                result_fasta.write(">"+contig_pfam+"\n")
                result_fasta.write(str(vc.seq) + "\n")
                break

result_fasta.close()

#        ID_list[line] = num

records = SeqIO.parse(fa_name+"_result.fa",'fasta')

final_fa = open(fa_name+"_final.fa", "w")

final_fa.write(">"+fa_name+"\n")
for rec in records:
    final_fa.write(str(rec.seq))
final_fa.write("\n")

final_fa.close()

#result_filename = "assembly"+sys.argv[1].split(".")[0]+".faa"
#SeqIO.write(records,result_filename,'fasta')

specie_contig_stat.write("You have extracted "+ str(contig_number) + " contigs for " + fa_name + ".\n")

specie_contig_stat.close()



