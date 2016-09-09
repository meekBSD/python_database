#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from glob import glob
from StringIO import StringIO
from Bio import AlignIO
import sys
import time

Eleven_gene_fragments = [ "PF00380.16" , "PF00410.16" , "PF00466.17", "PF00281.16", "PF00347.20"  , "PF00573.19" ,"PF00297.19", "PF00366.17" , "PF00411.16", "PF00750.16", "PF00416.19" ]

all_spe = [ 
         'AWOB01.1_Geoarchaeota' ,
         'Acaryochloris_marina_MBIC11017' ,
         'Acetobacter_pasteurianus_386B' ,
         'Acetobacterium_woodii_DSM_1030' ,
         'Acetohalobium_arabaticum_DSM_5501' ,
         'Achromobacter_xylosoxidans_A8' ,
         'Acidaminococcus_fermentans_DSM_20731' ,
         'Acidianus_hospitalis_W1' ,
         'Acidilobus_saccharovorans_345_15' ,
         'Acidimicrobidae_bacterium_YM16_304' ,
         'Acidimicrobium_ferrooxidans_DSM_10331' ,
         'Acidiphilium_cryptum_JF_5' ,
         'Acidobacterium_capsulatum_ATCC_51196' ,
         'Acidothermus_cellulolyticus_11B' ,
         'Aciduliprofundum_boonei_T469' ,
         'Aequorivita_sublithincola_DSM_14238' ,
         'Test_Spe_1' ,                                     # test specie 1
         'Aeropyrum_camini_SY1' ,
         'Akkermansia_muciniphila_ATCC_BAA_835' ,
         'Alistipes_finegoldii_DSM_17242' ,
         'Aminobacterium_colombiense_DSM_12261' ,
         'Anabaena_cylindrica_PCC_7122' ,
         'Anaerobaculum_mobile_DSM_13181' ,
         'Anaerolinea_thermophila_UNI_1' ,
         'Aquifex_aeolicus_VF5' ,
         'Archaeoglobus_fulgidus_DSM_4304' ,
         'Arthrospira_platensis_NIES_39' ,
         'Bacteroides_CF50' ,
         'Borrelia_afzelii_HLJ01' ,
         'Brachyspira_hyodysenteriae_WA1' ,
         'Caldiarchaeum_subterraneum' ,
         'Caldilinea_aerophila_DSM_14535___NBRC_104270' ,
         'Caldisericum_exile_AZM16c01' ,
         'Calditerrivibrio_nitroreducens_DSM_19672' ,
         'Chlamydia_muridarum_Nigg' ,
         'Chlamydophila_abortus_S26_3' ,
         'Chlorobaculum_parvum_NCIB_8327' ,
         'Chlorobium_chlorochromatii_CaD3' ,
         'Chloroflexus_aggregans_DSM_9485' ,
         'Chloroherpeton_thalassium_ATCC_35110' ,
         'Chthonomonas_calidirosea_T49' ,
         'Coraliomargarita_akajimensis_DSM_45221' ,
         'Deferribacter_desulfuricans_SSM1' ,
         'Deinococcus_deserti_VCD115' ,
         'Denitrovibrio_acetiphilus_DSM_12809' ,
         'Desulfurispirillum_indicum_S5' ,
         'Desulfurobacterium_thermolithotrophum_DSM_11699' ,
         'Dictyoglomus_thermophilum_H_6_12' ,
         'Elusimicrobium_minutum_Pei191' ,
         'Euryarchaeota_3' ,
         'Fervidobacterium_nodosum_Rt17_B1' ,
         'Fibrobacter_succinogenes_S85' ,
         'Firmicutes_cd12_1' ,
         'Fusobacterium_3_1_36A2' ,
         'GW2011_AR_AR10_43_0_Diapherotrites' ,
         'GW2011_AR_unknown_32_20_Woesearchaeota' ,
         'GWA2_AR18_30_20_Woesearchaeota' ,
         'GWA2_AR5_49_1_Aenigmarchaeota' ,
         'Gemmatimonas_aurantiaca_T_27' ,
         'Haloarcula_sp_CBA1115' ,
         'Halyomorpha_halys_symbiont' ,
         'Hydrogenobacter_thermophilus_TK_6' ,
         'Ignavibacterium_album_JCM_16511' ,
         'Ilyobacter_polytropus_DSM_2926' ,
         'Isosphaera_pallida_ATCC_43644' ,
         'Korarchaeum_cryptofilum_OPF8' ,
         'Koribacter_versatilis_Ellin345' ,
         'Kosmotoga_olearia_TBF_19_5_1' ,
         'LFWW01.1_Bathyarchaeota' ,
         'LKMY01.1' ,
         'LUCB01.1_Bathyarchaeota' ,
         'LUCE01.1_Bathyarchaeota' ,
         'Leptospirillum_ferriphilum_ML_04' ,
         'Leptotrichia_buccalis_C_1013_b' ,
         'Lokiarchaeum_sp_GC14_75' ,
         'Marinithermus_hydrothermalis_DSM_14884' ,
         'Marinitoga_piezophila_KA3' ,
         'Meiothermus_ruber_DSM_1279' ,
         'Melioribacter_roseus_P3M' ,
         'Metallosphaera_cuprina_Ar-4' ,
         'Methanocaldococcus_infernus_ME' ,
         'Methanothermobacter_thermautotrophicus_str_Delta_H' ,
         'Methylacidiphilum_infernorum_V4' ,
         'Methylomirabilis_oxyfera' ,
         'Micrarchaeum_acidiphilum_ARMAN-2' ,
         'Nanoarchaeota_archaeon_7A' ,
         'Nanoarchaeote_Nst1_Nst1_C1' ,
         'Nanoarchaeum_equitans_Kin4_M' ,
         'Nanosalina_sp_J07AB43' ,
         'Nanosalinarum_sp_J07AB56' ,
         'New_Euryarchaeota_1' ,
         'New_Firmicutes_1' ,
         'New_unclassified_archaea' ,
         'Nitrosopelagicus_brevis_strain_CN25' ,
         'Nitrosopumilus_AR2' ,
         'Nitrososphaera_gargensis_Ga9_2' ,
         'Nitrososphaera_viennensis_EN76' ,
         'Nitrosotenuis_cloacae_strain_SAT1' ,
         'Nitrospira_defluvii' ,
         'Parvarchaeum_acidophilus_ARMAN-5' ,
         'Phycisphaera_mikurensis_NBRC_102666' ,
         'Pirellula_staleyi_DSM_6068' ,
         'Protochlamydia_amoebophila_UWE25' ,
         'Pyrobaculum_sp_WP30' ,
         'SR1_bacterium_RAAC1_SR1_1' ,
         'Saccharibacteria_bacterium_RAAC3_TM7_1' ,
         'Saccharobacterium_alaburgensis' ,
         'Solibacter_usitatus_Ellin6076' ,
         'Sphaerochaeta_pleomorpha_Grapes' ,
         'Sulfolobus_islandicus' ,
         'Synergistetes_bacterium_SGP1' ,
         'Test_Spe_2' ,                                   # test specie 2
         'TA06_bacterium_DG_24_GCA' ,
         'TA06_bacterium_DG_26_GCA' ,
         'TA06_bacterium_DG_78_GCA' ,
         'TA06_bacterium_SM1_40_GCA' ,
         'TA06_bacterium_SM23_40_GCA' ,
         'Thermodesulfatator_indicus_DSM_15286' ,
         'Thermodesulfobacterium_OPB45' ,
         'Thermodesulfovibrio_yellowstonii_DSM_11347' ,
         'Thermofilum_sp_1910b' ,
         'WWE3_bacterium_RAAC2_WWE3_1' ,
         'archaeon_Mx1201' ,
         'halophilic_archaeon_DL31'             
          ]

all_fas = glob("*_result.fa")

spe_contain_380 = []
spe_contain_410 = []
spe_contain_466 = []
spe_contain_281 = []
spe_contain_347 = []
spe_contain_573 = []
spe_contain_297 = []
spe_contain_366 = []
spe_contain_411 = []
spe_contain_750 = []
spe_contain_416 = []

PF00380_16 = open("conserved_gene_PF00380.fa", "w")
PF00410_16 = open("conserved_gene_PF00410.fa", "w")
PF00466_17 = open("conserved_gene_PF00466.fa", "w")
PF00281_16 = open("conserved_gene_PF00281.fa", "w")
PF00347_20 = open("conserved_gene_PF00347.fa", "w")
PF00573_19 = open("conserved_gene_PF00573.fa", "w")
PF00297_19 = open("conserved_gene_PF00297.fa", "w")
PF00366_17 = open("conserved_gene_PF00366.fa", "w")
PF00411_16 = open("conserved_gene_PF00411.fa", "w")
PF00750_16 = open("conserved_gene_PF00750.fa", "w")
PF00416_19 = open("conserved_gene_PF00416.fa", "w")

for i in all_fas:
    for j in all_spe:
        if i.find(j) != -1:
            species_name = j
    records = SeqIO.parse(i, "fasta")
    for rec in records:
        if rec.id.find("PF00380.16") != -1:
            spe_contain_380.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00380_16, "fasta")
        elif rec.id.find("PF00410.16") != -1:
            spe_contain_410.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00410_16, "fasta")
        elif rec.id.find("PF00466.17") != -1:
            spe_contain_466.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00466_17, "fasta")
        elif rec.id.find("PF00281.16") != -1:
            spe_contain_281.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00281_16, "fasta")
        elif rec.id.find("PF00347.20") != -1:
            spe_contain_347.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00347_20, "fasta")
        elif rec.id.find("PF00573.19") != -1:
            spe_contain_573.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00573_19, "fasta")
        elif rec.id.find("PF00297.19") != -1:
            spe_contain_297.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00297_19, "fasta")
        elif rec.id.find("PF00366.17") != -1:
            spe_contain_366.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00366_17, "fasta")
        elif rec.id.find("PF00411.16") != -1:
            spe_contain_411.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00411_16, "fasta")
        elif rec.id.find("PF00750.16") != -1:
            spe_contain_750.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00750_16, "fasta")
        elif rec.id.find("PF00416.19") != -1:
            spe_contain_416.append(species_name)
            rec.id = "{0:04d}".format(all_spe.index(species_name)) + species_name+"_"+rec.id
            SeqIO.write(rec, PF00416_19, "fasta")

spe_not_contain_380 = set(all_spe) - set(spe_contain_380)
spe_not_contain_410 = set(all_spe) - set(spe_contain_410)
spe_not_contain_466 = set(all_spe) - set(spe_contain_466)
spe_not_contain_281 = set(all_spe) - set(spe_contain_281)
spe_not_contain_347 = set(all_spe) - set(spe_contain_347)
spe_not_contain_573 = set(all_spe) - set(spe_contain_573)
spe_not_contain_297 = set(all_spe) - set(spe_contain_297)
spe_not_contain_366 = set(all_spe) - set(spe_contain_366)
spe_not_contain_411 = set(all_spe) - set(spe_contain_411)
spe_not_contain_750 = set(all_spe) - set(spe_contain_750)
spe_not_contain_416 = set(all_spe) - set(spe_contain_416)


gap_seq = Seq("-"*50)


for rec in spe_not_contain_380:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00380")
    SeqIO.write(rec_generate_seq, PF00380_16, "fasta")
for rec in spe_not_contain_410:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00410")
    SeqIO.write(rec_generate_seq, PF00410_16, "fasta")
for rec in spe_not_contain_466:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00466")
    SeqIO.write(rec_generate_seq, PF00466_17, "fasta")
for rec in spe_not_contain_281:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00281")
    SeqIO.write(rec_generate_seq, PF00281_16, "fasta")
for rec in spe_not_contain_347:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00347")
    SeqIO.write(rec_generate_seq, PF00347_20, "fasta")
for rec in spe_not_contain_573:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00573")
    SeqIO.write(rec_generate_seq, PF00573_19, "fasta")
for rec in spe_not_contain_297:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00297")
    SeqIO.write(rec_generate_seq, PF00297_19, "fasta")
for rec in spe_not_contain_366:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00366")
    SeqIO.write(rec_generate_seq, PF00366_17, "fasta")
for rec in spe_not_contain_411:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00411")
    SeqIO.write(rec_generate_seq, PF00411_16, "fasta")
for rec in spe_not_contain_750:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00750")
    SeqIO.write(rec_generate_seq, PF00750_16, "fasta")
for rec in spe_not_contain_416:
    rec_generate_seq = SeqRecord(gap_seq, id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00416")
    SeqIO.write(rec_generate_seq, PF00416_19, "fasta")


PF00380_16.close()
PF00410_16.close()
PF00466_17.close()
PF00281_16.close()
PF00347_20.close()
PF00573_19.close()
PF00297_19.close()
PF00366_17.close()
PF00411_16.close()
PF00750_16.close()
PF00416_19.close()

con_fa =[  "conserved_gene_PF00380.fa",
           "conserved_gene_PF00410.fa",
           "conserved_gene_PF00466.fa",
           "conserved_gene_PF00281.fa",
           "conserved_gene_PF00347.fa",
           "conserved_gene_PF00573.fa",
           "conserved_gene_PF00297.fa",
           "conserved_gene_PF00366.fa",
           "conserved_gene_PF00411.fa",
           "conserved_gene_PF00750.fa",
           "conserved_gene_PF00416.fa"
        ]


time.sleep(5)
muscle = "/home/User1/software/muscle3.8.31_i86linux64"                        # set the executive application directory of Muscle
work_dir = "/home/User2/Python_Test/Genome/conserved_gene/result_fa"           # set the working dir of result fasta file
#print m_cline
#stdout, stderr = m_cline()

#align = AlignIO.read(StringIO(stdout), "clustal")
#align = AlignIO.read(StringIO(stdout), "fasta")
#print align

#result_handle = open("align_conserved_gene_PF00466.faa", "w")
#SeqIO.write(align, result_handle, "fasta")

aligned_pfam = []

for i in con_fa:
    FA_Conserved = i.split("_")[-1]
    aligned_pfam.append("align_" + FA_Conserved)
    m_cline = MuscleCommandline(muscle, input = work_dir + "/" + i, out = work_dir +"/align_"+FA_Conserved, clw=False)
    
    m_cline()

sequences = {}

# SN is the number of sequences in the alignment

#for i in aligned_pfam:
#    with open(i, "r") as seq_file:
#        for line in seq_file.readlines():
#            SN = line.rstrip()[0:4]
            
#            if line.startswith("0") and SN not in sequences:
#                sequences[SN] = line.rstrip()[4:]
#            elif line.startswith("0"):
#                sequences[SN] += line.rstrip().split()[1]

#final_align = open("all_align_new.fasta", "w")


#for num, seq in sequences.iteritems():
#    c = seq.split()
#    final_align.write(">"+ c[0]+"\n"+c[1]+ "\n")

#final_align.close()    


PF00380_16_al = open("alignGene_PF00380.fa", "w")
PF00410_16_al = open("alignGene_PF00410.fa", "w")
PF00466_17_al = open("alignGene_PF00466.fa", "w")
PF00281_16_al = open("alignGene_PF00281.fa", "w")
PF00347_20_al = open("alignGene_PF00347.fa", "w")
PF00573_19_al = open("alignGene_PF00573.fa", "w")
PF00297_19_al = open("alignGene_PF00297.fa", "w")
PF00366_17_al = open("alignGene_PF00366.fa", "w")
PF00411_16_al = open("alignGene_PF00411.fa", "w")
PF00750_16_al = open("alignGene_PF00750.fa", "w")
PF00416_19_al = open("alignGene_PF00416.fa", "w")

all_al = glob("align_*")

for i in all_al:
    records = SeqIO.parse(i, "fasta")  
    if i.find("PF00380") != -1:
        for rec in records:
            PF380_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00380_16_al, "fasta")
    elif i.find("PF00410") != -1:
        for rec in records:
            PF410_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00410_16_al, "fasta")
    elif i.find("PF00466") != -1:
        for rec in records:        
            PF466_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00466_17_al, "fasta")
    elif i.find("PF00281") != -1:
        for rec in records:
            PF281_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00281_16_al, "fasta")
    elif i.find("PF00347") != -1:
        for rec in records:    
            PF347_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00347_20_al, "fasta")
    elif i.find("PF00573") != -1:
        for rec in records:
            PF573_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00573_19_al, "fasta")
    elif i.find("PF00297") != -1:
        for rec in records:
            PF297_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00297_19_al, "fasta")
    elif i.find("PF00366") != -1:
        for rec in records:
            PF366_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00366_17_al, "fasta")
    elif i.find("PF00411") != -1:
        for rec in records:
            PF411_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00411_16_al, "fasta")
    elif i.find("PF00750") != -1:
        for rec in records:
            PF750_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00750_16_al, "fasta")
    elif i.find("PF00416") != -1:
        for rec in records:
            PF416_gap_length = len(rec.seq)
            SeqIO.write(rec, PF00416_19_al, "fasta")


for rec in spe_not_contain_380:
    rec_generate_seq = SeqRecord(Seq("-" * PF380_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00380")
    SeqIO.write(rec_generate_seq, PF00380_16_al, "fasta")
for rec in spe_not_contain_410:
    rec_generate_seq = SeqRecord(Seq("-" * PF410_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00410")
    SeqIO.write(rec_generate_seq, PF00410_16_al, "fasta")
for rec in spe_not_contain_466:
    rec_generate_seq = SeqRecord(Seq("-" * PF466_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00466")
    SeqIO.write(rec_generate_seq, PF00466_17_al, "fasta")
for rec in spe_not_contain_281:
    rec_generate_seq = SeqRecord(Seq("-" * PF281_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00281")
    SeqIO.write(rec_generate_seq, PF00281_16_al, "fasta")
for rec in spe_not_contain_347:
    rec_generate_seq = SeqRecord(Seq("-" * PF347_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00347")
    SeqIO.write(rec_generate_seq, PF00347_20_al, "fasta")
for rec in spe_not_contain_573:
    rec_generate_seq = SeqRecord(Seq("-" * PF573_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00573")
    SeqIO.write(rec_generate_seq, PF00573_19_al, "fasta")
for rec in spe_not_contain_297:
    rec_generate_seq = SeqRecord(Seq("-" * PF297_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00297")
    SeqIO.write(rec_generate_seq, PF00297_19_al, "fasta")
for rec in spe_not_contain_366:
    rec_generate_seq = SeqRecord(Seq("-" * PF366_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00366")
    SeqIO.write(rec_generate_seq, PF00366_17_al, "fasta")
for rec in spe_not_contain_411:
    rec_generate_seq = SeqRecord(Seq("-" * PF411_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00411")
    SeqIO.write(rec_generate_seq, PF00411_16_al, "fasta")
for rec in spe_not_contain_750:
    rec_generate_seq = SeqRecord(Seq("-" * PF750_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00750")
    SeqIO.write(rec_generate_seq, PF00750_16_al, "fasta")
for rec in spe_not_contain_416:
    rec_generate_seq = SeqRecord(Seq("-" * PF416_gap_length), id = "{0:04d}".format(all_spe.index(rec)) +rec, description ="#PF00416")
    SeqIO.write(rec_generate_seq, PF00416_19_al, "fasta")


PF00380_16_al.close()
PF00410_16_al.close()
PF00466_17_al.close()
PF00281_16_al.close()
PF00347_20_al.close()
PF00573_19_al.close()
PF00297_19_al.close()
PF00366_17_al.close()
PF00411_16_al.close()
PF00750_16_al.close()
PF00416_19_al.close()

new_align_pfam = []

for i in aligned_pfam:
    temp_al = "alignGene_" + i.split("_")[1] 
    new_align_pfam.append(temp_al)


for i in new_align_pfam:
    with open(i, "r") as seq_file:
        for line in seq_file.readlines():
            
            if line.startswith(">"):
                SN = line.rstrip()[1:5]
                sequences[SN] = line.rstrip()[5:33]

fa_content = {}

for i in new_align_pfam:
    for recs in SeqIO.parse(i , "fasta"):
        if recs.id[0:4] in sequences:
            seq_name = sequences.get(recs.id[0:4])
            if seq_name not in fa_content:
                fa_content[seq_name] = str(recs.seq)
            else:
                fa_content[seq_name] += str(recs.seq)

final_align = open("all_align_new.fasta", "w")


for name, seq in fa_content.iteritems():
    final_align.write(">"+ name+"\n"+seq+ "\n")

final_align.close()    

