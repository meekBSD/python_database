#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
import sys
import time

Eleven_gene_fragments = [ "PF00380.16" , "PF00410.16" , "PF00466.17", "PF00281.16", "PF00347.20"  , "PF00573.19" ,"PF00297.19", "PF00366.17" , "PF00411.16", "PF00750.16", "PF00416.19" ]


all_spe_phylum = { 
         'AWOB01.1_Geoarchaeota': "Geoarchaeota" ,
         'Acaryochloris_marina_MBIC11017': "Cyanobacteria" ,
         'Acetobacter_pasteurianus_386B': "Proteobacteria",
         'Acetobacterium_woodii_DSM_1030':"Firmicutes" ,
         'Acetohalobium_arabaticum_DSM_5501': "Firmicutes" ,
         'Achromobacter_xylosoxidans_A8': "Proteobacteria" ,
         'Acidaminococcus_fermentans_DSM_20731': "Firmicutes" ,
         'Acidianus_hospitalis_W1': "Crenarchaeota" ,
         'Acidilobus_saccharovorans_345_15': "Crenarchaeota" ,
         'Acidimicrobidae_bacterium_YM16_304': "Actinobacteria" ,
         'Acidimicrobium_ferrooxidans_DSM_10331': "Actinobacteria" ,
         'Acidiphilium_cryptum_JF_5': "Proteobacteria" ,
         'Acidobacterium_capsulatum_ATCC_51196': "Acidobacteria" ,
         'Acidothermus_cellulolyticus_11B': "Actinobacteria" ,
         'Aciduliprofundum_boonei_T469': "Euryarchaeota" ,
         'Aequorivita_sublithincola_DSM_14238': "Bacteroidetes" ,
         'AeroCD12-1': "CD12" ,
         'Aeropyrum_camini_SY1': "Crenarchaeota" ,
         'Akkermansia_muciniphila_ATCC_BAA_835' : "Verrucomicrobiae" ,
         'Alistipes_finegoldii_DSM_17242': "Bacteroidetes",
         'Aminobacterium_colombiense_DSM_12261': "Synergistetes" ,
         'Anabaena_cylindrica_PCC_7122' : "Cyanobacteria" ,
         'Anaerobaculum_mobile_DSM_13181': "Synergistetes" ,
         'Anaerolinea_thermophila_UNI_1': "Chloroflexi" ,
         'Aquifex_aeolicus_VF5': "Aquificae" ,
         'Archaeoglobus_fulgidus_DSM_4304': "Euryarchaeota" ,
         'Arthrospira_platensis_NIES_39' : "Cyanobacteria" ,
         'Bacteroides_CF50' : "Bacteroidetes",
         'Borrelia_afzelii_HLJ01' : "Spirochaetes" ,
         'Brachyspira_hyodysenteriae_WA1' : "Spirochaetes" ,
         'Caldiarchaeum_subterraneum' : "Thaumarchaeota" ,
         'Caldilinea_aerophila_DSM_14535___NBRC_104270': "Chloroflexi" ,
         'Caldisericum_exile_AZM16c01' : "Caldiserica" ,
         'Calditerrivibrio_nitroreducens_DSM_19672' : "Deferribacteres" ,
         'Chlamydia_muridarum_Nigg' : "Chlamydiae" ,
         'Chlamydophila_abortus_S26_3': "Chlamydiae" ,
         'Chlorobaculum_parvum_NCIB_8327': "Chlorobi" ,
         'Chlorobium_chlorochromatii_CaD3' : "Chlorobi" ,
         'Chloroflexus_aggregans_DSM_9485' : "Chloroflexi",
         'Chloroherpeton_thalassium_ATCC_35110' : "Chlorobi" ,
         'Chthonomonas_calidirosea_T49' : "Armatimonadetes",
         'Coraliomargarita_akajimensis_DSM_45221' : "Verrucomicrobia" ,
         'Deferribacter_desulfuricans_SSM1' : "Deferribacteres" ,
         'Deinococcus_deserti_VCD115': "Deinococcus-Thermus" ,
         'Denitrovibrio_acetiphilus_DSM_12809' : "Deferribacteres",
         'Desulfurispirillum_indicum_S5' : "Chrysiogenetes" ,
         'Desulfurobacterium_thermolithotrophum_DSM_11699' : "Aquificae" ,
         'Dictyoglomus_thermophilum_H_6_12' : "Dictyoglomi",
         'Elusimicrobium_minutum_Pei191' : "Elusimicrobia" ,
         'Euryarchaeota_3' : "Euryarchaeota",
         'Fervidobacterium_nodosum_Rt17_B1' : "Thermotogae" ,
         'Fibrobacter_succinogenes_S85': "Fibrobacteres" ,
         'Firmicutes_cd12_1' : "CD12" ,
         'Fusobacterium_3_1_36A2' : "Fusobacteria" ,
         'GW2011_AR_AR10_43_0_Diapherotrites' : "Diapherotrites" ,
         'GW2011_AR_unknown_32_20_Woesearchaeota' : "Woesearchaeota" ,
         'GWA2_AR18_30_20_Woesearchaeota' : "Woesearchaeota" ,
         'GWA2_AR5_49_1_Aenigmarchaeota': "Aenigmarchaeota" ,
         'Gemmatimonas_aurantiaca_T_27' : "Gemmatimonadetes" ,
         'Haloarcula_sp_CBA1115' : "Euryarchaeota",
         'Halyomorpha_halys_symbiont' : "Proteobacteria" ,
         'Hydrogenobacter_thermophilus_TK_6' : "Aquificae" ,
         'Ignavibacterium_album_JCM_16511' : "Ignavibacteriae" ,
         'Ilyobacter_polytropus_DSM_2926' : "Fusobacteria" ,
         'Isosphaera_pallida_ATCC_43644' : "Planctomycetes" ,
         'Korarchaeum_cryptofilum_OPF8' : "Korarchaeota" ,
         'Koribacter_versatilis_Ellin345' : "Acidobacteria",
         'Kosmotoga_olearia_TBF_19_5_1' : "Thermotogae" ,
         'LFWW01.1_Bathyarchaeota' : "Bathyarchaeota" ,
         'LKMY01.1' : "Nanohaloarchaeota" ,
         'LUCB01.1_Bathyarchaeota' : "Bathyarchaeota" ,
         'LUCE01.1_Bathyarchaeota' : "Bathyarchaeota",
         'Leptospirillum_ferriphilum_ML_04' : "Nitrospirae" ,
         'Leptotrichia_buccalis_C_1013_b' : "Fusobacteria" ,
         'Lokiarchaeum_sp_GC14_75' : "Lokiarchaeota" ,
         'Marinithermus_hydrothermalis_DSM_14884' : "Deinococcus-Thermus" ,
         'Marinitoga_piezophila_KA3' : "Thermotogae" ,
         'Meiothermus_ruber_DSM_1279' : "Deinococcus-Thermus" ,
         'Melioribacter_roseus_P3M' : "Ignavibacteriae",
         'Metallosphaera_cuprina_Ar-4' : "Crenarchaeota" ,
         'Methanocaldococcus_infernus_ME' : "Euryarchaeota" ,
         'Methanothermobacter_thermautotrophicus_str_Delta_H' : "Euryarchaeota" ,
         'Methylacidiphilum_infernorum_V4': "Verrucomicrobia" ,
         'Methylomirabilis_oxyfera' : "NC10" ,
         'Micrarchaeum_acidiphilum_ARMAN-2' : "Micrarchaeota" ,
         'Nanoarchaeota_archaeon_7A' : "Nanoarchaeota",
         'Nanoarchaeote_Nst1_Nst1_C1' : "Nanoarchaeota" ,
         'Nanoarchaeum_equitans_Kin4_M' : "Nanoarchaeota" ,
         'Nanosalina_sp_J07AB43' : "Nanohaloarchaeota" ,
         'Nanosalinarum_sp_J07AB56': "Nanohaloarchaeota" ,
         'New_Euryarchaeota_1' : "Euryarchaeota",
         'New_Firmicutes_1' : "CD12" ,
         'New_unclassified_archaea' : "unclassified_Archaea" ,
         'Nitrosopelagicus_brevis_strain_CN25' : "Thaumarchaeota" ,
         'Nitrosopumilus_AR2' : "Thaumarchaeota" ,
         'Nitrososphaera_gargensis_Ga9_2' : "Thaumarchaeota" ,
         'Nitrososphaera_viennensis_EN76' : "Thaumarchaeota" ,
         'Nitrosotenuis_cloacae_strain_SAT1' : "Thaumarchaeota" ,
         'Nitrospira_defluvii' : "Nitrospirae" ,
         'Parvarchaeum_acidophilus_ARMAN-5' : "Parvarchaeota" ,
         'Phycisphaera_mikurensis_NBRC_102666' : "Planctomycetes" ,
         'Pirellula_staleyi_DSM_6068' : "Planctomycetes" ,
         'Protochlamydia_amoebophila_UWE25' : "Chlamydiae" ,
         'Pyrobaculum_sp_WP30' : "Crenarchaeota" ,
         'SR1_bacterium_RAAC1_SR1_1' : "SR1" ,
         'Saccharibacteria_bacterium_RAAC3_TM7_1' : "Saccharibacteria" ,
         'Saccharobacterium_alaburgensis' : "Saccharibacteria" ,
         'Solibacter_usitatus_Ellin6076' : "Acidobacteria" ,
         'Sphaerochaeta_pleomorpha_Grapes' : "Spirochaetes" ,
         'Sulfolobus_islandicus' : "Crenarchaeota" ,
         'Synergistetes_bacterium_SGP1' : "Synergistetes" ,
         'TA06TA06TA06_TA06_TA06' : "TA06",
         'TA06_bacterium_DG_24_GCA' : "TA06" ,
         'TA06_bacterium_DG_26_GCA' : "TA06" ,
         'TA06_bacterium_DG_78_GCA' : "TA06" ,
         'TA06_bacterium_SM1_40_GCA' : "TA06" ,
         'TA06_bacterium_SM23_40_GCA' : "TA06" ,
         'Thermodesulfatator_indicus_DSM_15286' : "Thermodesulfobacteria" ,
         'Thermodesulfobacterium_OPB45' : "Thermodesulfobacteria" ,
         'Thermodesulfovibrio_yellowstonii_DSM_11347' : "Nitrospirae" ,
         'Thermofilum_sp_1910b' : "Crenarchaeota" ,
         'WWE3_bacterium_RAAC2_WWE3_1' : "WWE3" ,
         'archaeon_Mx1201' : "Euryarchaeota",
         'halophilic_archaeon_DL31' : "unclassified_Archaea"             
         }

#records = SeqIO.parse("clustal_muscle_mafft.fas", "fasta")
fa_dict = SeqIO.to_dict(SeqIO.parse("clustal_muscle_mafft.fas", "fasta" ))

result = open("result_clustal_muscle_mafft.fas", "w")

n = 1

for i in all_spe_phylum:
    if len(i) <28 :
        for kc, vc in fa_dict.iteritems():
            if kc.find(i) != -1:
                vc.id = str(n) +"_" + all_spe_phylum.get(i, "_")+ "_"+i
                SeqIO.write(vc, result, "fasta")
                n += 1
                break
    if len(i) >= 28:
        for kc, vc in fa_dict.iteritems():
            if i.find(kc) != -1:
                vc.id = str(n) +"_" + all_spe_phylum.get(i, "_")+ "_"+i
                SeqIO.write(vc, result, "fasta")
                n += 1
                break

result.close()

'''
for i in all_txts:
    if any(i.find(j)!= -1 for j in all_items):
        continue
    else:
       print i    

for i in all_spe:
    for j in all_txts:
        if j.find(i):
            break
    else:
        print i ,"not found"

l = sorted(all_spe)

for i in l:
    print repr(i), ","

'''


