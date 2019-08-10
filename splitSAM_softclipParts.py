#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import argparse

#print(os.path.getsize("Sample_19.sam"))

d_EM = {}

samA = open("Sample_19.sam", "r")

for line in samA:
    xA = line.strip().split("\t")
    ID = xA[0]
    readLen = len(xA[9])
    CIGAR = xA[5]
    mLs = CIGAR.find("M")
    mRs = CIGAR.rfind("M")
    ##### soft clipping #####
    sLs = CIGAR.find("S") 
    sRs = CIGAR.rfind("S")
    
    if sLs != sRs:

        #Lseqpos = int(CIGAR[:sLs])
        Lseqpos = CIGAR[:sLs]
        Rseqpos = int(CIGAR[mRs+1:sRs])    
        if Lseqpos < Rseqpos:
            d_EM[ID] = "r" + str(readLen - Rseqpos)
        else:
            d_EM[ID] = "l" + Lseqpos
            
    else:
        if sLs != -1 and sRs > mRs:
            Rsp = int(CIGAR[mRs+1:sRs])
            d_EM[ID] = "r" + str(readLen - Rsp)
        elif sLs != -1 and sLs < mLs:
            #Lsp = int(CIGAR[:sLs])
            Lsp = CIGAR[:sLs]
            d_EM[ID] = "l" + Lsp

    ##### hard clipping #####

    hLs = CIGAR.find("H")
    hRs = CIGAR.rfind("H")

    if hLs != hRs:

        #Lseqpos = int(CIGAR[:hLs])
        Lseqpos = CIGAR[:hLs]
        Rseqpos = int(CIGAR[mRs+1:hRs])
        if Lseqpos < Rseqpos:
            d_EM[ID] = "r" + str(readLen + int(Lseqpos) - Rseqpos)
        else:
            d_EM[ID] = "l" + Lseqpos

    else:
        if hLs != -1 and hRs > mRs:
            #Rsp = int(CIGAR[hRs+1:hRs])
            Rsp = CIGAR[hRs+1:hRs]
            d_EM[ID] = "r" + str(readLen)
        elif hLs != -1 and hLs < mLs:
            #Lsp = int(CIGAR[:hLs])
            Lsp = CIGAR[:hLs]
            d_EM[ID] = "l" + Lsp

samA.close()

RD = dict(list(filter(lambda x: x[1].startswith("r"), d_EM.items())))
LD = dict(list(filter(lambda x: x[1].startswith("l"), d_EM.items())))


ALSam = open("Sample_19AL_mapped.sam", "r")

OutTemp = open("pos_temp.txt", "w")

for a in ALSam:
    Xline = a.strip().split("\t")
    ID = Xline[0]
    readALen = len(Xline[9])

    if ID in d_EM:
        #print(ID)
        CIGAR = Xline[5]
        mLs = CIGAR.find("M")
        mRs = CIGAR.rfind("M")
        ##### soft clipping #####
        sLs = CIGAR.find("S")
        sRs = CIGAR.rfind("S")

        if sLs != sRs:

            Lseqpos = int(CIGAR[:sLs])
            #Lseqpos = CIGAR[:sLs]
            Rseqpos = int(CIGAR[mRs+1:sRs])    
            #Rseqpos = CIGAR[mRs+1:sRs]
            if Lseqpos < Rseqpos and ID in LD:
                #print("{0}\t{1}\t{2}".format(ID, LD[ID], Rseqpos))
                Rend_pos = readALen - Rseqpos
                peml4 = int(LD[ID].strip("l"))
                if peml4 -10 <= Rseqpos:
                    OutTemp.write("{0}\tR\t{1}\tL\t{2}\t{3}\t{4}\n".format(ID, peml4, Rend_pos, Xline[1], Xline[3]))
            elif ID in RD:
                #print("{0}\t{1}\t{2}".format(ID, RD[ID], Lseqpos))
                peml4 = int(RD[ID].strip("r"))
                if peml4 +10 >= Lseqpos:
                    OutTemp.write("{0}\tL\t{1}\tR\t{2}\t{3}\t{4}\n".format(ID, peml4, Lseqpos, Xline[1], Xline[3]))

        else:
            if sLs != -1 and sRs > mRs and ID in LD:
                Rsp = int(CIGAR[mRs+1:sRs])
                #Rsp = CIGAR[mRs+1:sRs]
                #print("{0}\t{1}\t{2}".format(ID, LD[ID], Rsp))
                Rend_pos = readALen - Rsp
                peml4 = int(LD[ID].strip("l"))
                if peml4 -10 <= Rsp:
                    OutTemp.write("{0}\tR\t{1}\tL\t{2}\t{3}\t{4}\n".format(ID, peml4, Rend_pos, Xline[1], Xline[3]))
            elif sLs != -1 and sLs < mLs and ID in RD:
                Lsp = int(CIGAR[:sLs])
                #Lsp = CIGAR[:sLs]
                #print("{0}\t{1}\t{2}".format(ID, RD[ID], Lsp))
                peml4 = int(RD[ID].strip("r"))
                if peml4 +10 >= Lsp:
                    OutTemp.write("{0}\tL\t{1}\tR\t{2}\t{3}\t{4}\n".format(ID, peml4, Lsp, Xline[1], Xline[3]))
        hLs = CIGAR.find("H")
        hRs = CIGAR.rfind("H")

        if hLs != hRs:
            Lseqpos = int(CIGAR[:hLs])
            #Lseqpos = CIGAR[:hLs]
            Rseqpos = int(CIGAR[mRs+1:hRs])
            #Rseqpos = CIGAR[mRs+1:hRs]
            if Lseqpos < Rseqpos and ID in LD:
                #print("{0}\t{1}\t{2}".format(ID, LD[ID], Rseqpos))
                HRend_pos = readALen + Lseqpos - Rseqpos
                peml4 = int(LD[ID].strip("l"))
                if peml4 -10 <= Rseqpos:
                    OutTemp.write("{0}\tR\t{1}\tL\t{2}\t{3}\t{4}\n".format(ID, peml4, HRend_pos, Xline[1], Xline[3]))
            elif ID in RD:
                #print("{0}\t{1}\t{2}".format(ID, RD[ID], Lseqpos))
                peml4 = int(RD[ID].strip("r"))
                if peml4 +10 >= Lseqpos:
                    OutTemp.write("{0}\tL\t{1}\tR\t{2}\t{3}\t{4}\n".format(ID, peml4, Lseqpos, Xline[1], Xline[3]))
        else:
            if hLs != -1 and hRs > mRs and ID in LD:
                Rsp = int(CIGAR[hRs+1:hRs])
                #Rsp = CIGAR[hRs+1:hRs]
                #print("{0}\t{1}\t{2}".format(ID, LD[ID], Rsp))
                HRend_pos = readALen
                peml4 = int(LD[ID].strip("l"))
                if peml4 -10 <= Rsp:
                    OutTemp.write("{0}\tR\t{1}\tL\t{2}\t{3}\t{4}\n".format(ID, peml4, HRend_pos, Xline[1], Xline[3]))
            elif hLs != -1 and hLs < mLs and ID in RD:
                Lsp = int(CIGAR[:hLs])
                #Lsp = CIGAR[:hLs]
                #print("{0}\t{1}\t{2}".format(ID, RD[ID], Lsp))
                peml4 = int(RD[ID].strip("r"))
                if peml4 +10 >= Lsp:
                    OutTemp.write("{0}\tR\t{1}\tL\t{2}\t{3}\t{4}\n".format(ID, peml4, Lsp, Xline[1], Xline[3]))

ALSam.close()
OutTemp.close()

fileH = open("pos_temp.txt", "r")
dH = {}

for line in fileH:
    s = line.strip().split("\t",1)
    dH[s[0]]= s[1]
fileH.close()


ReSA = open("Sample_19.sam", "r")

for line in ReSA:
    xs = line.strip().split("\t")
    readsname = xs[0]
    if readsname in dH:
        print("{0}\t{1}\t{2}\t{3}".format(readsname, xs[3], xs[1], dH[readsname]))

ReSA.close()


