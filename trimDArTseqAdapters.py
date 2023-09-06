#!/usr/bin/env python 

####################################################################################
# Script to trim DArTseq adapters from DArTseq reads                              #
#                                                                                 #
# Usage example:                                                                  #
# python trimDArTseqAdapters.py --reads f1.FASTQ.gz f2.FASTQ.gz fX.FASTQ.gz       #
#                                                                                 #
# Written by John Connell                                                         #
# john.connell@niab.com                                                           #
# NIAB                                                                            #
# 2023                                                                            #
###################################################################################

def sepReads(reads):
    counter = 0
    seq_list = []
    for x in reads:
        if (counter - 1) % 4 == 0:
            seq_list.append(x)
        counter += 1
        if counter == 10000:
            break
    return seq_list

def baseFrequency(seq):
    countList = []
    for x in range(len(seq[0])):
        baseCounts = {
            'N': 0,
            'A': 0, 
            'T': 0, 
            'G': 0, 
            'C': 0 
        }
        for y in seq:
            BASE = y[x]
            baseCounts[BASE] += 1
        mostCommonNucleotide = max(baseCounts, key=baseCounts.get)
        BF = baseCounts[mostCommonNucleotide] / len(seq)
        countList.append(BF)
    dartBases = 0
    for i in countList:
        if i >= 0.90:
            dartBases += 1
        else:
            break
    return dartBases

def writeRes(a,b):
    newFileName = os.path.splitext(a)[0] + "_barcodeTrimmed.gz"
    with gzip.open(newFileName, 'wt') as out:
        for x in b:
            out.write(f"{x}\n")

def trimBarcodes(reads):
    for file in reads:
        with gzip.open(file, 'rt') as fi:
            trimmedAdapters = []
            seqData=[x.strip() for x in fi.readlines()]
            position=baseFrequency(sepReads(seqData))
            counter = 0
            for x in seqData:
                if counter % 4 == 0:
                    trimmedAdapters.append(x) 
                elif (counter - 1) % 4 == 0:
                    trimmedAdapters.append(x[position:])
                elif (counter - 2) % 4 == 0:
                    trimmedAdapters.append(x)
                else:
                    trimmedAdapters.append(x)
                counter += 1
        writeRes(file,trimmedAdapters)
        
def main(): 
    ap=argparse.ArgumentParser()
    ap.add_argument(
        '--reads',
        type=str,
        required=True,
        nargs="+",
        help="DArTseq read file, comma seperated list"
        )
    parse=ap.parse_args()
    trimBarcodes(parse.reads)

import argparse
import os
import gzip
if __name__=="__main__":
    main()