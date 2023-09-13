#!/usr/bin/env python 

########################################################################
# Script to produce reads from a genome assembly, can be used for      #
# whole genome alignement where no short reads are available.          #
#                                                                      # 
# Usage example:                                                       #
# python genomeToReads.py --g <.fasta> --readL <300> --readS <50>      #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

def reverse_complement(x):
    complement={
            "A":"T", 
            "T":"A", 
            "C":"G", 
            "G":"C", 
            "N":"N"
    }
    return "".join(complement[i] for i in reversed(x))

def cut_genome(seq,read_length,step,header,forward,reverse):
    for i in range(0, len(seq) - (read_length - 1), step):
        region = f"{i}-{i + read_length}"
        headerName = f"@{header[1:]} {region}" 
        sequence = seq[i:i + read_length]
        rev_comp_seq = reverse_complement(sequence)
####Write forward reads         
        forward.write(" ".join([headerName,"TSFR\n"]))
        forward.write(f"{sequence}\n")
        forward.write(f"+\n{'F' * read_length}\n")
        forward.write(" ".join([headerName,"RCFR\n"]))
        forward.write(f"{rev_comp_seq}\n")
        forward.write(f"+\n{'F' * read_length}\n")
#####Write reverse         
        reverse.write(" ".join([headerName,"TSRR\n"]))
        reverse.write(f"{sequence}\n")
        reverse.write(f"+\n{'F' * read_length}\n")
        reverse.write(" ".join([headerName,"RCRR\n"]))
        reverse.write(f"{rev_comp_seq}\n")
        reverse.write(f"+\n{'F' * read_length}\n")

def extractRun(genome,length,step):
    forward=open(os.path.splitext(genome)[0] + "_1.fastq", 'w')
    reverse=open(os.path.splitext(genome)[0] + "_2.fastq", 'w')
    with open(genome,'r') as fi:
        data = [x.strip() for x in fi.readlines() if not x.startswith("#")]
    fasta = ""
    for x in data:
        if x.startswith(">"):
            if fasta:
                cut_genome(fasta,length,step,header,forward,reverse)
                header = x
                fasta = ""
            else: 
                header = x
        else:
            fasta += x
    cut_genome(fasta,length,step,header,forward,reverse)
    forward.close()
    reverse.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '--g',
        type = str,
        required = True,
        help = 'Input genone'
        )
    ap.add_argument(
        '--readL',
        type = int,
        required = True,
        help = 'Length of reads to cut'
        )
    ap.add_argument(
        '--readS',
        type = int,
        required = True,
        help = 'Lengh of read overlap'
        )
    parse = ap.parse_args()
    extractRun(
        parse.g,
        parse.readL,
        parse.readS
        )

if __name__=="__main__":
    import argparse,os
    main()