#!/usr/bin/env python

#######################################################################################
# Script to cut a genome into fragments that contains "N" characters                  #
#                                                                                     # 
# Options: --assembly --contig --step --window                                        #
#                                                                                     #
# Written by John Connell                                                             #
# john.connell@niab.com                                                               #
# NIAB                                                                                #
# 2024                                                                                #
#######################################################################################

def cutGenome(assembly,contig,step,window):
    found = False
    finished = False
    seq = ""
    lenCounter = 0 
    with open(assembly, 'r') as file, open(os.path.splitext(assembly)[0] + "_" + contig + "_" + str(window)+ "_" + str(step) + ".fasta", 'w') as out:
        for x in file:
            x = x.strip()
            if x.startswith(">"+contig):
                found = True
                continue 
            if found:
                for i in x:
                    lenCounter += 1
                    if lenCounter == start:
                        seq += i.upper()
                        if lenCounter == end:
                            finished = True
                            break
            if finished:
                break
        for i in range(0, len(seq), step):
            newEnd = min(i + window, len(seq))
            chunk = seq[i:newEnd]
            out.write(f">{assembly.split('.')[0].split('/')[-1]}_{contig}_{start + i-1}_{start + i + window-1}_{window}_{step}\n") 
            out.write("\n".join([chunk[o:o + 60] for o in range(0, len(chunk), 60)]) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '--assembly',
        type = str,
        required = True,
        help = 'Assembly to be cut'
    )
    ap.add_argument(
        '--step',
        type = int,
        required = True,
        help = 'Step length for cutting'
    )
    ap.add_argument(
        '--window',
        type = int,
        required = True,
        help = 'Length of read'
    )
    ap.add_argument(
        '--contig',
        type = str,
        required = True,
        help = 'Contig for cutting'
    )
    parse=ap.parse_args()
    cutGenome(
        parse.assembly,
        parse.contig,  
        parse.step,
        parse.window 
    )

if __name__ == '__main__':
    import os,argparse
    main()
