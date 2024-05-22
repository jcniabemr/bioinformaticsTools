#!/usr/bin/env python 

#######################################################################################
# Script to cut a genome into fragments between a defined region at defined intervals #
#                                                                                     # 
# Options: --assembly --start --end --step --window --contig                          #
#                                                                                     #
# Written by John Connell                                                             #
# john.connell@niab.com                                                               #
# NIAB                                                                                #
# 2024                                                                                #
#######################################################################################

def cutGenome(assembly,contig,start,end,step,window):
    found = False
    finished = False
    seq = ""
    lenCounter = 0 
    with open(assembly, 'r') as file, open(os.path.splitext(assembly)[0] + "_" + contig + "_" + str(start) + "_" + str(end) + "_" + str(window)+ "_" + str(step) + ".txt", 'w') as out:
        for x in file:
            x = x.strip()
            if x.startswith(">"+contig):
                found = True
                continue 
            if found:
                for i in x:
                    lenCounter += 1
                    if lenCounter >= start:
                        seq += i.upper()
                        if lenCounter == end:
                            finished = True
                            break
            if finished:
                break
        for i in range(0, len(seq), step):
            start = start + i
            newEnd = min(start, + i + window, len(seq))
            out.write(f">{assembly.split('.')[0].split('/')[-1]}_{contig}_{start}_{start + newEnd}_{window}_{step}\n") 
            out.write("\n".join([seq[o:o + 60] for o in range(i, newEnd, 60)]) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '--assembly',
        type = str,
        required = True,
        help = 'Assembly to be cut'
    )
    ap.add_argument(
        '--start',
        type = int,
        required = True,
        help = 'Start position for cutting'
    )
    ap.add_argument(
        '--end',
        type = int,
        required = True,
        help = 'End position for cutting'
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
        parse.start, 
        parse.end, 
        parse.step,
        parse.window 
    )

if __name__ == '__main__':
    import os,argparse
    main()


