#!/usr/bin/env python 

#############################################################################
# Script to identify MIMP sequences in Fusarium oxysporum genome assemblys  #   
#                                                                           #
# Usage example:                                                            #
# python findMimpSequences.py --fasta <.fasta>                              #
#                                                                           #
# Written by John Connell                                                   #
# john.connell@niab.com                                                     #
# NIAB                                                                      #
# 2023                                                                      #
#############################################################################

def writeGFF(sortedhits, filename):
    c = 0
    with open(os.path.splitext(filename)[0] + "_MIMP_hits.gff", 'w') as file:
        file.write("\t".join(["Contig", "Tool", "MIMP_motif", "Start", "End", "Score", "Strand", "Phase", "INFO"]) + "\n")
        for x in sortedhits:
            contig, start, stop, seq, strand = x.strip().split()
            c += 1
            file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, "findMimpSequences.py", seq, start, stop, ".", strand, ".", "ID=MIMP_" + str(c)))

def mimpSearch(contig,sequence):
    foundHits = []
    posHits = re.compile(r'CAGTGGG..GCAA[TA]AA').finditer(sequence)
    negHits = re.compile(r'TT[TA]TTGC..CCCACTG').finditer(sequence)
    for x in posHits:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"+"])))
    for x in negHits:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"-"])))
    return foundHits

def processFile(fasta):
    hits = []
    with open(fasta, 'r') as fi:
        genome = [x.strip() for x in fi.readlines() if not x.startswith("#")]
    contig = ''
    seq = ''
    for x in genome:
        if x.startswith(">"):
            if contig:
                hits.extend(mimpSearch(contig,seq))
                contig = x
                seq = ''
                continue 
            contig = x
            continue
        seq += x
    hits.extend(mimpSearch(contig,seq))
    writeGFF(sorted(hits, key = lambda i: (int(i.split()[0].split("_")[1]),int(i.split()[1]))),fasta)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument(
        '--fasta',
        type = str,
        required = True,
        help = 'Genome assembly'
        )
    parse = ap.parse_args()
    processFile(parse.fasta)

if __name__=="__main__":
    import argparse,re,os
    main()