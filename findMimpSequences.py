#!/usr/bin/env python 

#############################################################################
# Script to identify MIMP sequences in Fusarium oxysporum genome assemblys  #   
#                                                                           #
# Usage example:                                                            #
# python findMimpSequences.py --fasta <.fasta>                              #
#                                                                           #
# MIMP seq criteria taken from the paper:                                   #
# Genome-Wide Analysis of the Fusarium oxysporum mimp Family of MITEs       #
# and Mobilization of Both Native and De Novo Created mimps                 #
#                                                                           #
# Written by John Connell                                                   #
# john.connell@niab.com                                                     #
# NIAB                                                                      #
# 2023                                                                      #
#############################################################################

def writeGFF(sortedhits, filename):
    c = 0
    processedHits = set()
    with open(os.path.splitext(filename)[0] + "_MIMP_hits.gff", 'w') as file:
        file.write("\t".join(["#Contig","HitStrictness","MIMP_motif","Start","End","Score","Strand","Phase","INFO"]) + "\n")
        for x in sortedhits:
            contig,start,stop,seq,strand,strictness = x.strip().split()
            if "_".join([contig,start]) not in processedHits:
                processedHits.add("_".join([contig,start]))
                c += 1
                file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, strictness, seq, start, stop, ".", strand, ".", "ID=MIMP_" + str(c)))
            else:
                pass

def mimpSearch(contig,sequence):
    foundHits = []
    strictHitsPos = re.compile(r'CAGTGGGGTGCAATAAGTTTGAATACA').finditer(sequence)
    strictHitsNeg = re.compile(r'TGTATTCAAACTTATTGCACCCCACTG').finditer(sequence)
    moderateHitsPos = re.compile(r'CAGT[GA][GA]G[GA][TG]GCAA.AAGT[TA]T.AAT.C[AC]').finditer(sequence)
    moderateHitsNeg = re.compile(r'[GT]G.ATT.A[TA]ACTT.TTGC[CA][TC]C[TC][TC]ACTG').finditer(sequence)
    looseHitsPos = re.compile(r'CA.T..G..GCAA..A.T.T.......').finditer(sequence)
    looseHitsNeg = re.compile(r'.......A.A.T..TTGC..C..A.TG').finditer(sequence)
    for x in strictHitsPos:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"+","Strict"])))
    for x in strictHitsNeg:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"-","Strict"])))
    for x in moderateHitsPos:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"+","Moderate"])))
    for x in moderateHitsNeg:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"-","Moderate"])))
    for x in looseHitsPos:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"+","Loose"])))
    for x in looseHitsNeg:
        foundHits.append("\t".join(map(str,[contig[1:],x.start(),x.end(),x.group(),"-","Loose"])))
    return foundHits

def strictness_order(strictness):
    order = {'Strict': 1, 'Moderate': 2, 'Loose': 3}
    return order.get(strictness, 0)

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
    writeGFF(sorted(hits, key = lambda i: (int(i.split()[0].split("_")[1]),int(i.split()[1]), strictness_order(i[4]))),fasta)

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