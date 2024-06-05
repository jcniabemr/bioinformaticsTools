#!/usr/bin/env python 




def translate(DNAseq):
    translationTable={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "CAT": "H", "CAC": "H", 
        "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", 
        "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", 
        "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", 
        "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", 
        "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", 
        "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", 
        "GGG": "G", "TGA": "*", "TAG": "*", "TAA": "*",
    }
    proteinSeq=""
    cSeq=""
    for x in DNAseq:
        cSeq+=x
        if len(cSeq)==3:
            AA=translationTable.get(cSeq.upper(),"") 
            proteinSeq+=AA
            cSeq=""
    return proteinSeq

def main():
    ap = argparse.ArgumentParse()
    ap.add_argument(
        '--seq',
        required = True,
        type = str,
        help = 'seq to be analysed'
        )
    ap.parse_args()
    findAlternateSeqs(parse.seq)

if __name__ == '__main__':
    import argparse,os
    main()