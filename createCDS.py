#!/usr/bin/env python 

########################################################################
# Script to produce DNA and or protein CDS fasta from gff and assembly #
# Usage example:                                                       #
# python createCDS.py --gff <.gff> --fasta <.fasta> --strainName <x>   #
# ***Note Type field should be ordered as"gene,mRNA,CDS" descending*** #
#             ***This can be achieved using sortGFFs.py***             #
# Written by John Connell                                                #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

####Import functions 
import argparse 

####Parse Args 
ap=argparse.ArgumentParser()
ap.add_argument('--gff',type=str,required=True,help='gff file')
ap.add_argument('--fasta',type=str,required=True,help='DNA infile')
ap.add_argument('--strainName',type=str,required=True,help='Strain name for outfile')
parse=ap.parse_args()

####Open files 
with open(parse.gff,'r') as f1:
    gffData=[x.strip().split() for x in f1.readlines() if not x.startswith("#")]
with open(parse.fasta,'r') as f2:
    fastaData=[x.strip() for x in f2.readlines()]
DNA_CDS_Out=open("_".join([parse.strainName,"DNA_CDS.fasta"]),'w')
PROTEIN_CDS_out=open("_".join([parse.strainName,"PROTEIN_CDS.fasta"]),'w') 

####Define functions 
def cutCDS(contig,start,stop):
    cutFasta=""
    contigFound=False
    for x in fastaData:
        if x.startswith(">"):
            if contigFound:
                break
            if x.split(">")[1]==contig:
                contigFound=True
                continue
        if contigFound:
            cutFasta+=x
    return cutFasta[start-1:stop]

def reverseComplement(seq):
    baseComplement={'G':'C','C':'G','A':'T','T':'A'}
    reversedSeq=seq[::-1]
    reversedComplemented=''.join(baseComplement[x] for x in reversedSeq)
    return reversedComplemented

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
            AA=translationTable.get(cSeq,"") 
            proteinSeq+=AA
            cSeq=""
    return proteinSeq

####Run conversion 
DNASeq=""
geneNUM=0
cdsNUM=0

for x in gffData:
    if x[2]!="gene" and geneNUM==0:
        print("Error: Type field must start with 'gene'!")
        break 
    if x[2]=="gene":
        if geneNUM==0:
            geneNUM+=1
            orient=x[6]
            continue 
        if geneNUM!=0:
            if orient=="-":
                DNASeq=reverseComplement(DNASeq)
            DNA_CDS_Out.write(".".join([">g"+str(geneNUM),"t"+str(cdsNUM)])+"\n")
            PROTEIN_CDS_out.write(".".join([">g"+str(geneNUM),"t"+str(cdsNUM)])+"\n")
            DNA_CDS_Out.write("\n".join([DNASeq[i:i+60] for i in range(0,len(DNASeq),60)])+"\n")
            PROTEINSeq=translate(DNASeq)
            PROTEIN_CDS_out.write("\n".join([PROTEINSeq[i:i+60] for i in range(0,len(PROTEINSeq),60)])+"\n")
            DNASeq=""
            geneNUM+=1
            cdsNUM=0
            orient=x[6]
            continue
    if x[2]=="mRNA":
        if cdsNUM==0:
            cdsNUM+=1
            continue 
        else:
            if cdsNUM!=0:
                if orient=="-":
                    DNASeq=reverseComplement(DNASeq)
                DNA_CDS_Out.write(".".join([">g"+str(geneNUM),"t"+str(cdsNUM)])+"\n")
                PROTEIN_CDS_out.write(".".join([">g"+str(geneNUM),"t"+str(cdsNUM)])+"\n")
                DNA_CDS_Out.write("\n".join([DNASeq[i:i+60] for i in range(0,len(DNASeq),60)])+"\n")
                PROTEINSeq=translate(DNASeq)
                PROTEIN_CDS_out.write("\n".join([PROTEINSeq[i:i+60] for i in range(0,len(PROTEINSeq),60)])+"\n")
                DNASeq=""
                cdsNUM+=1
                continue 
    if x[2]=="CDS":
            DNASeq+=cutCDS(x[0],int(x[3]),int(x[4]))
            continue 
if orient=="-":
    DNASeq=reverseComplement(DNASeq)
DNA_CDS_Out.write(".".join([">g"+str(geneNUM),"t"+str(cdsNUM)])+"\n")
PROTEIN_CDS_out.write(".".join([">g"+str(geneNUM),"t"+str(cdsNUM)])+"\n")
DNA_CDS_Out.write("\n".join([DNASeq[i:i+60] for i in range(0,len(DNASeq),60)])+"\n")
PROTEINSeq=translate(DNASeq)
PROTEIN_CDS_out.write("\n".join([PROTEINSeq[i:i+60] for i in range(0,len(PROTEINSeq),60)])+"\n")

DNA_CDS_Out.close()
PROTEIN_CDS_out.close()