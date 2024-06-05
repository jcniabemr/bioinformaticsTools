#!/usr/bin/env python 

###########################################################################################
# Script to identify hrp box moteifs in selected genes                                    #   
#                                                                                         #
# Usage example:                                                                          #
# python findHRP.py --assembly <.fasta> --gff <.gff>                                      #
#                                                                                         # 
# hrp boxes have a conserved GGCACC-16N-CCAC seq starting 35bp upstream of ATG            # 
# hrp seqs from https://onlinelibrary.wiley.com/doi/10.1046/j.1365-2958.2002.02964.x      # 
#                                                                                         #
# Written by John Connell                                                                 #
# john.connell@niab.com                                                                   #
# NIAB                                                                                    # 
# 2024                                                                                    #
###########################################################################################

def searchHRP(assembly, contig, start, gene, orient):
	if orient == "+":
		begin = max(1, start - 150)
		end = start
	else:
		begin = start 
		end = start + 150
	seq = ""
	seqCounter = 0
	found = False
	foundHits = []
	with open(assembly, 'r') as genomeAssembly:
		for line in genomeAssembly:
			line = line.strip()
			if line.startswith(">" + contig):
				found = True
				continue
			if found:
				for i in line:
					seqCounter += 1
					if seqCounter >= begin:
						seq += i
						if seqCounter == end:
							hrpA_Pos = re.compile(r'TGGAACC.{16}CCACCTA', re.IGNORECASE).finditer(seq)
							hrpA_Neg = re.compile(r'TAGGTGG.{16}GGTTCCA', re.IGNORECASE).finditer(seq)
							for x in hrpA_Pos:
								foundHits.append("\t".join(map(str,[contig,x.group(),"hrpA",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","+","-",gene])))
							for x in hrpA_Neg:
								foundHits.append("\t".join(map(str,[contig,x.group(),"hrpA",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","-","-",gene])))
							hrpJ_Pos = re.compile(r'GGGAACT.{16}CCACTCA', re.IGNORECASE).finditer(seq)
							hrpJ_Neg = re.compile(r'TGAGTGG.{16}AGTTCCC', re.IGNORECASE).finditer(seq)
							for x in hrpJ_Pos:
								foundHits.append("\t".join(map(str,[contig,x.group(),"hrpJ",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","+","-",gene])))
							for x in hrpJ_Neg:
								foundHits.append("\t".join(map(str,[contig,x.group(),"hrpJ",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","-","-",gene])))
							hrpW_Pos = re.compile(r'GGGAACC.{15}CCACTCA', re.IGNORECASE).finditer(seq)
							hrpW_Neg = re.compile(r'TGAGTGG.{15}GGTTCCC', re.IGNORECASE).finditer(seq)
							for x in hrpW_Pos:
								foundHits.append("\t".join(map(str,[contig,x.group(),"hrpW",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","+","-",gene])))
							for x in hrpW_Neg:
								foundHits.append("\t".join(map(str,[contig,x.group(),"hrpW",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","-","-",gene])))
							Consensus_Pos = re.compile(r'GGCACC.{16}CCAC', re.IGNORECASE).finditer(seq)
							Consensus_Neg = re.compile(r'GTGG.{16}GGTGCC', re.IGNORECASE).finditer(seq)
							for x in Consensus_Pos:
								foundHits.append("\t".join(map(str,[contig,x.group(),"Consensus",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","+","-",gene])))
							for x in Consensus_Neg:
								foundHits.append("\t".join(map(str,[contig,x.group(),"Consensus",str(start - 150 + x.start()),str(start - 150 + x.end()-1),"-","-","-",gene])))
							found = False 
							break		
	return foundHits

def extractRegions(assembly, gff):
	seq = ""
	found = False
	hits = []
	hits.append("##gff-version 3")
	with open(assembly, 'r') as genomeAssembly, open(gff, 'r') as gffFile, open(os.path.splitext(assembly)[0] + "_hrp_hits.gff", 'w') as out:
		gffData = [x.strip().split() for x in gffFile.readlines() if not x.startswith("#")]
		for line in gffData:
			if line[2] == "gene":
				contig = line[0]
				gene = line[8]
				found = True
				if line[6] == "+":
					hits.extend(searchHRP(assembly,contig,int(line[3]),gene,"+"))
				else:
					hits.extend(searchHRP(assembly,contig,int(line[4]),gene,"-"))
			else:
				continue
		for x in hits:
			out.write(f"{x}\n")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--assembly',
		type = str,
		required = True,
		help = 'genome assembly'
	)
	ap.add_argument(
		'--gff',
		type = str,
		required = True,
		help = 'gff file'
	)
	parse = ap.parse_args()
	extractRegions(parse.assembly,parse.gff)

if __name__ == '__main__':
	import os,argparse,re
	main()