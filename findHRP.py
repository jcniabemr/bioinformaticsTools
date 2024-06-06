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
		begin = max(1, start - 320)
		end = start
	else:
		begin = start 
		end = start + 320
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
							hrpA_Pos = re.compile(r'TGGAACC.{15,16}CCACCTA', re.IGNORECASE).finditer(seq)
							hrpA_Neg = re.compile(r'TAGGTGG.{15,16}GGTTCCA', re.IGNORECASE).finditer(seq)
							for x in hrpA_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpA"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpA"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpA_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpA"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpA"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							hrpF_Pos = re.compile(r'TGGAACC.{15,16}CCACTCA', re.IGNORECASE).finditer(seq)
							hrpF_Neg = re.compile(r'TGAGTGG.{15,16}GGTTCCA', re.IGNORECASE).finditer(seq)
							for x in hrpF_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpF"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpF"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpF_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpF"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpF"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							hrpP_Pos = re.compile(r'TGGAACT.{15,16}CCACTTA', re.IGNORECASE).finditer(seq)
							hrpP_Neg = re.compile(r'TAAGTGG.{15,16}AGTTCCA', re.IGNORECASE).finditer(seq)
							for x in hrpP_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpP"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpP"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpP_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpP"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpP"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							hrpJ_Pos = re.compile(r'GGGAACT.{15,16}CCACTCA', re.IGNORECASE).finditer(seq)
							hrpJ_Neg = re.compile(r'TGAGTGG.{15,16}AGTTCCC', re.IGNORECASE).finditer(seq)
							for x in hrpJ_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpJ"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpJ"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpJ_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpJ"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpJ"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							hrpW_Pos = re.compile(r'GGGAACC.{15,16}CCACTCA', re.IGNORECASE).finditer(seq)
							hrpW_Neg = re.compile(r'TGAGTGG.{15,16}GGTTCCC', re.IGNORECASE).finditer(seq)
							for x in hrpW_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpW"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpW"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpW_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpW"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpW"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							hrpK_Pos = re.compile(r'TGGAACC.{15,16}CCACACA', re.IGNORECASE).finditer(seq)
							hrpK_Neg = re.compile(r'TGTGTGG.{15,16}GGTTCCA', re.IGNORECASE).finditer(seq)
							for x in hrpK_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpK"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpK"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpK_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpK"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpK"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							hrpU_Pos = re.compile(r'GTGGAAC.{15,16}CCACTTA', re.IGNORECASE).finditer(seq)
							hrpU_Neg = re.compile(r'TAAGTGG.{15,16}GTTCCAC', re.IGNORECASE).finditer(seq)
							for x in hrpU_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpU"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpU"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in hrpU_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpU"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"hrpU"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							Consensus_Pos = re.compile(r'GGAAC[CT].{15,17}CCAC]', re.IGNORECASE).finditer(seq)
							Consensus_Neg = re.compile(r'GTGG.{15,17}[AG]GTGCC', re.IGNORECASE).finditer(seq)
							for x in Consensus_Pos:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"Consensus"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","+","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"Consensus"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","+","-",gene])))
							for x in Consensus_Neg:
								if orient == "+":
									foundHits.append("\t".join(map(str,[contig,x.group(),"Consensus"+":"+str(320 - x.start()),str(start - 320 + x.start()),str(start - 320 + x.end()-1),"-","-","-",gene])))
								else:
									foundHits.append("\t".join(map(str,[contig,x.group(),"Consensus"+":"+str(x.start()),str(start + x.start()),str(start + x.end()-1),"-","-","-",gene])))
							found = False 
							break		
	return foundHits

def extractRegions(assembly, gff):
	seq = ""
	found = False
	hits = []
	hits.append("##gff-version 3")
	with open(assembly, 'r') as genomeAssembly, open(gff, 'r') as gffFile, open(os.path.splitext(assembly)[0] + "_hrp_hits.gff", 'w') as out:
		gffData = [x.strip().split() for x in gffFile.readlines() if x.startswith("NODE")]
		for line in gffData:
			if line[2] == "CDS":
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