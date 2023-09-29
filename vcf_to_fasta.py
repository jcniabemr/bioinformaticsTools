#!/usr/bin/env python 

########################################################################
# Script to produce fasta for each sample in a VCF file for use in     #
# SNP based phylogeny                                                  #
#                                                                      #
# Usage example:                                                       #
# python vcf_to_fasta.py --vcf <.vcf> --ploidy <1/2>                   #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

def writeFile(res,vcf):
	with open(os.path.splitext(vcf)[0] + ".fasta", 'w') as out:
		for x in res:
			out.write(f"{x}\n")

def defineCodes(seq):
	iupac = {
	'G/G': 'G', 'C/C': 'C', 'T/T': 'T', 'A/A': 'A',
	'G/T': 'K', 'T/G': 'K', 'A/C': 'M', 'C/A': 'M',
	'C/G': 'S', 'G/C': 'S', 'A/G': 'R', 'G/A': 'R',
	'A/T': 'W', 'T/A': 'W', 'C/T': 'Y', 'T/C': 'Y'
	}
	return iupac.get(seq,"-")

def createFastaDiploid(vcf):
	with open(vcf, 'r') as fi:
		data=[x.strip().split() for x in fi.readlines() if not x.startswith("##")]
		fastaFile = []
		for i in range(9, len(data[0])):
			fastaFile.append("".join([">",data[0][i]]))
			seq = ""
			for x in data[1:]:
				if len(x[3]) > 1 or len(x[4]) > 1:
					continue
				ref = x[3]
				alt = x[4]
				if x[i][0] == ".":
					seq += "-"
				if x[i][0] == "0" and x[i][2] == "0":
					seq += defineCodes(ref+"/"+ref)
				if x[i][0] == "0" and x[i][2] == "1":
					seq += defineCodes(ref+"/"+alt)
				if x[i][0] == "1" and x[i][2] == "1":
					seq += defineCodes(alt+"/"+alt)
			fastaFile.append(seq)
	writeFile(fastaFile,vcf)

def createFastaHaploid(vcf):
	with open(vcf, 'r') as fi:
		data=[x.strip().split() for x in fi.readlines() if not x.startswith("##")]
		fastaFile = [] 
		for i in range(9, len(data[0])):
			fastaFile.append("".join([">",data[0][i]]))
			seq = ""
			for x in data[1:]:
				if len(x[3]) > 1 or len(x[4]) > 1:
					continue 
				if x[i][0] == ".":
					seq += "-"
				elif x[i][0] == "0":
					seq += x[3]
				elif x[i][0] == "1":
					seq += x[4]
			fastaFile.append(seq)
	writeFile(fastaFile,vcf)

def main():
	ap=argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		required = True,
		help = 'VCF file'
		)
	ap.add_argument(
		'--ploidy',
		type = str,
		required = True,
		help = 'Organism ploidy'
		)
	parse=ap.parse_args()
	if parse.ploidy	== "1":
		createFastaHaploid(parse.vcf)
	else:
		createFastaDiploid(parse.vcf)

if __name__ == "__main__":
	import argparse,os
	main()
