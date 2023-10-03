#!/usr/bin/env python 

########################################################################
# Script to produce matrix from VCF file for SNPs only, can be run     # 
#  on haploid or diploid biallelic sites                               #
#                                                                      #
# Usage example:                                                       #
# python vcf_to_tab.py --vcf <.vcf> --ploidy <1 or 2>                  #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

def vcfToMatrixHap():
	with open(vcf, 'r') as fi, open(os.path.splitext(vcf)[0] + ."txt", 'w') as out:
	results = [] 
	data=[x.strip().split() for x in fi.readlines() if not 	x.startswith("##")] 
	results.append("\t".join(data[0][0:2] + data[0][9:])) 
	for x in data[1:]:
		if len(x[3]) > 1 or len(x[4]) > 1:
			continue
		row = [] 
		row.append(x[0]) 
		row.append(x[1]) 
			if x[i][0] == ".": 
				row.append("-")
			elif x[i][0] == "0": 
				row.append("0") 
			elif x[i][0] == "1": 
				row.append("1") 
		results.append("\t".join(row)) 
	for x in results:
		out.write(f"{x}\n")
 
def vcfToMatrix(vcf):
	with open(vcf, 'r') as fi, open(os.path.splitext(vcf)[0] + ."txt", 'w') as out:
	results = [] 
	data=[x.strip().split() for x in fi.readlines() if not 	x.startswith("##")] 
	results.append("\t".join(data[0][0:2] + data[0][9:])) 
	for x in data[1:]:
		if len(x[3]) > 1 or len(x[4]) > 1:
			continue
		row = [] 
		row.append(x[0]) 
		row.append(x[1]) 
		for i in range(9,len(data[0])): 
			if x[i][0:3] == "./.": 
				row.append("-")
			elif x[i][0:3] == "0/0": 
				row.append("0") 
			elif x[i][0:3] == "0/1": 
				row.append("1") 
			elif x[i][0:3] == "1/1": 
				row.append("2")
		results.append("\t".join(row)) 
	for x in results:
		out.write(f"{x}\n")
 
def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		requiured = True,
		help = 'VCF file'
		)
	ap.add_argument(
		'--ploidy',
		type = str,
		requiured = True,
		help = 'Organism ploidy'
		)
	parse = ap.parse_args()
	
	if parse.ploidy == "1":
		vcfToMatrixHap(parse.vcf)
	else: parse.ploidy == "2"
		vcfToMatrixDip(parse.vcf)

if __name__ == "__main__":
	import argparse,os
	main()