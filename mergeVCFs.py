#!/usr/bin/env python 

########################################################################
# Script to merge multipe VCF files into a single VCF                  #
# This script is undoubtely amazing af but i wrote it for fun          #
# Use bcftools/vcftools if you actually want to merge vcf files        #
#                                                                      # 
# Usage example:                                                       #
# python mergeVCFs.py --vcf <.vcf> <.vcf> <.vcf>                       #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

import argparse,os
from collections import defaultdict

""" 
DICT POSOTIONS 
i[0] == #CHROM
i[1] == POS
i[2] == ID
i[3] == REF
i[4] == ALT
i[5] == QUAL
i[6] == FILTER
i[7] == FORMAT 
"""

headerSet = set()
headerDict = defaultdict(list)
def sortHeaders(header):
	for x in header:
		if x in headerSet:
			continue 
		headerType = x.split("=")[0]
		headerDict[headerType].append(x)
		headerSet.add(x) 

def mergeVCF(vcf):
	vcfDict = defaultdict(list)
	for x in vcf:
		with open(x, 'r') as file:
			vcfData = [i.strip() for i in file.readlines()]
			sortHeaders([i for i in vcfData if i.startswith("##")])

			for i in vcfData:
				if i.startswith("#"):
					continue
				i = i.split()
				key = "".join(i[0:2])

				if key in vcfDict:
					ID REF ALT QUAL FILTER FORMAT = i[3:9]
					
####ID 					
					if ID != vcfDict[key][0] and VCF[key][0] == ".":
						vcfDict[key][0] = ID
####REF				
					if REF != vcfDict[key][1]:
						vcfDict[key][1] = ",".join([vcfDict[key][1], REF])
####ALT 			
					if ALT != vcfDict[key][2]:
						vcfDict[key][2] = ",".join([vcfDict[key][2], ALT])
####QUAL
					if float(QUAL) > float(vcfDict[key][3]):
						vcfDict[key][3] = QUAL
####FILTER 
####FORMAT
					



					
				else:
					vcfDict["".join(i[0:2])].extend(i[3:])




	for x,y in vcfDict.items():
		print(x,y)



	





def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		required = True,
		nargs = '+',
		help = 'space sep list of VCF files to be merged'
		)
	parse = ap.parse_args()
	mergeVCF(parse.vcf)

if __name__ == '__main__':
	main()