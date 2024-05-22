#!/usr/bin/env python 

##############################################
# Script to convert a vcf to tsv format      #
#                                            #
# Usage example:                             #
# python vcf_to_tsv.py --vcf <.vcf>          #
#                                            #
# Written by John Connell                    #
# john.connell@niab.com                      #
# NIAB                                       #
# 2024                                       #
##############################################


def convertTSV(vcf):
	vcfDict = defaultdict(list)
	data = []
	with open(vcf, 'r') as f, open(os.path.splitext(vcf)[0] + ".tsv", 'w') as out:
		for x in f:
			if x.startswith("##"):
				continue
			#	out.write(f"{x}")
			else:
				data.append(x)
		splitData = [x.strip().split() for x in data]
		for x in splitData[1:]:
			for i in range(9, len(splitData[0])):
				if x[i][0] == "1":
					vcfDict["_".join(x[0:2] + x[3:5])].append(splitData[0][i])
		for x,y in vcfDict.items():
			j = "\t".join(map(str, y))
			out.write(f"{x}\t{j}\n")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		required = True,
		type = str,
		help = 'Input VCF file'
		)
	parse = ap.parse_args()
	convertTSV(parse.vcf)

if __name__ == "__main__":
	import os,argparse
	from collections import defaultdict 
	main()
