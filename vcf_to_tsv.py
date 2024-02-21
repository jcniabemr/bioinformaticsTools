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
	with open(vcf, 'r') as f, open(os.path.splitext(vcf)[0] + ".tsv", 'w') as out:
		for x in f.readlines().strip():
			if x.startswith("##"):
				out.write(f"{x}\n")
				continue
			if x.startwith("#CHROM"):
				x=x.split()
				out.write("\t".join([x[0], x[9:], x[2], x[3:5], x[]]))


def main():
	import os,argparse 
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
	main()
