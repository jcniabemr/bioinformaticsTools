#!/usr/bin/env python 

########################################################################
# Script to produce similarity matrix of similair sites for all        #
# samples in a VCF file                                                #
#                                                                      #
# Usage example:                                                       #
# python vcf_similarity_matrix.py --vcf <.vcf>                         #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################



def similarityMatrix(vcf):
	with open(vcf, 'r') as fi:
		data = [x.strip().split() for x in fi.readlines() if not x.startswith("##")]
		matResults = []
		matResults.append("\t".join([""] + data[0][9:]))
		# for x in data[1:]:
		# 	for i in range(9, len(data[0]):




def main():
	ap=argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		required = True,
		type = str,
		help = 'VCF file'
		)
	parse = ap.parse_args()
	similarityMatrix(parse.vcf)

if __name__ == "__main__":
	import argparse, os
	main()