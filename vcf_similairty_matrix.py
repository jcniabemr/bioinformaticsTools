#!/usr/bin/env python 

########################################################################
# Script to produce similarity matrix for all samples vs all samples   #
# in a VCF file                                                        #
#                                                                      #
# Usage example:                                                       #
# python vcf_similarity_matrix.py --vcf <.vcf>                         #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################


def similarityMatrixHap(vcf):
	matResults = []
	with open(vcf, 'r') as fi, open(os.path.splitext(vcf)[0] + "_similairty_matrix.txt", 'w') as fii:
		data = [x.strip().split() for x in fi.readlines() if not x.startswith("##")]
		matResults.append("\t".join([""] + data[0][9:]))
		for i in range(9, len(data[0])):
			row = []
			row.append(data[0][i])
			for c in range(9, len(data[0])):
				total = 0
				match = 0
				for x in data[1:]:
					total += 1
					if x[i][0] == x[c][0]:
						match += 1
				row.append(str(int(match)/int(total)*100))
			matResults.append("\t".join(row))
		for x in matResults:
			fii.write(f"{x}\n")

def similarityMatrixDip(vcf):
	matResults = []
	with open(vcf, 'r') as fi, open(os.path.splitext(vcf)[0] + "_similarity_matrix.txt", 'w') as fii:
		data = [x.strip().split() for x in fi.readlines() if not x.startswith("##")]
		matResults.append("\t".join([""] + data[0][9:]))
		for i in range(9, len(data[0])):
			row = []
			row.append(data[0][i])
			for c in range(9, len(data[0])):
				total = 0
				match = 0
				for x in data[1:]:
					total += 1
					if x[i][0:3] == x[c][0:3]:
						match += 1
				row.append(str(int(match)/int(total)*100))
			matResults.append("\t".join(row))
		for x in matResults:
		 	fii.write(f"{x}\n")

def main():
	ap=argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		required = True,
		type = str,
		help = 'VCF file'
		)
	ap.add_argument(
		'--ploidy',
		required =True,
		type=str,
		help = 'Oranism ploidy'
		)
	parse = ap.parse_args()
	if parse.ploidy == "1":
		similarityMatrixHap(parse.vcf)
	else:
		similarityMatrixDip(parse.vcf)

if __name__ == "__main__":
	import argparse, os
	main()