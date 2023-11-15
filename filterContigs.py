#!/usr/bin/env python 

######################################################################################
# Script to filter contigs in a genome assembly by contig length,                    #
#  rename contigs and wrap sequences                                                 #
#                                                                                    #
# Usage example:                                                                     #
# python filterContifs.py --fasta <.fasta> --minLengh <500>  --contigName <contig>   #
#                                                                                    #
# Written by John Connell                                                            #
# john.connell@niab.com                                                              #
# NIAB                                                                               #
# 2023                                                                               #
######################################################################################

def filterAssembly(fasta,size,name):
	contigDict = {}
	i = 0
	contigFound = False
	with open(fasta, 'r') as fi, open(os.path.splitext(fasta)[0] + "".join(["_filtered_contigs_min_",size,"bp.fa"]), 'w') as fii:
		for x in fi:
			x = x.strip()
			if x.startswith(">"):
				if contigFound:
					if len(seq) > int(size):
						contigDict["".join([">",name,"_",str(i)])] = seq
						i+=1
					seq = ""
					continue
				if not contigFound:
					contigFound = True
					i+=1
					seq = ""
					continue
			seq += x 
		for x,y in contigDict.items():
			fii.write(x + "\n")
			for k in range(0, len(y), 60):
				fii.write(y[k:k+60] + "\n")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--fasta',
		required = True,
		type = str,
		help = "Genome assembly"
		)
	ap.add_argument(
		'--minLength',
		required = True,
		type = str,
		help = "Minimum contig length"
		)
	ap.add_argument(
		'--contigName',
		required = True,
		type = str,
		help = "New contig name"
		)
	parse = ap.parse_args()
	filterAssembly(
		parse.fasta, 
		parse.minLength,
		parse.contigName
		)

if __name__ == "__main__":
	import os,argparse
	main()