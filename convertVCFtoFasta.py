#!/usr/bin/env python 

########################################################################
# Script to produce sequence from a VCF, can be used for ref guided    #
# seq using BCFtools output all sites                                  #
#                                                                      # 
# Usage example:                                                       #
# python convertVCFtoFasta.py --vcf <file.vcf>                         #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2024                                                                 #
########################################################################

def createFasta(vcf):
	start = False
	with open(vcf, 'r') as fi:
		data = [x.strip().split() for x in fi.readlines() if not x.startswith("##")]
		for i,c in enumerate(data[0][9:], start = 9):
			seq = ""
			seqCounter = 1
			for x in data[1:]:
				if not start:
					start = True
					if int(x[1]) != 1:
						seq += "N" * int(x[1])
						seqCounter = int(x[1]) 
				if start: 
					if int(x[1]) == seqCounter:
						seqCounter += 1
					else:
						seq += "N" * (int(x[1]) - (seqCounter -1))
						seqCounter = int(x[1]) + 1
					if x[i][0:3] == "0/0":
						seq += x[3]
					else:
						seq += x[4]
			with open(os.path.splitext(vcf)[0] + "_" + c + "_fasta.fa", 'w') as out:
				out.write(">contig_1" + "\n")
				out.write("\n".join([seq[k:k + 60] for k in range(0, len(seq), 60)]))
			break 

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		required = True,
		help = 'input vcf with all sites called'
	)
	parse = ap.parse_args()
	createFasta(parse.vcf)

if __name__ == '__main__':
	import os,argparse
	main()