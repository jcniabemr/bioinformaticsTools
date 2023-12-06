#!/usr/bin/env python 

###########################################################################
# Script to estimate the sequencing coverage of a genome                  #
# 				                                                          # 
# Options: -i read.fastq.gz -s genome size Mb	                          #
# Usage example:                                                          #
# python calculateSequencingCoverage.py -i <.fasta1> <.fastaX> -s <30>    #
#                                                                         #
# Written by John Connell                                                 #
# john.connell@niab.com                                                   #
# NIAB                                                                    #
# 2023                                                                    #
###########################################################################

def countBases(f):
	rc = 0
	bc = 0
	for x in f:
		x = x.strip()
		if (rc - 1) % 4 == 0:
			bc += len(x)
		rc += 1
	return bc

def calculateCoverage(f,s):
	c = 0
	for i in f:
		if os.path.splitext(i)[1] == ".gz":
			with gzip.open(i, 'rt') as file:
				c += countBases(file)
		else:
			with open(i, 'r') as file:
				c += countBases(file)	
	print(str(round(int(c)/(int(1000000)*s))) + "x coverage")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'-i',
		type = str,
		required = True,
		nargs = '+',
		help = 'Space sep list of fastq files'
		)
	ap.add_argument(
		'-s',
		type = int,
		required = True,
		help = 'Genome size in Mb'
		)
	parse = ap.parse_args()
	calculateCoverage(parse.i,parse.s)

if __name__ == '__main__':
	import os,argparse,gzip
	main()