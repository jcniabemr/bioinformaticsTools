#!/usr/bin/env python 

########################################################################
# Script to caluclate contig lengths                                   #
#                                                                      # 
# Usage example:                                                       #
# python calculateContigLengths.py --assembly <assembly.fasta>         #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2024                                                                 #
########################################################################


def calculateLengths(g):
	contigFound = False
	contigLengthDict = {}
	seqLenCounter = 0 
	with open(g, 'r') as fi, open(os.path.splitext(g)[0] + '_contigLengths.txt', 'w') as out:
		data = [x.strip() for x in fi.readlines()]
		for x in data:
			if x.startswith(">"):
				if not contigFound:
					contigFound = True
					contigName = x[1:]
					continue
				if contigFound:
					contigLengthDict[contigName] = str(seqLenCounter)
					contigName = x[1:]
					seqLenCounter = 0
				continue
			seqLenCounter += len(x)
		contigLengthDict[contigName] = str(seqLenCounter)
		contigLengthDict["Total len"] = sum([int(i) for i in contigLengthDict.values()])
		for x,y in contigLengthDict.items():
			out.write(f"{x}\t{y}\n")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--assembly',
		type = str,
		required = True,
		help = 'Genome assembly'
		)
	parse = ap.parse_args()
	calculateLengths(parse.assembly)

if __name__ == '__main__':
	import os,argparse
	main()