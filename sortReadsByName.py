#!/usr/bin/env python 

########################################################################
# Script to sort reads in fastq files by name                          #
# Usage example:                                                       #
# python sortReadsByName.py --readF <_1.fastq> --readR <_2.fastq>      #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2024                                                                 #
########################################################################

def sortReads(a,b):
	with open(a, 'r') as fi, open(b, 'r') as fii, open(os.path.splitext(a)[0] + '_matchedLines.fastq' , 'w') as forwardOut, open(os.path.splitext(b)[0] + '_matchedLines.fastq', 'w') as reverseOut:
		forwardRead = [x.strip() for x in fi.readlines()]
		reverseRead = [x.strip() for x in fii.readlines()]
		seqHeaders = list(forwardRead[0::4])
		for x in seqHeaders:
			if x in reverseRead:
				fVal = [i for i, v in enumerate(forwardRead) if v == x]
				forwardOut.write('\n'.join(forwardRead[fVal[0]:fVal[0] + 4]) + '\n')
				forwardRead = forwardRead[:fVal[0]] + forwardRead[fVal[0] + 4:]
				rVal = [i for i, v in enumerate(reverseRead) if v == x]
				reverseOut.write('\n'.join(reverseRead[rVal[0]:rVal[0] + 4]) + '\n')
				reverseRead = reverseRead[:rVal[0]] + reverseRead[rVal[0] + 4:]

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--readF',
		type = str,
		required = True,
		help = 'Forward Read'
		)
	ap.add_argument(
		'--readR',
		type = str,
		required = True,
		help = 'Reverse Read'
		)
	parse = ap.parse_args()
	sortReads(parse.readF, parse.readR)

if __name__ == '__main__':
	import os,argparse
	main()