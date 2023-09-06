#!/usr/bin/env python 

########################################################################
# Script to produce tiles of SNP density from VCF file                 #
# Usage example:                                                       #
# python createSNPtiles.py --vcf <VCF file>                            #
#                                                                      # 
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

#def writeRes():


def openSort(vcf):
	with open(vcf, 'r') as fi:
		vcfData=[x.strip().split() for x in fi.readlines() if not x.startswith("#")]
		sortedData=sorted(
			vcfData, 
			key = lambda i: (int(i[0].split("_")[1]), int(i[1])))
	return sortedData

def createWindows(vcf,window):
	sortedData=openSort(vcf)
	contigTracker = ""
	lenTracker = 0
	snpCount = 0
	for x in sortedData:
		if contigTracker == "":
			contigTracker = x[0]
			lenTracker = int(x[1])
			snpCount += 1
			windowTracker = 0
			continue
		if x[0] == contigTracker:
			if windowTracker + int(x[1]) - int(lenTracker) < int(windowSize) -1:
				windowTracker += int(x[1]) - int(lenTracker)
				snpCount += 1
				lenTracker = int(x[1])
				continue
			else:
				results.write("\t".join(map(str,[contigTracker,snpCount,"\n"])))
				snpCount = 1
				windowTracker = 0
				lenTracker = int(x[1])
				continue
		if x[0] != "" and x[0] != contigTracker:
			results.write("\t".join(map(str,[contigTracker,snpCount,"\n"])))
			contigTracker=x[0]
			snpCount=1
			windowTracker=0
			lenTracker=int(x[1])
			continue
	results.write("\t".join(map(str,[contigTracker,snpCount,"\n"])))
	results.close()









def main():
	ap=argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type=str,
		required=True,
		help='VCF file'
		)
	ap.add_argument(
		'--windowSize',
		type=int,
		required=True,
		help='Window size to partion by'
		)
	parse=ap.parse_args()
	createWindows(parse.vcf,parse.windowSize)

if __name__=="__main__":
	import argparse
	main()