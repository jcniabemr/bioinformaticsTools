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

def writeRes(i,c):
	resultFile = os.path.splitext(c)[0] + ".SNP_windows.txt"
	with open(resultFile,'w') as out:
		out.write("\t".join(map(str,["Contig","Number of SNPs","Window number"]))+"\n\n")
		for x in i:
			out.write(f"{x}\n") 

def openSort(vcf):
	with open(vcf, 'r') as fi:
		vcfData = [x.strip().split() for x in fi.readlines() if not x.startswith("#")]
		sortedData = sorted(
			vcfData, 
			key = lambda i: (int(i[0].split("_")[1]), int(i[1])))
	return sortedData

def createWindows(vcf,windowSize):
	sortedData=openSort(vcf)
	contigTracker = ""
	lenTracker = 0
	snpCount = 0
	windowNumber = 0
	windowResults = []
	for x in sortedData:
		if contigTracker == "":
			contigTracker = x[0]
			lenTracker = int(x[1])
			snpCount += 1
			windowNumber += 1
			windowTracker = 0
			continue
		if x[0] == contigTracker:
			if windowTracker + int(x[1]) - int(lenTracker) < int(windowSize) -1:
				windowTracker += int(x[1]) - int(lenTracker)
				snpCount += 1
				lenTracker = int(x[1])
				continue
			else:
				windowResults.append("\t".join(map(str,[contigTracker,snpCount,windowNumber])))
				snpCount = 1
				windowNumber += 1
				windowTracker = 0
				lenTracker = int(x[1])
				continue
		if x[0] != "" and x[0] != contigTracker:
			windowResults.append("\t".join(map(str,[contigTracker,snpCount,windowNumber])))
			contigTracker = x[0]
			snpCount = 1
			windowNumber = 1
			windowTracker = 0
			lenTracker = int(x[1])
			continue
	windowResults.append("\t".join(map(str,[contigTracker,snpCount,windowNumber])))
	writeRes(windowResults,vcf)

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
	import argparse,os
	main()