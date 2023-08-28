#!/usr/bin/env python 

########################################################################
# Script to combine one or more gff files into a single file,          # 
# remove duplicated features, take the longest isoform,                #
# sort and rename features uniformly                                   #
#                                                                      #
# Usage example:                                                       #
# python sortGFFs.py --gff <.gff1 .gff2 .gff3 .gffX>                   #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

####Open, extract and order gff features. 
def parseFiles(gff):
	catData = []
	for x in gff:
		with open(x,'r') as f:
			data = [x.strip().split() for x in f if not x.startswith("#")]
			catData.extend(data)
	sortedCatData = sorted(
		catData,
		key = lambda i: (int(i[0].split("_")[1]),int(i[3])))
	return sortedCatData

####Extract gene groups ascending. 
def createGeneGroups(i):
	from collections import OrderedDict
	geneGroup = OrderedDict()
	for x in i:
		if x[2] == "gene":
			gene="_".join([x[2],x[3],x[4]])
			geneGroup[gene] = []
			for y in i:
				if int(y[3]) >= int(x[3]) and int(y[4]) <= int(x[4]) and y[0] == x[0]:
					geneGroup[gene].append(y)
			print(geneGroup[gene])
	return geneGroup

####Criteria to order features within gene groups. 
def orderCriteria(feature):
	import re
	transcriptNumber = re.match(r".*\.t(\d+)",feature[-1])
	if transcriptNumber:
		tn = int(transcriptNumber.group(1))
	else:
		tn = 1
	if feature[2] == "gene":
		return (0,1)
	elif feature[2] == "mRNA":
		return (tn,2)
	else:
		return (tn,3,feature[3])

####Order features.
def orderFeatures(geneDict):
	from collections import OrderedDict
	orderedGeneDict = OrderedDict()
	for gene, feature in geneDict.items():
		print(geneDict[gene])
		orderedFeatures = sorted(
			feature, 
			key = orderCriteria)
		orderedGeneDict[gene] = orderedFeatures
	return orderedGeneDict

####Remove dups and take longes isoform.
def RemoveDupsLongestIsoform(orderedDict):
	processedFeature = set()
	gff = []
	for x in orderedDict.values():
		for y in x:
			if tuple(y) not in processedFeature:
				gff.append(y)
				processedFeature.add(tuple(y))
	return gff

def locateFeatures(gff):
	renamedGff = []
	geneCounter = 0
	transcriptCounter = 0 
	exonCounter = 0 
	for x in gff:
		if x[2] == "gene": 
			geneCounter += 1
			transcriptCounter = 0
			exonCounter = 0
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter])))
			renamedGff.append("\t".join(x[0:9]))
			continue
		elif x[2] == "mRNA":
			transcriptCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,";Parent=g",geneCounter])))
			renamedGff.append("\t".join(x[0:9]))
			continue 
		elif x[2] == "stop_codon" or x[2] == "start_codon" or x[2] == "intron":
			x[8] = x[8].replace(x[8],"".join(map(str,["Parent=g",geneCounter,".t",transcriptCounter])))
			renamedGff.append("\t".join(x[0:9])) 
			continue 
		elif x[2] == "exon":
			exonCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".exon",exonCounter,";Parent=g",geneCounter,".t",transcriptCounter])))
			renamedGff.append("\t".join(x[0:9]))
			continue 
		elif x[2] == "CDS":
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".cds;Parent=g",geneCounter,".t",transcriptCounter])))
			renamedGff.append("\t".join(x[0:9]))
			continue 
		else:
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,";Parent=g",geneCounter])))
			renamedGff.append("\t".join(x[0:9]))
	return renamedGff

####Write results. 
def writeRes(file):
	with open("Sorted_renamed_genes.gff",'w') as outfile:
		for x in file:
			outfile.write(f"{x}\n")

####Run script
def main():
	import argparse 
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--gff',
		type=str,
		required=True,
		nargs='+',
		help='Space seperated list of gff files to process'
		)
	parse = ap.parse_args()
	writeRes(
	locateFeatures(
	RemoveDupsLongestIsoform(
	orderFeatures(
	createGeneGroups(
	parseFiles(
	parse.gff)
	)))))

if __name__=="__main__":
	main()