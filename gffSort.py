#!/usr/bin/env python 

########################################################################
# Script to combine one or more gff files into a single file,          # 
# remove duplicated features, sort and rename features uniformly.      #
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
			data = [x.strip().split() for x in f.readlines() if not x.startswith("#")]
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
				if y[0] == x[0] and int(y[3]) >= int(x[3]) and int(y[4]) <= int(x[4]) and x[6] == y[6]:
					geneGroup[gene].append(y)
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
		return (tn,3,int(feature[3]))

####Order features.
def orderFeatures(geneDict):
	from collections import OrderedDict
	orderedGeneDict = OrderedDict()
	for gene, feature in geneDict.items():
		orderedFeatures = sorted(
			feature, 
			key = orderCriteria)
		orderedGeneDict[gene] = orderedFeatures
	return orderedGeneDict

####Remove duplicated features.
def RemoveDups(orderedDict):
	processedFeature = set()
	gff = []
	for x in orderedDict.values():
		for y in x:
			if tuple(y) not in processedFeature:
				gff.append(y)
				processedFeature.add(tuple(y))
	return gff

####Rename all features according to gene and feature count. 
def locateFeatures(gff):
	renamedGff = []
	geneCounter = 0
	for x in gff:
		if x[2] == "gene": 
			geneCounter += 1
			transcriptCounter = 0 
			exonCounter = 0 
			cdsCounter = 0
			startCounter = 0 
			stopCounter = 0
			intronCounter = 0
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,";"])))
			renamedGff.append("\t".join(x[0:9]))
			continue
		elif x[2] == "mRNA":
			transcriptCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,";Parent=g",geneCounter,";"])))
			renamedGff.append("\t".join(x[0:9]))
		elif x[2] == "tRNA":
			transcriptCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,";Parent=g",geneCounter,";"])))
			renamedGff.append("\t".join(x[0:9]))
			continue 
		elif x[2] == "start_codon":
			startCounter += 1 
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".start",startCounter,";Parent=g",geneCounter,".t",transcriptCounter,";"])))
			renamedGff.append("\t".join(x[0:9])) 
			continue 
		elif x[2] == "stop_codon": 
			stopCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".stop",stopCounter,";Parent=g",geneCounter,".t",transcriptCounter,";"])))
			renamedGff.append("\t".join(x[0:9])) 
			continue 
		elif x[2] == "intron":
			intronCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".intron",intronCounter,";Parent=g",geneCounter,".t",transcriptCounter,";"])))
			renamedGff.append("\t".join(x[0:9])) 
			continue 
		elif x[2] == "exon":
			exonCounter += 1
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".exon",exonCounter,";Parent=g",geneCounter,".t",transcriptCounter,";"])))
			renamedGff.append("\t".join(x[0:9]))
			continue 
		elif x[2] == "CDS":
			cdsCounter += 1 
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".CDS",cdsCounter,";Parent=g",geneCounter,".t",transcriptCounter,";"])))
			renamedGff.append("\t".join(x[0:9]))
			continue 
		else:
			x[8] = x[8].replace(x[8],"".join(map(str,["ID=g",geneCounter,".t",transcriptCounter,".",x[2],";Parent=g",geneCounter,".t",transcriptCounter,";"])))
			renamedGff.append("\t".join(x[0:9]))
	return renamedGff

####Write results. 
def writeRes(file):
	with open("Sorted_renamed_genes.gff",'w') as outfile:
		for x in file:
			outfile.write(f"{x}\n")

####Main func.
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
	RemoveDups(
	orderFeatures(
	createGeneGroups(
	parseFiles(
	parse.gff)
	)))))

####Run prog.
if __name__=="__main__":
	main()