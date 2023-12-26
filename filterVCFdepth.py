#!/usr/bin/env python 

###############################################################################
# Script to filter a VCF per sample by depth                                  #
# 				                                                              # 
# Usage example:                                                              #
# python filterVCFdepth.py --vcf file.vcf --minDP <integer> --maxDP <integer> #   
#                                                                             #
# Written by John Connell                                                     #
# john.connell@niab.com                                                       #
# NIAB                                                                        #
# 2023                                                                        #
###############################################################################

def autoFilterQual(vcf):
	with open(vcf, 'r') as file, open(os.path.splitext(vcf)[0] + "_AutoQualFilter.vcf", 'w') as out:
		qualList = []
		fi = list(file)
		out.write("".join(map(str, [x for x in fi if x.startswith("#")])))
		data = [x.strip().split() for x in fi if not x.startswith('#')]
		for x in data:
			qualList.append(float(x[5]))
		percentile=np.percentile(qualList, 25)
		for x in data:
			if float(x[5]) > float(percentile):
				out.write("\t".join(x) + "\n")

def autoFilterDepth(vcf):
	with open(vcf, 'r') as file, open(os.path.splitext(vcf)[0] + "_AutoDepthFilter.vcf", 'w') as out:
		filterDict = {}
		fi = list(file)
		out.write("".join(map(str, [x for x in fi if x.startswith("#")])))
		data = [x.strip().split() for x in fi if not x.startswith('##')]
		for i in range(9, len(data[0])):
			strain = data[0][i]
			depthList = []
			for x in data[1:]:
				if x[i].split(":")[2] != ".":
					depthList.append(int(x[i].split(":")[2]))
	#			print("strain is ",strain," val ",i," depth list has ",len(depthList)," values", "95 percentile is ",np.percentile(depthList, 95)," mode is ", mode(depthList))
			filterDict[strain] = float(np.percentile(depthList, 95))
		for x in data[1:]:
			row = []
			for i in range(0,9):
				row.append(x[i])
			for i in range(9, len(data[0])):
				strain = data[0][i]
				if x[i].split(":")[2] == ".":
					row.append(x[i])
					continue 
				if int(filterDict[strain]) > int(x[i].split(":")[2]) > int(3):
					row.append(x[i])
					continue
				else:
					row.append(":".join(["."] + x[i].split(":")[1:]))
			out.write("\t".join(row) + "\n")


def filterVCF(vcf, minDP, maxDP):
	with open(vcf, 'r') as file, open(os.path.splitext(vcf)[0] + "_depthFilter.vcf", 'w') as out:
		fi = list(file)
		out.write("".join(map(str, [x for x in fi if x.startswith("#")])))
		data = [x.strip().split() for x in fi if not x.startswith("##")]
		for x in data[1:]:
			row = []
			for i in range(0,9):
				row.append(x[i])
			for i in range(9, len(data[0])):
				if x[i].split(":")[2] == ".":
					row.append(x[i])
					continue 
				if maxDP > int(x[i].split(":")[2]) > minDP:
					row.append(x[i])
					continue 
				else:
					row.append(":".join(["."] + x[i].split(":")[1:]))
			out.write("\t".join(row) + "\n")


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		required = True,
		help  = 'vcf file'
		)
	ap.add_argument(
		'--qual',
		action = 'store_true',
		required = False, 
		help = 'include to perform quality filter'
		)
	ap.add_argument(
		'--depth',
		action = 'store_true',
		required = False, 
		help = 'include to perform depth filter'
		)
	ap.add_argument(
		'--minDP',
		type = int,
		required = False,
		help = 'min depth'
		)
	ap.add_argument(
		'--maxDP',
		type = int,
		required = False,
		help = 'max depth'
		)
	ap.add_argument(
		'--auto',
		action = 'store_true', 
		required = False,
		help = 'Auto filter all samples by depth mode +/- mode *1.8'
		)
	parse = ap.parse_args()
	if parse.qual:
		autoFilterQual(parse.vcf)
	elif parse.depth:
		if parse.auto:
			autoFilterDepth(parse.vcf)
		else:
			filterDepth(parse.vcf,parse.minDP,parse.maxDP)

if __name__ == '__main__':
	import os,argparse
	from statistics import mode 
	import numpy as np 
	main()
