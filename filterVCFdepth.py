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

def autoFilterVCF(vcf):
	with open(vcf, 'r') as file, open(os.path.splitext(vcf)[0] + "_AutoDepthFilter.vcf", 'w') as out:
		modeDict = {}
		fi = list(file)
		out.write("".join(map(str, [x for x in fi if x.startswith("##")])))
		data = [x.strip().split() for x in fi if not x.startswith('##')]
		for i in range(9, len(data[0])):
			strain = data[0][i]
			depthList = []
			modeDict[strain] = []
			for x in data[1:]:
				if x[i].split(":")[2] != ".":
					depthList.append(x[i].split(":")[2])
			print(strain,i,len(depthList))
			#modeDict[strain].append(float(mode(depthList))*1.99)

	# for x,y in modeDict.items():
	# 	print(x,y)



def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		required = True,
		help  = 'vcf file'
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
		action='store_true', 
		required = False,
		help = 'Auto filter all samples by depth mode +/- mode *1.8'
		)
	parse = ap.parse_args()
	if parse.auto:
		autoFilterVCF(parse.vcf)
	else:
		filterVCF(parse.vcf,parse.minDP,parse.maxDP)

if __name__ == '__main__':
	import os,argparse
	from statistics import mode 
	main()
