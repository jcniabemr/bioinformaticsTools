#!/usr/bin/env python 

########################################################################
# Script to rename contigs in a genome assembly fasta file             # 
#                                                                      #
# Usage example:                                                       #
# python renameContigs.py --assembly <.fasta> --contigName <contig>    #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

def writeRes(i,c):
	newName = os.path.splitext(c)[0] + "_renamed.fa"
	with open(newName,'w') as out:
		for x in i:
			out.write(f"{x}\n")

def renameContigs(assembly,name):
	results = []
	contigCouter = 0
	with open(assembly,'r') as fi:
		data=[x.strip() for x in fi.readlines() if not x.startswith("#")]
		for x in data:
			if x.startswith(">"):
				contigCouter += 1
				results.append(x.replace(x,"".join(map(str,[">",name,"_",contigCouter]))))
			else:
				results.append(x)
	writeRes(results,assembly)

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--assembly',
		required = True,
		type = str,
		help = 'Genome assembly file'
		)
	ap.add_argument(
		'--contigName',
		required = True,
		type = str,
		help = 'New contig prefix'
		)
	parse = ap.parse_args()
	renameContigs(parse.assembly,parse.contigName)

if __name__=="__main__":
	import argparse,os
	main()