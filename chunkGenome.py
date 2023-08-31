#!/usr/bin/env python 

########################################################################
# Script to cut genome or proteome into chunks for tools               #
# such as SignalP where a maximum file size is limit exists            #
#                                                                      #
# Usage example:                                                       #
# python chunkGenome.py --infile <.fasta> --chunkSize <1000>           #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

def writeFile(file,counter,chunkFile):
	filename="".join([file.split(".")[0],"_","partition","_",str(counter),".",file.split(".")[1]])
	with open(filename,'w') as out:
		for x in chunkFile:
			out.write(f"{x}\n")

####Open extract info.
def parseInfo(file,chunkSize):
	chunkFile = []
	contigFound = False
	contigCounter = 0
	fileCounter = 0
	seq = ""
	with open(file,'r') as f:
		data=[x.strip() for x in f.readlines() if not x.startswith("#")]
	for x in data:
		if x.startswith(">"):
			if contigFound:
				chunkFile.append("\n".join([seq[i:i+60] for i in range(0,len(seq),60)])+"\n")
				seq = ""
				contigFound = False 
			if int(contigCounter) == chunkSize:
				fileCounter += 1
				writeFile(file,(int(contigCounter)*int(fileCounter)),chunkFile)
				chunkFile = []
				contigCounter = 0
			contigFound = True
			contigCounter += 1
			chunkFile.append(x)
			continue 
		if contigFound:
			seq += x
	fileCounter += 1
	chunkFile.append("\n".join([seq[i:i+60] for i in range(0,len(seq),60)])+"\n")
	writeFile(file,(chunkSize * int(fileCounter)),chunkFile)

####Main func.
def main():
	import argparse 
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--infile',
		type=str,
		required=True,
		help='Proteome or genome assembly'
		)
	ap.add_argument(
		'--chunkSize',
		type=int,
		required=True,
		help='Number of elements to seperate by'
		)
	parse = ap.parse_args()
	parseInfo(parse.infile,parse.chunkSize)
	

####Run prog.
if __name__=="__main__":
	main()