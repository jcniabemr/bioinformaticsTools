#!/usr/bin/env python

#######################################################################################
# Script to cut a genome into fragments that contains "N" characters                  #
#                                                                                     # 
# Options: --assembly --step --window                                                 #
#                                                                                     #
# Written by John Connell                                                             #
# john.connell@niab.com                                                               #
# NIAB                                                                                #
# 2024                                                                                #
#######################################################################################

def cutGenome(assembly,step,window):
	found = False 
	written = False
	header = ""
	with open(assembly, 'r') as file, open(os.path.splitext(assembly)[0] + "_split_into_" + str(window)+ "_at_" + str(step) + "_steps.fasta", 'w') as out:
		for x in file:
			x = x.strip()
			if x.startswith(">"):
				if found:
					for n in range(0, len(seq), step):
						fragmentCounter += 1 
						newEnd = min(n + window, len(seq))
						chunk = seq[n:newEnd]
						out.write(f"{header}_{n}_{newEnd}_{window}_{step}_fragment{fragmentCounter}\n") 
						out.write("\n".join([chunk[o:o + 60] for o in range(0, len(chunk), 60)]) + "\n")
				found = True
				seq = ""
				header = x
				fragmentCounter = 0
				continue
			for i in x:
				if i not in ['N', 'n']:
					seq += i.upper()
					written = False 
				else:
					if not written:
						for n in range(0, len(seq), step):
							fragmentCounter += 1 
							newEnd = min(n + window, len(seq))
							chunk = seq[n:newEnd]
							out.write(f"{header}_{n}_{newEnd}_{window}_{step}_fragment{fragmentCounter}\n") 
							out.write("\n".join([chunk[o:o + 60] for o in range(0, len(chunk), 60)]) + "\n")
						seq = ""
						written = True
					else:
						continue 
		for n in range(0, len(seq), step):
			fragmentCounter += 1 
			newEnd = min(n + window, len(seq))
			chunk = seq[n:newEnd]
			out.write(f"{header}_{n}_{newEnd}_{window}_{step}_fragment{fragmentCounter}\n") 
			out.write("\n".join([chunk[o:o + 60] for o in range(0, len(chunk), 60)]) + "\n")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--assembly',
		type = str,
		required = True,
		help = 'Assembly to be cut'
	)
	ap.add_argument(
		'--step',
		type = int,
		required = True,
		help = 'Step length for cutting'
	)
	ap.add_argument(
		'--window',
		type = int,
		required = True,
		help = 'Length of read'
	)
	parse=ap.parse_args()
	cutGenome(
		parse.assembly,
		parse.step,
		parse.window
	)

if __name__ == '__main__':
	import os,argparse
	main()
