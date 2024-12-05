#!/usr/bin/env python 

import sys, os 

with open(sys.argv[1], 'r') as file: # open(os.path.splitext(sys.argv[1])[0] + "_unwrapped.aln", 'w') as out:
	data = [x.strip() for x in file.readlines()] 
	seq = ""
	for x in data:
		if x.startswith(">"):
			if seq:
				print(seq)
				seq = ""
			print(x)
			continue
		else:
			seq+=x
	print(seq)