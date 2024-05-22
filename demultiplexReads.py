#!/usr/bin/env python 

def demultiplex(reads):
	print(reads)
	


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--reads',
		type = str,
		required = True,
		help = 'reads .fastq infile'
		)
	demultiplex(ap.parse_args().reads)

if __name__ == '__main__':
	import argparse
	main()