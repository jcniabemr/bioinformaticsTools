#!/usr/bin/env python 

def plot(data, plotType, file, name):
	savePath = os.path.join(os.path.expanduser('~'), 'wheatPlots', plotType, file, plotType + '_of_' + name + '.jpeg')
	os.makedirs(os.path.dirname(savePath), exist_ok = True)
	sns.histplot(data, bins = 2000)
	plt.title(f'{plotType} of {name}')
	plt.ylabel('Frequency')
	plt.xlabel(f'{plotType} value')
	plt.savefig(savePath, format = 'jpeg')
	plt.clf()

def extractInfo(vcf):
	with open(vcf, 'r') as fi:
		data = [x.strip().split() for x in fi.readlines() if not x.startswith('##')]
		#val = data[1][8].split(':').index('DP')
		val = [n for n,c in enumerate(data[1][8].split(":")) if c == "DP"]
		for i,c in enumerate(data[0][9:], start = 9):
			depthList = []
			for x in data[1:]:
				depth = x[i].split(':')[val[0]]
				if depth != '.':
					depthList.append(depth)
			plot(depthList, "Depth", os.path.splitext(vcf)[0], c)
		for x in data[1:]:
			qualList = []
			qualList.append(x[5])
			plot(qualList, "Qual", os.path.splitext(vcf)[0], "overallQual")

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument(
		'--vcf',
		type = str,
		required = True,
		help = 'vcf file',
	)
	parse = ap.parse_args()
	extractInfo(parse.vcf)

if __name__ == '__main__':
	import os,argparse
	import seaborn as sns 
	import matplotlib.pyplot as plt
	main()