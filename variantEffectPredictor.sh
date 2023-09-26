#!/usr/bin/env bash 
#SBATCH -J VEP
#SBATCH -p long
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10

###################################################################################
# Script to predict the effect of genetic variants and produce a gff              #
#  of effected genes containing variants with predicted effect                    #
#                                                                                 #
# Usage example:                                                                  #
# python variantEffectorPredictor.sh --gff <.gff> --fasta <.fasta> --vcf <.vcf>   #
#                                                                                 #
# Written by John Connell                                                         #
# john.connell@niab.com                                                           #
# NIAB                                                                            #
# 2023                                                                            #
###################################################################################

#vcfFile=$1 gffFile=$2 assembly=$3

cp ${1} ${TMPDIR}/VCF
cd ${TMPDIR}

alleles=$(
python << nestedVariants
variantDict={}
with open("VCF", 'r') as fi:
	data = [x.strip().split() for x in fi.readlines() if not x.startswith("##")]
	for i in range(9, len(data[0])):
		strain = data[0][i]
		if strain not in variantDict:
			variantDict[strain] = {}
		for x in data[1:]:
			if x[i][0] == "0":
				variantDict[strain][x[1]] = x[3]
			elif  x[i][0] == "1":
				variantDict[strain][x[1]] = x[4]
			else:
				pass
for x,y in variantDict.items():
	print(x,y)
nestedVariants
)


while read -r x; do 
	echo $x
done <<< ${alleles}



# #def writeGFF():


# def processVariants(vcf):
# 	header=[x for x in vcf if x.startswith("#")]
# 	#for x in vcf:
# 	#	samples = x[9:]





# def readFiles(gff,fasta,vcf):
# 	with open(gff, 'r') as fi:
# 		gffData = [x.strip().split() for x in fi.readlines() if not x.startswith("#")]
# 	with open(fasta, 'r') as fii:
# 		genome = [x.strip() for x in fii.readlines() if not x.startswith("#")]
# 	with open(vcf, 'r') as fiii:
# 		vcfData = [x.strip().split() for x in fiii.readlines() if not x.startswith("##")]
# 	return gffData,genome,vcfData

# def main():
# 	ap=argparse.ArgumentParser()
# 	ap.add_argument(
# 		'--gff',
# 		required = True,
# 		type = str,
# 		help = 'gff file'
# 		)
# 	ap.add_argument(
# 		'--fasta',
# 		required = True,
# 		type = str,
# 		help = 'Reference genome'
# 		)
# 	ap.add_argument(
# 		'--vcf',
# 		required = True,
# 		type = str,
# 		help = 'Variant call file'
# 		)
# 	parse=ap.parse_args()
# 	processVariants(
# 	readFiles(
# 	parse.gff,
# 	parse.fasta,
# 	parse.vcf
# 	))


# if __name__=="__main__":
# 	import argparse
# 	main()
