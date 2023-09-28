#!/usr/bin/env bash 
#SBATCH -J VEP
#SBATCH -p long
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10

###################################################################################
# Script to predict the effect of genetic variants from a VCF                     #
#                                                                                 #
# Usage example:                                                                  #
# python variantEffectorPredictor.sh --gff <.gff> --fasta <.fasta> --vcf <.vcf>   #
#                                                                                 #
# Written by John Connell                                                         #
# john.connell@niab.com                                                           #
# NIAB                                                                            #
# 2023                                                                            #
###################################################################################


# ####Check and assign args 
# while [[ $# -gt 0 ]]; do
#   case ${1} in
#     --vcf)
#       vcfFile=${2}
#       shift 
#       shift 
#       ;;
#     --gff)
#       gffFile=${2}
#       shift 
#       shift 
#       ;;
#     --fasta)
#       fastaFile=${2}
#       shift 
#       shift 
#       ;;
#     *)
#       echo "Unknown option: ${1}"
#       exit 1
#       ;;
#   esac
# done

# if [ -z ${gffFile} ] || [ -z ${vcfFile} ] || [ -z ${fastaFile} ]; then
# 	echo "Error: Please provide GFF, VCF, and FASTA files using --gff, --vcf, --fasta flags."
# 	exit 1
# else
# 	echo "Using GFF file: ${gffFile}"
# 	echo "Using VCF file: ${vcfFile}"
# 	echo "Using FASTA file: ${fastaFile}"
# fi




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
			if x[i][0] == "1" or x[i][2] == "1":
				variantDict[strain]["_".join([x[0],x[1]])] = x[4]
#			elif  x[i][0] == "0" and x[i][2] == "1":
#				variantDict[strain]["_".join([x[0],x[1]])] = x[4]
			else:
				pass
for x,y in variantDict.items():
	for y,z in y.items():
		print(x,y,z)
nestedVariants
)

while read -r x; do 
	strain=$(echo ${x} | awk '{print $1}')
	contig=$(echo ${x} | awk '{split($2,a,"_");{print a[1]"_"a[2]}}')
	position=$(echo ${x} | awk '{split($2,a,"_");{print a[3]}}')
	allele=$(echo $x | awk '{print $3}')
	#echo ${strain} ${contig} ${position} ${allele}

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
