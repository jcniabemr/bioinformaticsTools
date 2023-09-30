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


####Check and assign input args 
while [[ $# -gt 0 ]]; do
  case ${1} in
    --vcf)
      vcfFile=${2}
      shift 
      shift 
      ;;
    --gff)
      gffFile=${2}
      shift 
      shift 
      ;;
    --fasta)
      fastaFile=${2}
      shift 
      shift 
      ;;
    --ploidy)
	  ploidy=${2}
	  shift
	  shift
	  ;;
    *)
      echo "Unknown option: ${1}"
      exit 1
      ;;
  esac
done

if [ -z ${gffFile} ] || [ -z ${vcfFile} ] || [ -z ${fastaFile} ] || [ -z ${ploidy} ]; then
	echo "Error: Please provide GFF, VCF, and FASTA files using --gff, --vcf, --fasta, --ploidy flags."
	exit 1
else
	echo "Using GFF file: ${gffFile}"
	echo "Using VCF file: ${vcfFile}"
	echo "Using FASTA file: ${fastaFile}"
	echo "Ploidy is ${ploidy}"
fi

####Copy required infiles 
cp ${vcfFile} ${TMPDIR}/VCF
cp ${gffFile} ${TMPDIR}/GFF
cp ${fastaFile} ${TMPDIR}/FASTA
cd ${TMPDIR}

####Funcion process diploid variants
processDipVariants() { 
python  << nestedVariants
variantDict={}
with open("$1", 'r') as fi:
	data = [x.strip().split() for x in fi.readlines() if not x.startswith("##")]
	for i in range(9, len(data[0])):
		strain = data[0][i]
		if strain not in variantDict:
			variantDict[strain] = {}
		for x in data[1:]:
			if x[i][0] == "1" and x[i][2] == "1":
				variantDict[strain]["_".join(x[0:2] + x[3:5])] = [x[4] + "/" + x[4]]
			elif  x[i][0] == "0" and x[i][2] == "1":
				variantDict[strain]["_".join(x[0:2] + x[3:5])] = [x[3] + "/" + x[4]]
			else:
				pass
	for x,y in variantDict.items():
		for y,z in y.items():
			print(x,y,z)
nestedVariants
}

####Determine run type 
if [ ${ploidy} -eq 1 ]; then
	alleles=$(processHapVariants VCF)
elif [ ${ploidy} -eq 2 ]; then 
	alleles=$(processDipVariants VCF)
else
	echo "Ploidy ${ploidy} is not aloud, please set a ploidy of either 1 or 2"
	exit 1
fi


while read -r x; do 
	echo $x
	# strain=$(echo ${x} | awk '{print $1}')
	# contig=$(echo ${x} | awk '{split($2,a,"_");{print a[1]"_"a[2]}}')
	# position=$(echo ${x} | awk '{split($2,a,"_");{print a[3]}}')
	# allele=$(echo $x | awk '{print $3}')
	# echo ${strain} ${contig} ${position} ${allele}
done <<< ${alleles}



