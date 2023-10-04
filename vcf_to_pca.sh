#!/usr/bin/env bash 

########################################################################
# Script to produce similarity matrix of similair sites for all        #
# samples in a VCF file                                                #
#                                                                      #
# Usage example:                                                       #
# bash vcf_to_pca.sh --vcf <.vcf>                                      #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################

createVCF(){
Rscript - <<PCA

install.packages("adegenet")
install.packages("ggplot2")
library(adegenet)
library(ggplot2)

vcf_data <- read.vcf("$1", verbose = FALSE)

genotypes <- genind2genpop(genotype(vcf_data))
pca_result <- dudi.pca(genotypes, scannf = FALSE, nf = 2)

pca_df <- as.data.frame(pca_result$li)
colnames(pca_df) <- c("PC1", "PC2")

ggplot(pca_df, aes(x = PC1, y = PC2)) +
geom_point() +
labs(title = "PCA Plot", x = "Principal Component 1", y = "Principal Component 2")
ggsave("pca_plot.png", plot = last_plot(), width = 6, height = 4, dpi = 300)

PCA
}


while [[ $# -gt 0 ]]; do
	case ${1} in
    	--vcf)
    		vcfFile=${2}
    	  	shift 
      		shift 
    		;;
    	*)
      		echo "Unknown option: ${1}"
      		exit 1
      		;;
  	esac
done

if [ -z ${vcfFile} ]; then
	echo "Error: Please provide a ploidy value and a GFF, VCF, and FASTA file using arguments: --gff, --vcf, --fasta, --ploidy flags."
	exit 1
else
	echo "Using VCF file: ${vcfFile}"
	createVCF ${vcfFile}
fi