#!/usr/bin/env python 

########################################################################
# Script to produce similarity matrix of similair sites for all        #
# samples in a VCF file                                                #
#                                                                      #
# Usage example:                                                       #
# python vcf_to_pca.py --vcf <.vcf>                                    #
#                                                                      #
# Written by John Connell                                              #
# john.connell@niab.com                                                #
# NIAB                                                                 #
# 2023                                                                 #
########################################################################


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument(
    '--vcf',
    type = str,
    required = True,
    help = 'VCF file'
    )
  ap.add_argument(
    '--ploidy',
    type = str,
    required = True,
    help = 'organism ploidy'
    )
  parse = ap.parse_args()
  if parse.ploidy == "1":
    createPCA(parse.vcf)
  elif parse.ploidy == "2":
    createPCA(parse.vcf)
  else:
    print("Ploidy argumnent :",parse.ploidy," is not aloud, please use either 1 or 2")
    exit(1)


if __name__ == "__main__":
  import os,argparse
  main()




