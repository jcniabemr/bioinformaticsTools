#!/usr/bin/env python 

####Script to rename genes in a gff where the ID field is not continous. This may occur when trying to merge multiple annotation files
####To run script please provide a list of files <infile.gff infile2.gtf infile3> containing annotations to be sorted, renamed and concatenated 
####Example run python order_genes.py -i gff1 gff2 gff3 -o outfile
####Written by John Connell  

####Import functions 
import argparse 

####Parse files 
ap = argparse.ArgumentParser()
ap.add_argument('-i',type=str,required=True,nargs='+',help='space seperated list of infiles to be sorted')
ap.add_argument('-o', type=str,help='outfile',default='genes_renamed_sorted')
args = ap.parse_args()
infiles = (args.i)
outname = (args.o)

####Set dicts 
gene_info = {}
features = {}

####Set lists 
all_data = []
contig_list = []
new_gff=[]
out_file_info = []

####Open and extract info from all files 
def open_extract (file):
	data = open(file)
	for y in data:
		y = y.replace("\n","")
		all_data.append(y)

for x in infiles:
	open_extract(x)

####Sort data into dicts and lists for later use
for y in all_data:
	if y.startswith("#"):
		continue 
	x = y.split("\t")
	if x[2] == "gene":
		contig = x[0] 
		start = int(x[3])
		contig_list.append(contig + "_" + str(start))
		if contig not in gene_info:
			gene_info[contig] = [start]
		else:
			gene_info[contig].append(start)
		contig_gene_info = ("_".join([contig, str(start)]))
	if contig_gene_info not in features:
		features[contig_gene_info] = [y]
	else:
		features[contig_gene_info].append(y)
			
###Sort data
for contig in sorted(contig_list, key = lambda x: (int(x.split("_")[1]), int(x.split("_")[2]))):
	new_gff.extend(features[contig])

####Rename genes
i=0
for x in new_gff:
	x=x.strip("\n")
	line=x.split()
	renamecol = line[8]
	if line[2] == "gene":
		gene_id = renamecol.split("=")[1]
		gene_id = gene_id.replace(";","")
		i+=1
		new_id = "g" + str(i)
	renamecol = renamecol.replace(gene_id, new_id)
	line[8] = renamecol
	out_file_info.append("\t".join(line))

outfile = (outname + ".gff")
with open(outfile, 'w') as f:
    for x in out_file_info:
        f.write(f"{x}\n")