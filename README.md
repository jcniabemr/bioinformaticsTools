Bioinformatics tools that i have written in python for various tasks  


1.) gffSort.py: For combining, sorting, duplicat removal and renaming features for 1+ gff files.
	gffSort.py --gff file1.gff file2.gff file3.gff fileX.gff


2.) createCDS.py: For extracting CDS from genome assembly and translating to protein. Requires gff and genome assembly and will produce proten and DNA CDS
	createCDS.py --gff file.gff --fasta file.fasta --strainName output prefix name
    
