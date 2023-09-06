Bioinformatics tools that i have written in python for various tasks  


1.) gffSort.py: For combining, sorting, duplicate removal and renaming features for 1+ gff files. Use bcftools intersect prior to combining multiple gff files 

	python gffSort.py --gff file1.gff file2.gff file3.gff fileX.gff


2.) createCDS.py: For extracting CDS from genome assembly and translating to protein. Requires gff and genome assembly and will produce proten and DNA CDS

	python createCDS.py --gff file.gff --fasta file.fasta --strainName output prefix name
    
3.) chunkGenome.py: For cutting up DNA genome or proteome into chunks. 

	python chunkGenome.py --infile genome.fasta --chunkSize 1000 

4.) trimDArTseqAdapters.py: For trimming the DarTseq barcode region from sequencing reads. Read files can be listed space seperated. 

	python trimDArTseqAdapters.py --reads file1.gz file2.gz file3.gz filex.gz 

4.) createSNPtiles.py: For counting tiles of variants given a tile size. 
	
	python createSNPtiles.py --vcf <vcf file> --windowSize <numberic window size>
