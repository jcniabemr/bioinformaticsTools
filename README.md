# Bioinformatics tools that i have written in python for various tasks  


1.) gffSort.py: For combining, sorting, duplicate removal and renaming features for 1+ gff files. Use bcftools intersect prior to combining multiple gff files 

	python gffSort.py --gff file1.gff file2.gff file3.gff fileX.gff


2.) createCDS.py: For extracting CDS from genome assembly and translating to protein. Requires gff and genome assembly and will produce proten and DNA CDS

	python createCDS.py --gff file.gff --fasta file.fasta --strainName output prefix name
    
3.) chunkGenome.py: For cutting up DNA genome or proteome into chunks. 

	python chunkGenome.py --infile genome.fasta --chunkSize 1000 

4.) trimDArTseqAdapters.py: For trimming the DarTseq barcode region from sequencing reads. Read files can be listed space seperated. 

	python trimDArTseqAdapters.py --reads file1.gz file2.gz file3.gz filex.gz 

5.) createSNPtiles.py: For counting the number of variants in a given tile size. 
	
	python createSNPtiles.py --vcf x.vcf --windowSize 1000

6.) renameContigs.py: For renaming contigs in a fasta file. Contigs will be renamed to any string given adter the --contigName flag

	python renameConigs.py --assembly x.fasta --contigName contig

7.) genomeToReads.py: For cutting an assmebled genome into reads of specified length and step. 

	python genomeToReads.py --g genome.fasta --readL 300 --readS 50

8.) findMimpSequences.py: For identifying MIMP sequences in Fusarium genomes. 

	python findMimpSequences.py --fasta assembly.fasta

9.) variantEffectPredictor.sh: For analysis of SNPs from VCF file to predict alterations to amino acid sequence, SNP type and severity of each SNP
	
	bash variantEffectPredictor.sh --vcf file.vcf --gff file.gff --fasta genomeAssembly.fasta 

10.) vcf_to_fasta.py: For producing a fasta sequence for each sample in a VCF file.
	
	python vcf_to_fasta.py --vcf file.vcf --ploidy 1 OR 2 

11.) vcf_to_matrix.py: For converting a VCF file into a matrix.

	python vcf_to_matrix.py --vcf <file.vcf> --ploidy <1 OR 2>

12.) vcf_similairty_matrix.py: For creation of a data matrix detailing an all vs all % sample similarity 
	
	python vcf_similairty_matrix.py --vcf file.vcf --ploidy 1 OR 2

13.) filterContifs.py: For filtering out small contigs from a genome assembly, renaming contigs and wrapping sequences
	
	python filterContigs.py --fasta <.fasta> --minLengh <500>  --contigName <contig> 

14.) calculateSequencingCoverage.py: For estimating sequencing coverge using fastq reads.

	python calculateSequencingCoverage.py -i file.fastq file2.fastq fileX.fastq -s <size Mb>

15.) mergeVCFs.py: For merging VCF files. 
	
	python mergeVCFs.py --vcf <.vcf1> <.vcf2> <.vcfx>  

16.) calculateContigLenghts.py: For creating a file of contig lengths for a genome assembly. 
	
	python calculateContigLengths.py --assembly <genomeAssembly.fa>

17.) filterVCFdepth.py: For filtering genotype depth in VCF file 

	python filterVCFdepth.py --vcf file.vcf --minDP <integer> --maxDP <integer>

18.) sortReadsByName.py: For sorting sequencing reads in fastq files by name as opposed to coordiante 

	python sortReadsByName.py --readF <_1.fastq> --readR <_2.fastq>

19.) cutGenome.py: Cut a genome into fragments between a defined region at defined intervals

	python cutGenome.py --assembly --start --end --step --window --contig
