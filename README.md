#CEGMA_CDS_extract

A script that maps similar reads to a reference sequence file and extracts genes by annotation (such as the [CEGMA](http://korflab.ucdavis.edu/datasets/cegma/) output) with the advantage that each individual gene remains in frame with the genetic code. This allows for subsequent analyses for example to split the data by codon position or to select for individual genes. Although we used and tested the program specificlly for the CEGMA output each gff annotated and reference sequence should work.

The script offers an option to build a RAxML phylogeny of each individual alignment. However, individual alignments can also be merged into a supermatrix for further analysis.

If you use this code please cite:  
Steven D. Leavitt, Felix Grewe, Todd Widhelm, Lucia Muggia, Brian Wray, and H. Thorsten Lumbsch **Robust phylogenies of lichen-forming fungi inferred from genomes - overkill in the era of phylogenomics or reassurance?** *In prep*  

##Requirements:

Installed globally  
- Perl (with BioPerl)  
- bowtie2  
- samtools   
- bam2consensus (from the bambam package)  
- muscle  
- raxmlHPC (optional)  

-------------------------

##To run the script:

###Illumina sequence data
You can use either single-end or paired-end reads in this analysis but all reads will be handled as single-ended. We recommend to trim your reads before use; e.g. with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).  

###Have you files ready  
You need several files to get going:  
1. A comma separated table that points to your sequence read files (see [seq_loc.csv](https://github.com/felixgrewe/CEGMA_CDS_extract/blob/master/seq_loc.csv) as an example)  
2. Reference fasta file; we used the CEGMA standard output file *output.cegma.dna*  
3. Annotation of the reference fasta file (see 2) in gff format; we used the CEGMA standard output file *output.cegma.local.gff*  

###Run the script 
~~~
CEGMA_CDS_extract.pl  -t [-f -g -r -p -y -h]

 -table -t <.csv>			Comma seperated table with NGS sequence read location (name,se_or_pe1_reads[,se_or_pe2_reads])
 
 [-fasta -f <.fasta>]		reference sequence file (default=output.cegma.dna)
 [-gff -g <.gff>]			annotation of reference sequence file (default=output.cegma.local.gff)
 [-reference -r <STRING>]	Name of reference sequence (default=Ref)
 [-threads -p <INT>]		Number of available phreads (default=8)
 [-run_phylo -y]			If you want to run raxML phylogeny for individual genes
 [-help -h]					Displays this basic usage information 
~~~

>Please note that you chose if you want to run RAxML analyses for each individual gene alignment.

###Output
The script will produce two output folders. The *OUTPUT_alignments* folder holds all alignments of mapping across whole fasta reference sequences and the *OUTPUT_good_CDS_alignments* folder contains all good gene alignments which are in frame to the open reading frame of the reference annotation. 
If you chose to analyze phylogenetic inferences for each individual good gene alignment (option -y), you'll also see a *OUTPUT_phylogenies* folder that contains all RAxML trees.

###Further Processing

####Supermatrix
You can use the program [FASconCAT](https://www.zfmk.de/en/research/research-centres-and-groups/fasconcat) (Kueck & Meusemann, 2010) to build a supermatrix. Place the FASconCAT perl scrip into the *OUTPUT_good_CDS_alignment* folder and execute with
~~~
perl FASconCAT_v1.0.pl
~~~
and follow the instructions.

####Combining all trees
If you chose the -f option before and see the *OUTPUT_phylogenies* folder, you can combine all individual gene trees into one newick file with:  
~~~
cat ./OUTPUT_phylogenies/*\.bipartitions\.*RAXML .
~~~

