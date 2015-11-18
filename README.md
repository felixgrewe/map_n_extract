#map_n_extract

A script that maps similar reads to a reference sequence file and extracts genes by annotation (such as the [CEGMA](http://korflab.ucdavis.edu/datasets/cegma/) output) with the advantage that each individual gene sequence remains in-frame. This allows for a wide range of subsequent analyses; for example, splitting the data by codon position or selecting individual genes. Although we used and tested the program specifically for using the CEGMA output, any gff annotation and reference fasta sequence should work.

The script offers an option to build RAxML phylogenies from each individual alignment. However, individual alignments can also be merged into a supermatrix for further analysis.

If you use this code please cite:  
Steven D. Leavitt, Felix Grewe, Todd Widhelm, Lucia Muggia, Brian Wray, and H. Thorsten Lumbsch **Robust phylogenies of lichen-forming fungi inferred from genomes - overkill in the era of phylogenomics or reassurance?** *In prep*  

##Requirements:

Installed globally  
- Perl (with [BioPerl](http://www.bioperl.org/wiki/Main_Page))  
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)   
- [bcftools (with vcfutils.pl) and samtools](http://samtools.github.io/bcftools/)
- [BamTools](http://sourceforge.net/projects/bamtools/) (with vcfutils.pl)   
- [muscle](http://www.drive5.com/muscle/)  
- [raxmlHPC](http://sco.h-its.org/exelixis/web/software/raxml/index.html) (optional)  

additional Perl modules
- File::Which
- Getopt::Long
- Pod::Usage

-------------------------

##To run the script:

###Illumina sequence data
You can use either single-end or paired-end reads in this analysis, but all reads will be handled as single-ended. We recommend that you trim your reads before use; we used [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).  

###Have your files ready  
You need several files:  
1. A comma separated table that points to your sequence read files (see [seq_loc.csv](https://github.com/felixgrewe/CEGMA_CDS_extract/blob/master/seq_loc.csv) as an example)  
2. Reference fasta file; we used the CEGMA standard output file **output.cegma.dna**  
3. Annotation of the reference .fasta file in .gff format; we used the CEGMA standard output file **output.cegma.local.gff**  

###Run the script 
~~~
perl map_n_extract.pl  -t [-f -g -r -p -y -h]

 -table -t <.csv>			Comma seperated table with NGS sequence read location (name,se_or_pe1_reads[,se_or_pe2_reads])
 
 [-fasta -f <.fasta>]		reference sequence file (default=output.cegma.dna)
 [-gff -g <.gff>]			annotation of reference sequence file (default=output.cegma.local.gff)
 [-reference -r <STRING>]	Name of reference sequence (default=Ref)
 [-threads -p <INT>]		Number of available phreads (default=8)
 [-run_phylo -y]			If you want to run raxML phylogeny for individual genes
 [-help -h]					Displays this basic usage information 
~~~

>Please note that you can chose to run RAxML analyses for each individual gene alignment with option -y.

###Output
The script will produce two output folders. The **OUTPUT_alignments** folder holds all alignments of mapping across whole fasta reference sequences and the **OUTPUT_good_CDS_alignments** folder contains only gene alignments which are in frame with the open reading frame of the reference annotation.  
If you chose to analyze phylogenetic inferences for each individual gene alignment (option -y), you'll also see a **OUTPUT_phylogenies** folder that contains all RAxML trees.

###Further Processing

####Supermatrix
You can use the program [FASconCAT](https://www.zfmk.de/en/research/research-centres-and-groups/fasconcat) (Kueck & Meusemann, 2010) to build a supermatrix. Place the FASconCAT perl scrip into the *OUTPUT_good_CDS_alignment* folder, execute with
~~~
perl FASconCAT_v1.0.pl
~~~
and follow the instructions.

####Combining all trees
If you chose the -y option before and see the **OUTPUT_phylogenies** folder, you can combine all individual gene trees into one newick file (**all_tree.tre**) with:  
~~~
cat ./OUTPUT_phylogenies/*bipartitions\.*RAXML >> all_tree.tre
~~~

