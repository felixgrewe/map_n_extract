#!/usr/bin/perl -w
#
# map_n_extract.pl
#
# Discription  : Maps NGS reads to gene regions and extracts CDS for phylogenetic analyses.
#
# OS prog req  : bowtie2, samtools, bam2consensus, muscle, raxmlHPC
#
# files req    : csv table, output.cegma.dna, output.cegma.local.gff
# 
# Author       : Felix Grewe (fgrewe@fieldmuseum.org), www.felixgrewe.de
# Date created : 05/29/2015
#

use strict;

use File::Which;
use Getopt::Long;
use Pod::Usage;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Align::Utilities qw(:all);
use Bio::Seq;

##################
# Get input data #
##################

### setup my defaults
our ($fasta, $gff, $table, $reference, $cores, $run_phylo, $help);
$fasta = "output.cegma.dna";
$gff = "output.cegma.local.gff";
$reference = "Ref";
$run_phylo = undef;
$cores = 8;

GetOptions(	'fasta|f=s' 	=> \$fasta,
			'gff|g=s'		=> \$gff,
			'table|t=s'		=> \$table,
			'reference|r=s'	=> \$reference,
			'threads|p=i'	=> \$cores,
			'run_phylo|y'	=> \$run_phylo,
			'help|h|?'		=> \$help,
);

pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1) unless $table;

######################
# Check dependencies #
######################

which("bowtie2") or die "Seems like bowtie2 is not installed.\n";
which("samtools") or die "Seems like samtools is not installed.\n";
which("bcftools") or die "Seems like bcftools is not installed.\n";
which("vcfutils.pl") or die "Seems like vcfutils.pl is not installed.\n";
which("muscle") or die "Seems like muscle is not installed.\n";
if (defined($run_phylo)) {
	which("raxmlHPC") or die "Seems like raxmlHPC is not installed.\n";
}
#############################
# Load sequence information #
#############################

### Reads reference cegma fasta

my @ref_genes = ();
my $ref = 0;
my $id = 0;
my %seq = ();

open (INPUT_CEGMA, "$fasta") or die "Can't open cegma_fasta with reference reads!!!\n";

while (my $line = readline(INPUT_CEGMA)) { # reads the fasta
    chomp ($line);
    if ($line =~ m/^>/) {      
        #$fasta_name[$seq_cnt++] = $line; # reads fasta name and creats array - sugg by JJR      
        $id = $line; # reads fasta and creates hash - sugg by JJR
        $seq{$id} = ""; # empties values for key   
    } else {     
        #$fasta_seq[$seq_cnt] .= $line; # reads fasta sequence and fills array - sugg by JJ   
        $seq{$id} .= $line; # reads fasta sequence and fills hash - sugg by JJR
    }
}

close INPUT_CEGMA or die "Can't close cegma_fasta!!!\n";

### Sequence location

my @taxon = ();
my @pe_path = ();
my @pe_path2 = ();

open (INPUT_TABLE, "$table") or die "Can't open table with sequencing information!!!\n";

while (my $line = readline(INPUT_TABLE)) {
    chomp $line;
    
    if ($line =~ m/^#/) { # doesn't read comments in the table               
        next;       
    } else {        
        #writing into arrays
        my @table = split ",", $line; 	# splitting the blast output and reading in arrays
		push (@taxon, $table[0]);			# name of taxon
        push (@pe_path, $table[1]);    		# absolute path to sequence location
		push (@pe_path2, $table[2]);    	# absolute path to other sequence location
    }
}  

close INPUT_TABLE or die "Can't close table!!!\n";


### CEGMA annotation

my @cegma_scaffold = ();
my @cegma_exon_pos = ();
my @cegma_start = ();
my @cegma_stop = ();
my @cegma_gene = ();

my @cegma_aln_start = ();
my @cegma_aln_stop = ();

my %fas_loc = ();

open (INPUT_GFF, "$gff") or die "Can't open cegma gff file!!!\n";
while (my $line = readline(INPUT_GFF)) {
    chomp $line;
    
    if ($line =~ m/^#/) { # doesn't read comments in the table               
        next;       
    } else {        
        #writing into arrays
        my @table = split (/\t/, $line); 
#        push (@cegma_scaffold, $table[0]);
		push (@cegma_exon_pos, $table[2]);
		push (@cegma_start, $table[3]);
		push (@cegma_stop, $table[4]);
		push (@cegma_gene, $table[0]);
    }
}  

close INPUT_GFF or die "Can't close gff!!!\n";

### output Folder
my $ali_folder = "OUTPUT_alignments";
mkdir $ali_folder;
my $CDS_folder = "OUTPUT_good_CDS_alignments";
mkdir $CDS_folder;

##################
# Bowtie mapping #
##################

for (my $i = 0; $i < @taxon; $i++) {
	my $fastq_reads;
	
	mkdir $taxon[$i];
	chdir ($taxon[$i]) or die "Can't change into $taxon[$i]";
		
	if (exists($pe_path2[$i])){
		$fastq_reads = $pe_path[$i].",".$pe_path2[$i];
		print "Found multiple sequence files for: $taxon[$i]\n";
	} else {
		$fastq_reads = $pe_path[$i];
		print "Found single sequence files for: $taxon[$i]\n";
	}

	print "###Start Bowtie2 for $taxon[$i]\n\n";
	
	my $bowtieBuild = 'bowtie2-build ../'.$fasta.' '.$taxon[$i];
	my $bowtie2 = 'bowtie2 -x '.$taxon[$i].' -U '.$fastq_reads.' -N 1 -S '.$taxon[$i].'.sam -p '.$cores;
	my $samtoolsView = 'samtools view -bS '.$taxon[$i].'.sam > '.$taxon[$i].'.bam';
	my $samtoolsSort = 'samtools sort '.$taxon[$i].'.bam '.$taxon[$i].'.sorted';
	my $samtoolsIndex = 'samtools index '.$taxon[$i].'.sorted.bam';
	my $samtoolsMpileupUF = 'samtools mpileup -uf ../'.$fasta.' '.$taxon[$i].'.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > '.$taxon[$i].'_vcfcons.fastq';
	my $copy_vcf = 'cp '.$taxon[$i].'_vcfcons.fastq ../'.$ali_folder;
	
	system ($bowtieBuild);
	system ($bowtie2);
	system ($samtoolsView);
	system ($samtoolsSort);
	system ($samtoolsIndex);
	system ($samtoolsMpileupUF);
	system ($copy_vcf);
	
	print "###Bowtie2 for $taxon[$i] finished\n\n";
	
	chdir ("..");
}

###################
# Build alignment #
###################

chdir ($ali_folder) or die "Can't change into $ali_folder";

print "Building alignments\n\n";

foreach my $key (sort(keys %seq)) {
	my ($gene_name) = $key =~ m/^\>(.*)/;
#	print "Gene name:\n";
#	print $gene_name."\n";
	open (OUTPUT_FASTA, ">".$gene_name."_region.fas") or die "Can't open output fasta file!!!\n";
	print OUTPUT_FASTA ">$reference\n";
	print OUTPUT_FASTA $seq{$key}."\n";
	close OUTPUT_FASTA or die "Can't close output fasta from reference round!!!\n";
}

for (my $i = 0; $i < @taxon; $i++) {

	my $mapping_id = 0;
	my %mapping_seq = ();

	my $inseq = Bio::SeqIO->new(-file => "<".$taxon[$i]."_vcfcons.fastq", -format => "fastq");
	
	while (my $seq = $inseq->next_seq) {
		$mapping_id = $seq->id;
		$mapping_seq{$mapping_id} = $seq->seq();
	}
	
	foreach my $mapping_key (sort(keys %mapping_seq)) {
		my $gene_name = $mapping_key;
		
		open (OUTPUT_FASTA, ">>".$gene_name."_region.fas") or die "Can't open output fasta file!!!\n";
		print OUTPUT_FASTA ">$taxon[$i]\n";
		print OUTPUT_FASTA $mapping_seq{$mapping_key}."\n";
		close OUTPUT_FASTA or die "Can't close output fasta from mapping round!!!\n";
	}
}

print "All alignmtes are done\n\n";

########################################
# run muscle on all CEGAM region files #
########################################

foreach my $key (sort(keys %seq)) {
	my ($gene_name) = $key =~ m/^\>(.*)/;
	my $file_name = $gene_name."_region.fas";
	my $file_out = $gene_name."_regions_muscle.fas";
	my $muscle = 'muscle -in '.$file_name.' -out '.$file_out;
	
	system ($muscle);
}
	
########################
# Alignment processing #
########################

# chdir ($ali_folder) or die "Can't change into $ali_folder";

my %sliced_alis = ();
my %combined_alis = ();
my %sliced_no = ();
my $current_gene = $cegma_gene[0];

for (my $i = 0; $i < @cegma_gene; $i++) {
	my $to_open = $cegma_gene[$i]."_regions_muscle.fas";
	my $str = Bio::AlignIO->new(-file => $to_open, -format => "fasta");
	my $aln = $str->next_aln();
	my $pos_start = $aln->column_from_residue_number($reference, $cegma_start[$i]);
	push (@cegma_aln_start, $pos_start); 			#saving new start position in array
	my $pos_stop = $aln->column_from_residue_number($reference, $cegma_stop[$i]);
	push (@cegma_aln_stop, $pos_stop); 				#saving new stop positions in array
	my $sliced = $aln->slice($pos_start, $pos_stop, 1);
	$sliced_no{$cegma_gene[$i]}++;
	$sliced_alis{$cegma_gene[$i]}{$sliced_no{$cegma_gene[$i]}} = $sliced;
	$combined_alis{$cegma_gene[$i]} = "";
}

my $ali123 = 0;
my $Ref_seq = 0;
my $Ref_prot = 0;

foreach my $gene (sort keys %sliced_alis) {

	my $count = 0;
	foreach my $number (sort keys %{ $sliced_alis{$gene} }) {
		if ($number == 1) {
			$ali123 = $sliced_alis{$gene}{$number};
		} elsif ($number > 1) {
			$ali123 = cat($ali123, $sliced_alis{$gene}{$number});
		}
	}
	
	### Getting rid of all other insertions in all other sequences
	my $ali123_ref = $ali123->set_new_reference($reference);
	my $ali123_nogap = $ali123_ref->remove_gaps;
	$ali123_nogap->set_displayname_flat();
	
	### Check if alignment is complete by ORF
	foreach my $seq_obj ($ali123_nogap->each_seq_with_id($reference)) {
		$Ref_seq = $seq_obj->seq();
		$Ref_prot = ($seq_obj->translate)->seq();
	}
	$count = ($Ref_prot =~ tr/\*//);
	
	if ($Ref_seq !~ /^ATG/) {
		print "$gene does not start with ATG and was excluded from further analysis\n";
	} elsif ($Ref_seq !~ /TGA$|TAA$|TAG$/) {
		print "$gene does not stop with TGA, TAA, or TAG and was excluded from further analysis\n"; 
	} elsif ($count > 1) {
		print "$gene has $count stop codons in translation and was excluded from further analysis\n";
	} elsif ((length($Ref_seq)%3) != 0) {
		print "$gene is not divisible by 3 and was excluded from further analysis\n";
	} else {
		my $to_close_fas = "../".$CDS_folder."/".$gene."_CDS.fas";
		my $out_fas = Bio::AlignIO->new(-file => ">".$to_close_fas, -format => "fasta");
		$out_fas->write_aln($ali123_nogap); 

		my $to_close_phy = "../".$CDS_folder."/".$gene."_CDS.phy";
		my $out_phy = Bio::AlignIO->new(-file => ">".$to_close_phy, -format => "phylip");		
		$out_phy->write_aln($ali123_nogap); 
	}
}
chdir ("..");

############################################
# Build RaxML phylogenies on all CEGMA CDS #
############################################

if (defined($run_phylo)) {
	my $phy_folder = "OUTPUT_phylogenies";
	mkdir $phy_folder;
	chdir ($phy_folder) or die "Can't change into $phy_folder";
	
	my $move = 'mv ../'.$CDS_folder.'/*_CDS.phy .';
	system ($move);

	foreach my $gene_name (sort(keys %sliced_alis)) {
		
		my $file_out = $gene_name."_CDS.phy";
		my $rax_out = $gene_name."_CDS_RAXML";

		my $raxml = 'raxmlHPC -s '.$file_out.' -n '.$rax_out.' -m GTRGAMMA -f a -p 194955 -x 12345 -# 100 -T '.$cores;
		
		system ($raxml);
	}
	chdir ("..");
} else {
	my $del = 'rm ./'.$CDS_folder.'/*_CDS.phy';
	system ($del);
}

print "All done. Enjoy!\n";


__END__

=head1 NAME

CEGMA_CDS_extract

=head1 COPYRIGHT

   copyright (C) 2015 Felix Grewe

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Maps NGS reads to annotated gene regions and extracts CDS for phylogenetic analyses.

=head1 SYNOPSIS

perl map_n_extract.pl  -t [-f -g -r -p -y -h]

 -table -t <.csv>		Comma seperated table with NGS sequence read location (name,se_or_pe1_reads[,se_or_pe2_reads])
 
 [-fasta -f <.fasta>]		reference sequence file (default=output.cegma.dna)
 [-gff -g <.gff>]		annotation of reference sequence file (default=output.cegma.local.gff)
 [-reference -r <STRING>]	Name of reference sequence (default=Ref)
 [-threads -p <INT>]		Number of available phreads (default=8)
 [-run_phylo -y]		If you want to run RAxML phylogeny for individual genes
 [-help -h]			Displays this basic usage information 
 


