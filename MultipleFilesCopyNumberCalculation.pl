#!/usr/bin/perl
use strict;
use warnings;

my $PROGNAME = "MultipleFilesCopyNumberCalculation.pl";
my $VERSION = "1.01";
my $USAGE=<<"USAGE" ;



#############################################################



 $PROGNAME $VERSION



 Part of the doctoral thesis "Ecology of marine microorganisms: biodiversity, genomics and metagenomics". 02/04/2020 Original version,

                            (c) Marta Cobo-Simón, CNB-CSIC.
 
 This script calculates the copy number of KOs in complete genomes from PATRIC database. It requires FASTA files of the complete genomes.
 

 Usage

	perl $PROGNAME <pathway>
	
	<pathway>: Pathway of FASTA files. 
 
 Example:
 
	perl $PROGNAME /home/genomes




#############################################################



USAGE



#############################################################

if (( $ARGV[0] eq "" ) or ( $ARGV[0] eq "-h") or ($ARGV[0] eq "--help")) { print $USAGE;  exit; }
else {



my $pathway = $ARGV[0];
my @names = split(/\//, $pathway);
my $foldername = $names[-1]; # Folder where FASTA files are placed. 

if (( $ARGV[0] eq "-h") or ( $ARGV[0] eq "--help") or ( $ARGV[0] eq "")) {
	print "Usage: perl buclecalculacopynumberrealkegg.pl <pathway>\n";
}
else {

	my $command="prodigal -a $pathway/$foldername.proteins.fasta -d $pathway/$foldername.cds.fasta -i $pathway/$foldername.fna -f gff -o $pathway/$foldername.gff"; # Gene prediction with Prodigal.
	system $command;
	
	my $command2="diamond blastp -q $pathway/$foldername.proteins.fasta -p 12 -d /home/fpuente/db/keggdb.dmnd -e 1e-03 --id 30 -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $pathway/$foldername.KEGG.m8 --comp-based-stats 0 --masking 0"; # Search of homologies with Diamond. 
	system $command2;

	my $command3="perl fun3assignadaptado.pl $pathway/$foldername.KEGG.m8 $pathway/$foldername.fun3kegg.txt"; # Functional assignment of proteins.
	system $command3;
	
	my $command4 = "perl CopyNumberCalculation.pl $pathway/$foldername.fun3kegg.txt $pathway/copynumber_bestaverage.txt"; # Calculation of copy number.
	system $command4;
}


