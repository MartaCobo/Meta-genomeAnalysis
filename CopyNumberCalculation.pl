#!/usr/bin/perl
use strict;
use warnings;

my $PROGNAME = "CopyNumberCalculation.pl";
my $VERSION = "1.01";
my $USAGE=<<"USAGE" ;



#############################################################



 $PROGNAME $VERSION



This script calculates the copy number of different functions (KO, COG, PFAM) of complete genomes from PATRIC database. It requires the output of the program fun3.pl

 
 File format (example):
 
 protein	BESTHIT	BESTAVERAGE
 p1_1	K02010	K02010
 p1_2	K02011	K02011
 

 Usage

2- Indicate by arguments:

	First parameter: Input file (output of the program fun3.pl)
	Second parameter: Output file.
	
	perl $PROGNAME <input> <output>
	
 
 Example:
 
	perl $PROGNAME </home/Ecoli.fun3.txt> </home/Ecoli.copynumber.txt>
	
	
Author:

Marta Cobo Simón



#############################################################



USAGE



#############################################################

if (( $ARGV[0] eq "" ) or ( $ARGV[0] eq "-h") or ($ARGV[0] eq "--help")) { print $USAGE;  exit; }
else {


my $input = $ARGV[0];
my $output = $ARGV[1];


my %kegg_proteina;
open (my $infile, "$input" ) || die "I cannot open the file\n";
open (my $outfile, ">$output") || die "I cannot open the output\n";
print $outfile "#--Created by $0, ", scalar localtime, "\n";  # Write in the file the script that produced it and the date of execution. 
while(<$infile>) {
	unless ( $_ =~ /\#/) {
        chomp;
        my @columnas = split (/\t/, $_);
        my $proteina = $columnas[0];
		my $kegg = $columnas[2];
		if ( defined $kegg ) { $kegg_proteina{$kegg}{$proteina} = 1; } # Hash which related the protein and the KO. 
	}	
}
close $infile;

foreach my $kegg ( sort keys %kegg_proteina ) {
	my $copynumber = 0;
	foreach my $proteina ( sort keys %{$kegg_proteina{$kegg}} ) { # Calculation of copy number. 
		$copynumber++;
	}
	print $outfile "$kegg\t$copynumber\n";
}
}

