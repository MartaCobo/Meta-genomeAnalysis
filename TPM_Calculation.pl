
#!/usr/bin/perl
use strict;
use warnings;


my $PROGNAME = "TPM_Calculation.pl";
my $VERSION = "1.0";
my $USAGE=<<"USAGE" ;



#############################################################



 $PROGNAME $VERSION



 This script calculates the TPM (Transcripts Per Million of reads) of KOs from metagenomes per sample and per taxon. The reason for using TPM instead of RPKM is that the sum of all TPMs in each sample are the same. This makes it easier to compare the proportion of reads that mapped to a gene in each sample. In contrast, with RPKM and FPKM, the sum of the normalized reads in each sample may be different, and this makes it harder to compare samples directly.
 
 It requires: 
	- The following output files from SqueezeMeta (Tamames & Puente-Sánchez, 2019):
		- The file .orftable
		- The file .contiglog
	- A file with the selected samples to analyze, one sample in each row. Recommendation: use the samples where the selected taxa is abundant, for example, higher than 0.7% of the sample. Example:
		MP0311
		MP0313
		MP0321
		
	The result is a file called TPM_<selected_taxon>.txt, consisting on a table in long format saved in the output path selected, which contains the KO, the sample and the TPM. 
		
		KO	Sample	TPM
		K02010	MP0311	3000
		K02011	MP0311	4000
		K02010	MP0313	3500
		K02011	MP0313	4500

 Usage

	perl $PROGNAME <selected_taxon> <selected_samples_file> <orftable_file> <contiglog_file> <output_path> <samples_file> <keggcog>
	
	selected_taxon: Taxon selected to calculate the TPM. 
	selected_samples_file: File with the selected samples to do the calculation, each sample in a different row. 
	orftable_file: File .orftable from SqueezeMeta (Tamames & Puente-Sánchez, 2019).
	contiglog_file: File .contiglog from SqueezeMeta (Tamames & Puente-Sánchez, 2019).
	output_path: The path chosen to save the TPM table. 
	samples_file: File .samples from SqueezeMeta (Tamames & Puente-Sánchez, 2019).
	keggcog: Functional annotation based on kegg, cog or pfam. Options: kegg/cog/pfam.
 
 Example:
 
	perl $PROGNAME Pelagibacterales SamplesPelagibacterales.txt 12.coassembly.orftable 08.coassembly.contiglog /home/mcobo province10.samples kegg
	
	
Author:

Marta Cobo Simón, 2020

Doctoral thesis "Ecology of marine microorganisms: biodiversity, genomics and metagenomics".



#############################################################



USAGE



#############################################################

if (( $ARGV[0] eq "" ) or ( $ARGV[0] eq "-h") or ($ARGV[0] eq "--help")) { print $USAGE;  exit; }
 
my $generoescogido = $ARGV[0];
my $ficheroabundancias = $ARGV[1];
my $orftable = $ARGV[2];
my $contigstable = $ARGV[3];
my $outputpath = $ARGV[4];
my $ficheromuestras = $ARGV[5];
my $keggcogpfam = $ARGV[6];


# Selection of the samples.

my %muestras; 

open (my $infile, "$ficheroabundancias") || die "I cannot open the selected_samples file\n";
while(<$infile>) {
	chomp;
	$muestras{$_} = 1;	
}
close $infile;


# Selection from .contiglog of those contigs assigned to the selected taxon. They will be used to recover genes whose contig is assigned to the chosen taxon although the gene was not assigned. 

my %contigs;
open (my $infile2, "$contigstable") || die "I cannot open the .contiglog file\n";
while(<$infile2>) {
	unless  ( $_ =~ /Created/) {
		chomp;
		my ($contig, $taxonomiacompleta) = split(/\t/, $_);
		if ( $taxonomiacompleta =~ /$generoescogido/) {
			$contigs{$contig} = 1;
		}
	}
}  

close $infile2;

# From the .orftable, selection of the raw counts (number of reads mapped against one particular gene predicted from the assembly of the metagenome(s)). Select the genes assigned to the selected taxon or the genes whose contig is assigned to the selected taxon. 

my $numerofilas = 0;
open (my $infile3, "$ficheromuestras") || die "I cannot open the .samples file\n";
while(<$infile3>) {
	$numerofilas++;
}
my $numeromuestras = $numerofilas / 2; # Number of samples in the .orftable.
	
close $infile3;

my %SampleKEGGTaxaGeneRawcounts;
my %sample_kegg_genero_orf_length;
			
for ( my $i = -$numeromuestras; $i <= -1; $i++ ) {
	open (my $infile4, "$orftable") || die "I cannot open the .orftable file\n";
	while(<$infile4>) {
		if  ( $_ =~ /^ORF/) {
			chomp;
			my @lista = split(/\t/, $_);
			my $string = $lista[$i];
			my @strings = split (/\s/, $string );
			my $sample = $strings[-1];
			print "Analyzing $sample\n"; # Print the sample that is being analyzed.
			if ( defined $muestras{$sample} ) {
				while(<$infile4>) {
					unless (( $_ =~ /^ORF/) or ( $_ =~ /#/)) {
						chomp;
						my @lista2 = split(/\t/, $_);
						my $keggstring = '';
						if ( $keggcogpfam eq "kegg") { 
							if ( $lista2[6] =~ /K/) {
								$keggstring = $lista2[6]; 
							}
						}
						elsif ( $keggcogpfam eq "cog") {
							if ( $lista2[9] =~ /COG/) {
								$keggstring = $lista2[9];
							}
						}
						elsif ( $keggcogpfam eq "pfam") {
							if ( $lista2[12] =~ /PF/) {
								$keggstring = $lista2[12];
							}
						}
						my $lengthorf = $lista2[2];
						my @keggs = split(/;/, $keggstring);
						my $taxonomiacompleta = $lista2[5];
						my $orf = $lista2[0];
						my $contig = $lista2[1];
						my $rawcounts = $lista2[$i]; 
						if (( $keggstring =~ /\*/) && ( $rawcounts ne "") && ($rawcounts ne "0") && ($rawcounts ne "0.000") && ( @keggs) && ( $sample ne "" ) && ( $lengthorf ne "" ) && ( $taxonomiacompleta ne "" ) && ( $orf ne "")) { # best average
							my $numero = scalar (@keggs);
							my $newrawcounts = $rawcounts / scalar (@keggs); # Since one gene can be assigned to two or more KOs, the raw reads are divided by the number of KOs. 
							if (( $taxonomiacompleta =~ /$generoescogido/) or ( defined $contigs{$contig})) {
								my $genero = "$generoescogido";
								foreach my $kegg ( @keggs ) {
									$kegg =~ s/\*//; # Only KOs, COGs or PFAMs assigned with best average algorithm (SqueezeMeta), marked by '*', are recovered. The '*' is removed. 
									$SampleKEGGTaxaGeneRawcounts{$sample}{$genero}{$kegg}{$orf} = $newrawcounts;
									my $lengthorfgen = 3 * $lengthorf; # Calculation of gene length from protein length.
									$sample_kegg_genero_orf_length{$sample}{$kegg}{$genero}{$orf} = $lengthorfgen;
								}
							}
						}
					}
				}
			}
		}
	}
	close $infile4;
}

my %SampleGeneroKEGGRawcounts;

# Calculation of the abundance of each function (KO, COG or PFAM) in each sample for the selected taxon. 


foreach my $sample ( sort keys %SampleKEGGTaxaGeneRawcounts ) {
	foreach my $genero ( sort keys %{$SampleKEGGTaxaGeneRawcounts{$sample}} ) {
		foreach my $kegg ( sort keys %{$SampleKEGGTaxaGeneRawcounts{$sample}{$genero}} ) {
			my $keggreads = 0;
			foreach my $orf ( sort keys %{$SampleKEGGTaxaGeneRawcounts{$sample}{$genero}{$kegg}} ) {
				$keggreads += $SampleKEGGTaxaGeneRawcounts{$sample}{$genero}{$kegg}{$orf};
			}
			$SampleGeneroKEGGRawcounts{$sample}{$genero}{$kegg} = $keggreads;
		}
	}
}


# Calculation of the median of length of genes assigned to a function (KO, COG or PFAM) per sample. 

my %sample_kegg_genero_medianalongitudes;

foreach my $sample ( sort keys %sample_kegg_genero_orf_length ) {
	foreach my $kegg ( sort keys %{$sample_kegg_genero_orf_length{$sample}} ) {
        foreach my $genero ( sort keys %{$sample_kegg_genero_orf_length{$sample}{$kegg}} ) {
            my @longitudes = ();
            foreach my $orf ( sort keys %{$sample_kegg_genero_orf_length{$sample}{$kegg}{$genero}} ) {
                push @longitudes, $sample_kegg_genero_orf_length{$sample}{$kegg}{$genero}{$orf};
            }
            my $medianalongitudes = median (@longitudes);
            $sample_kegg_genero_medianalongitudes{$sample}{$kegg}{$genero} = $medianalongitudes;
		}
    }
}


# Calculation of TPM per sample and per taxon (to remove the bias of taxon abundance).

my %sample_genero_kegg_tpm;
my %sample_genero_kegg_abundancianormalizadalongitud;
my %sample_genero_total;

# Abundancies normalized by KO, COG or PFAM length.

foreach my $sample ( sort keys %SampleGeneroKEGGRawcounts ) {
	if ( defined $muestras{$sample} ) {
        foreach my $genero ( sort keys %{$SampleGeneroKEGGRawcounts{$sample}} ) {
            foreach my $kegg ( sort keys %{$SampleGeneroKEGGRawcounts{$sample}{$genero}} ) {
                my $abundancianormalizadalongitud = $SampleGeneroKEGGRawcounts{$sample}{$genero}{$kegg} * 1000 / $sample_kegg_genero_medianalongitudes{$sample}{$kegg}{$genero}; # longitud kegg por muestra
				$sample_genero_kegg_abundancianormalizadalongitud{$sample}{$genero}{$kegg} = $abundancianormalizadalongitud;
            }
        }
	}
}

# Sum of the abundances normalized by length.					

foreach my $sample ( sort keys %sample_genero_kegg_abundancianormalizadalongitud ) {
	foreach my $genero ( sort keys %{$sample_genero_kegg_abundancianormalizadalongitud{$sample}} ) {
		my $total = 0;
		foreach my $kegg (  sort keys %{$sample_genero_kegg_abundancianormalizadalongitud{$sample}{$genero}} ) {
			$total += $sample_genero_kegg_abundancianormalizadalongitud{$sample}{$genero}{$kegg};
		}
		$sample_genero_total{$sample}{$genero} = $total;
	}
}


# Calculation of TPM.

foreach my $sample ( sort keys %sample_genero_kegg_abundancianormalizadalongitud ) {
        foreach my $genero ( sort keys %{$sample_genero_kegg_abundancianormalizadalongitud{$sample}} ) {
		foreach my $kegg (  sort keys %{$sample_genero_kegg_abundancianormalizadalongitud{$sample}{$genero}} ) {
			my $tpm = $sample_genero_kegg_abundancianormalizadalongitud{$sample}{$genero}{$kegg} * 1000000 / $sample_genero_total{$sample}{$genero};
			$sample_genero_kegg_tpm{$genero}{$kegg}{$sample} = $tpm;
		}
	}
}

# Print of the long-format table.

foreach my $genero ( sort keys %sample_genero_kegg_tpm) {
	open( my $outfile, ">$outputpath/TPM_$genero.txt") || die "No puedo abrir el fichero output2\n";
	print $outfile "#--Created by $0, ", scalar localtime, "\n";  # Write in the file the script that produced it and the date of execution. 
	print $outfile "KO\tSample\tTPM\n";
	foreach my $kegg ( sort keys %{$sample_genero_kegg_tpm{$genero}} ) {
		foreach my $sample ( sort keys %{$sample_genero_kegg_tpm{$genero}{$kegg}} ) {
			print $outfile "$kegg\t$sample\t$sample_genero_kegg_tpm{$genero}{$kegg}{$sample}\n";
		}
	}
	close $outfile;
}

# Function to calculate the median.

sub median
	{
    		my @vals = sort {$a <=> $b} @_;
    		my $len = @vals;
    		if($len%2) #odd?
    		{
        		return $vals[int($len/2)];
   		}
    		else #even
    		{
        		return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    		}
	}
