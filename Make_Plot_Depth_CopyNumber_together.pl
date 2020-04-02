#!/usr/bin/perl
use strict;
use warnings;

my $PROGNAME = "Make_plot_Depth_CopyNumber_together.pl";
my $VERSION = "1.01";
my $USAGE=<<"USAGE" ;



#############################################################



 $PROGNAME $VERSION



 Part of the doctoral thesis "Ecology of marine microorganisms: biodiversity, genomics and metagenomics". 02/04/2020 Original version,

                            (c) Marta Cobo-Simón, CNB-CSIC.
 
 This script creates scatterplots to relate the average copy number per genome of different functions from metagenomes to the sample depth. All the plots are located to the same file. Files must contain the function, the sample and the copy number. Files names must contain the functional class (nutrient in this case) and the station. 

 File name: ironlowaffinityestacion19
 
 File format (example):
 
 kegg	sample	copynumber
 K02010	MP0311	3.1
 K02011	MP0311	2.0
 K02010	MP0323	1.1
 K02011	MP0323	2.1
 

 Usage

1- Modify the path of the files (path_files) and the plots (path_plots).

2- Indicate by arguments:

	First parameter: functional class (irondepthvariation/ironhemeenterobactin/irondicitratetonB/ironlowaffinity)
	Second parameter: Maximum of x axis
	Third parameter: Interval used in x axis
	Fourth parameter: Number of rows in the file to organize the plots.
	Fifth parameter: Number of columns in the file to organize the plots. 
	Sixth parameter: legend position in the file (outside the plot) (top/right/left/bottom/center/topleft/topright/bottomleft/bottomright)

	perl $PROGNAME <functional_class> <xlim> <interval> <rows> <columns> <legend_position>
 
 Example:
 
	perl $PROGNAME irondepthvariation 6 2 4 3 bottom




#############################################################



USAGE



#############################################################

if (( $ARGV[0] eq "" ) or ( $ARGV[0] eq "-h") or ($ARGV[0] eq "--help")) { print $USAGE;  exit; }
else {


my $nutriente = $ARGV[0];
my $xlim = $ARGV[1];
my $interval =$ARGV[2];
my $rows = $ARGV[3];
my $columns = $ARGV[4];
my $legend_position = $ARGV[5];

my $height = 4 * $rows;  # Heigth of the file containing the plots.
my $width = 4 * $columns; # Width of the file containing the plots. 



my $path_files = "/media/Backup/disk1/marta/coassemblylonghurstjuntos/Pelagibacterales/Masmuestras0.7nuevosilva/kegg/transportadores/estacionlista4";  # Path of the infile files.
my $path_plots = "/media/Backup/disk1/marta/coassemblylonghurstjuntos/Pelagibacterales/Masmuestras0.7nuevosilva/kegg/transportadores/plotestacionlista4"; # Path of the output files.

opendir(my $indir,$path_files); 
my @dirs=readdir $indir;  # Read the files in the directory indicated by the path and introduce them in the array @dirs. 
closedir $indir;

# This script uses R to make the plots. To do that, it creates a file where all the commands for R will be written and, at the end, executes the commands in R. 

open (my $outfile, ">$path_plots/matrix") || die "I cannot open the outfile\n"; # Create a file where the commands for R will be written.  
print $outfile "#--Created by $0, ", scalar localtime, "\n";  # Write in the file the script that produced it and the date of execution. 
print $outfile "svg(\"$path_plots/$nutriente.svg\", width=$width, height=$height)\n"; # File with the plots, the width and the height are established. 
print $outfile "par(mfrow=c($rows,$columns))\n"; # Number of rows and columns of the plots file. 
print $outfile "par(oma = c(5, 5, 1, 1))\n"; # make room (i.e. the 4's) for the overall x and y axis titles
print $outfile "par(mar = c(2, 2, 1, 1))\n"; # make the plots be closer together
my $numero = length($nutriente);
my $numerototal = $numero + 8;
foreach my $mydir(sort { substr(my $a, $numerototal) <=> substr(my $b, $numerototal)  } @dirs) { # Sort the files by the station number
	if ( $mydir =~ /$nutriente/ ) { 
		my @nombre = split(/$nutriente/, $mydir);
		my $nombrebonito = $nombre[1];
		$nombrebonito =~ s/estacion/Station /g; # Set the name of the plot. 
		$newpath="$path_files/$mydir";  
		open (my $infile, "$newpath") || die "I cannot open the file\n"; # First plot, corresponding to one station. 
		print $outfile "lista<-read.table(\"$newpath\",header=TRUE, sep=\"\t\")\n";
		print $outfile "colores <- c(\'red\',\'green\',\'blue\',\'purple\')\n";
		print $outfile "keggnuevo <- subset(lista, subset=lista\$kegg==lista\$kegg[1])\n"; #Select the first function of the file.
		print $outfile "copynumberkegg <- keggnuevo\$copynumber\n";
		print $outfile "names(copynumberkegg) = keggnuevo\$sample\n";
		print $outfile "par(mai=c(0.3,0.6,0.5,0.1))\n";
		print $outfile "plot(keggnuevo\$copynumber, keggnuevo\$sample, lwd = 3, col=colores[1], type=\"b\", ylim = rev(c(0, 4000)), xlim = c(0,$xlim), main = \"$nombrebonito\", xaxt=\'n\', las=1, cex.axis=1.5, ylab = \"\", xlab = \"\", cex.main=2)\n"; # First line in the plot, corresponding to the first function, without x axis.
		print $outfile "labels<-seq(from = 0, to = $xlim, by = $interval)\n";
		print $outfile "axis(1, at=labels,labels=labels, cex.axis=1.5)\n"; # Set the x axis.
		print $outfile "for (i in 1:length(unique(lista\$kegg))) {\n";
		print $outfile "par(new=T)\n";
		print $outfile "keggnuevo <- subset(lista, subset=lista\$kegg==lista\$kegg[i])\n";
		print $outfile "copynumberkegg <- keggnuevo\$copynumber\n";
		print $outfile "names(copynumberkegg) = keggnuevo\$sample\n";
		print $outfile "plot(keggnuevo\$copynumber, keggnuevo\$sample, lwd = 3, col=colores[i], type=\"b\", ylim = rev(c(0, 4000)), xlim = c(0,$xlim), xaxt='n', yaxt='n', xlab=\"\", ylab=\"\")\n"; # Add new lines to the same plot, corresponding to the rest of the functions. 
		print $outfile "}\n";
	}
}

print $outfile "mtext(\"Average copy number per genome\", side = 1, outer = TRUE, line = 1.5, cex=1.5)\n"; # Text for the global x axis (whole file).
print $outfile "mtext(\"Depth (m)\", side = 2, outer = TRUE, line = 1.5, cex=1.5)\n"; # Test for the global y axis (whole file).
print $outfile "par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)\n"; # Set the position of the plots in the file. 
print $outfile "plot(0, 0, type = \"n\", bty = \"n\", xaxt = \"n\", yaxt = \"n\")\n";
print $outfile "legend(\"$legend_position\", legend = unique(lista\$kegg), xpd = TRUE, inset = c(9, 0.08), bty = \"n\", col = colores, cex = 2, title=\"Genes\", lwd=3, ncol=1)\n"; # Legend
print $outfile "dev.off()\n";
close $outfile;
system("R --slave --no-save --no-restore --no-environ --silent < $path_plots/matrix > $path_plots/plot"); # Send the file with the commands to R. 

}

