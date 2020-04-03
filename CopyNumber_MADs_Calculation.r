# Part of the doctoral thesis "Ecology of marine microorganisms: biodiversity, genomics and metagenomics". 02/04/2020 Original version,

#                            (c) Marta Cobo-Simón, CNB-CSIC.

# This script calculates the average copy number per genome and highly-variable genes. 
# It requires a table of TPMs whose columns are samples and the rows are KOs or COGs. 


# Usage: source("CopyNumber_MADs_Calculation.R"). Choose the Universal Single Copy Genes (singleCopy) in KO or COG. 

# In the folder where the script is run, there must be the table with TPMs. The following outputs will be saved in the same folder:

# CopyNumber_longformat.txt -> Average copy number per genome of different functions at different samples in long format.

# CopyNumber_wideformat.txt -> Average copy number per genome of different functions at different samples in wide format.

# Highly-Variable.txt -> Highly-Variable genes defined as the ones whose Median Absolute Deviation (MAD) is higher than twice the 
# standard deviation.

# Highly-Variable.png -> Variation in copy number of Universal Single-Copy Genes (USiCGs) represented in boxplots for each sample in 
# comparison to highly-variable genes. 

# The variable kruskalwalliscopynumberusicgs contains the Kruskal-Wallis test for USiCGs. It must be non-significant. 

# The variable kruskalwalliscopynumbertodos contains the Kruskal-Wallis test for all the genes. 

# The variable foldchangeusicgsintersample contains the log2(foldchange) of USiCGs among different samples.

# The variable foldchangeusicgsintrasample contains the log2(foldchange) of USiCGs in each sample. 

# Test: test_CopyNumber_MADs_Calculation.txt as table of TPMs.


# Read the table of TPMs.
peleae10 = read.table("test_CopyNumber_MADs_Calculation.txt",header=T,sep="\t",row.names=1)


# Select which group of USiCGs must be used.

# Set of 15 OGs that are widely conserved (present in most species), rarely occur as duplicate genes in known genomes, are not 
# subject to horizontal gene transfer, and do not scale with genome size as marker genes  (Raes et al., 2007).
singleCopy = c('K02863','K02867','K02871','K02874','K02876','K02886','K02931','K02933','K02946','K02948','K02950','K02952','K02986',
	       'K02992','K02994') 
#singleCopy = c('COG0081','COG0080','COG0102','COG0093','COG0200','COG0090','COG0094','COG0097','COG0051','COG0100','COG0048','COG0099',
#               'COG0522','COG0049','COG0096')  

usicgspelae10reducidotabla <- subset(peleae10, rownames(peleae10) %in% singleCopy) # Extract only rows with USiCGs.
medianasusicgs <- colMeans(usicgspelae10reducidotabla) # Median of TPM USiCGs in each sample.  
copynumber <- t(t(peleae10) / medianasusicgs) # Calculation of Copy Number. 

write.table(copynumber, file="CopyNumber_wideformat.txt", sep="\t", quote=F) # Table in wide format with the copy number. 


# Calculation of MADs and selection of highly-variable genes. 

library(matrixStats)
MADs <- as.data.frame(rowMedians(as.matrix(abs(copynumber - rowMedians(as.matrix(copynumber))))), row.names=rownames(copynumber))
hipervariablesMADs <- subset(copynumber, subset=MADs > (2 * colSds(as.matrix(MADs)) + colMeans(as.matrix(MADs))))
MADstablatodos <- as.data.frame(as.matrix(abs(copynumber - rowMedians(as.matrix(copynumber))), row.names=rownames(copynumber)))


# Table format from wide format to long format. 

library(reshape2)
listacopynumber <- melt(copynumber)
colnames(listacopynumber) <- c("cog", "muestra", "rpkm") # Change the name of the columns.
listausicgspelae10reducido <- subset(listacopynumber, listacopynumber$cog %in% singleCopy) # Extract only rows with USiCGs. 
kruskalwalliscopynumberusicgs <- oneway.test(rpkm ~ muestra, data = listausicgspelae10reducido)
kruskalwalliscopynumbertodos <- oneway.test(rpkm ~ muestra, data = listacopynumber)
write.table(listacopynumber, file="listacopynumber.txt", sep="\t", quote=F)


# Highly-variable genes.

# Highly-variable genes and their TPMs in long format
hipervariables <- subset(listacopynumber, subset=listacopynumber$cog%in%rownames(hipervariablesMADs))

# Highly-variable genes and their TPMs in wide format.
hipervariablesimprimir <- subset(copynumber, subset=rownames(copynumber)%in%rownames(hipervariablesMADs)) 

write.table(hipervariablesimprimir, file="tablahipervariablesMADstpmcopynumber.txt", sep="\t", quote=F)


# Calculation of log(foldchange) of USiCGs.

# Extract USiCGs in wide format table. 
usicgspeleae10 <- subset(copynumber, subset=rownames(copynumber) %in% singleCopy)

# Calculation of log2(foldchange) of USiCGs among samples. 
foldchangeusicgsintersample <- apply(usicgspeleae10, 1, function(x) log2(max(x+0.00001)/min(x+0.00001))) 
				     
# Calculation of log2(foldchange) of USiCGs in each sample.				   
foldchangeusicgsintrasample <- apply(usicgspeleae10, 2, function(x) log2(max(x+0.00001)/min(x+0.00001))) 


# Plot with variation in copy number of USiCGs in each sample and copy number of highly-variable genes. 

png("Highly-Variable.png", width = 7*300,height = 7*300,res = 300,pointsize = 6)
colores <- c('#97B065', '#B93838', '#83AAF0', '#E2AA48', '#A9A3D2','#797900', '#F4BCE6', '#84F4F6', '#9AD63B', '#A285F5','#D2B48C', 
	'#4EA24E', '#465569','#9C669C', '#6495ED')                                				   
cogshipervariables <- hipervariables$cog
cogshipervariablesordenado <- sort(cogshipervariables)
cogshipervariablesordenadosinduplicados <- cogshipervariablesordenado[!duplicated(cogshipervariablesordenado)]
par(mar=c(9,9,9,9)+0.1,mgp=c(6,1,0))
plot(listausicgspelae10reducido$rpkm ~ listausicgspelae10reducido$muestra, ylim=c(0,9), lwd = 1, col="cornflowerblue", xlab= "Samples", 
	 ylab = "Average copy number per genome", tcl = 1, xaxt='n', cex.axis=2, cex.lab=2, las=1)
axis(1, at=1:length(unique(hipervariables$muestra)), labels = FALSE)
labels<-unique(hipervariables$muestra)
text(x=1:length(unique(hipervariables$muestra)), y=par("usr")[3]-0.2,labels=labels, srt=45, adj=1, xpd=TRUE, cex=1.5)
for (i in 1:length(unique(cogshipervariables))) {
	par(new=T)
	cog <- cogshipervariablesordenadosinduplicados[i]
	cognuevo <- subset(hipervariables, subset=hipervariables$cog == cogshipervariablesordenadosinduplicados[i]) 
	logrpkm = cognuevo$rpkm
   	names(logrpkm) = cognuevo$muestra
		plot(1:length(logrpkm), logrpkm, lwd = 2.5, ylim=c(0,9), type='l', col=colores[i], ylab = "", axes=F, xlab= "", tcl = 1, 
		xaxt='n', yaxt='n')
    
}
legend("topright", legend=cogshipervariablesordenadosinduplicados, col=colores, lwd=2.5, ncol=1, cex=2)

dev.off()
				     
