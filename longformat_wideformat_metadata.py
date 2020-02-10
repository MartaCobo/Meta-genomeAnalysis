#!/usr/bin/python3.5


# This program join abundances of KOs in different samples and metadata from the samples.

# Usage: Change the paths of abundances table, metadata and outfile. 

# Author: Marta Cobo Simón


 
MuestraKEGGAbundancia={}
muestras = [];

# Read the table of KOs abundances in different samples in wide format and introduce the data into a dictionary.

with open('/home/marta/Desktop/Perl/Metagenomas_tesis/Pipeline/hipervariablesMADsdataframe.txt') as infile: # Rows: genes. Columns: samples.
    for line in infile:
        if 'MP' in line:
            muestras = line.strip().split('\t')
            muestrasordenadas = muestras.sort() 
        else:
            line = line.strip().split('\t')
            kegg = line.pop(0)
            for j in range(0, len(line)):
                muestra = muestras[j]
                if kegg not in MuestraKEGGAbundancia:
                    MuestraKEGGAbundancia[kegg] = { muestra: line[j] }
                else:
                    if muestra not in MuestraKEGGAbundancia[kegg]:
                        MuestraKEGGAbundancia[kegg][muestra] = line[j]
						
# Read the table of metadata in different samples in long format and introduce the information into a dictionary.

MuestraDepth = {}
MuestraMetadata = {}
with open ('/home/marta/Desktop/Perl/Metagenomas_tesis/Metagenoma_mock/Pelagibacter/metadatos.txt') as infile:
    for line in infile:
        if 'Code' not in line:
            line=line.strip().split('\t')
            muestra = 'MP' + line[0]
            depth = line[1]
            metadata = line[2]
            MuestraDepth[muestra] = depth
            MuestraMetadata[muestra] = metadata

# Join both informations and print the table in long format. 
			
for kegg in MuestraKEGGAbundancia:
    outputfile=open('/home/marta/Desktop/Perl/Metagenomas_tesis/Pipeline/geneshipervariablesMADs/' + kegg + '.txt','w')
    outputfile.write('{}\t{}\t{}\t{}\n'.format('muestra','cuenta','profundidad','metadato'))
    for muestra in MuestraKEGGAbundancia[kegg]:
        outputfile.write('{}\t{}\t{}\t{}\n'.format(muestra,MuestraKEGGAbundancia[kegg][muestra],MuestraDepth[muestra],MuestraMetadata[muestra]))
        #print (muestra + '\t' + kegg + '\t' + MuestraKEGGAbundancia[kegg][muestra] + '\t' + MuestraDepth[muestra] + '\t' + MuestraMetadata[muestra]) 
        

