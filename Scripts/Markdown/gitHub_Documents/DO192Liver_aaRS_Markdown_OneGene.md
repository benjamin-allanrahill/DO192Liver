DO192Liver\_aaRS\_Markdown
================
Ben Allan-Rahill
6/25/2018

Running QTL scans, with covariates, for one aaRS gene
=====================================================

### Install & load packages

Here we will install the packages needed for running the QTL mapping and plotting

``` r
# Install packages
# This is all from the Manhattan/AE plot script.
install.packages("devtools")
install.github("dmgatti/DOQTL")
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Mmusculus.UCSC.mm10")
install.github("simecek/intermediate")
install.github("kbroman/qtlcharts")
install.packages(tidyverse)
```

Load packages:

``` r
library(devtools)
library(qtl2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
```

### Load the data files

``` r
load("/Users/c-allanb/Desktop/Benjamin_Allan-Rahill/DO192Liver/Data/Input/DO192LiverData_Formattedfor_rQTL2.Rdata")
```

Scan
====

If running the script on the HPC, you will need to specify the number of cores needed for the job.

``` r
## Establish global variables that do not need to be the for loop
  
# Set the number of processor cores to be used. Use 1 on the HPC.
cores.num = 1
```

### Select gene name, index, and chromosome

For our purposes, we will select the first value (Aars) from the list of aaRS genes. We will then retrieve the index of the genes from the annotations dataframe. This will allow us to subset based off the index later on.

``` r
gene_name <- "Aars"

# Get the gene's index in the protein and rna data.
index.protein <- which(annotations.protein.192$Associated.Gene.Name == gene_name)
index.rna <- which(annotations.rna.192$Gene == gene_name)

# Get the chromosome that the gene is on
gene_chr <- annotations.protein.192[index.protein,]$Chromosome.Name  
```

### Kinship matricies

Here we will establish the kinship matricies that will be used by the scan funtions later on. We will need to subset based off the mice that have data for protein and RNA we are looking at.

``` r
# Only subset mice with protein expression data. I don't remember what issue this prevented.
mice.with.data.r <- !is.na(expr.rna.192[,index.rna])
mice.with.data.p <- !is.na(expr.protein.192[,index.protein])
    
# Set up kinship matrices for the mice that have data for your given RNA/protein. Does it need to be for all chromosomes?
K.LOCO.192.chr.r <- K.LOCO.qtl2
K.LOCO.192.chr.p <- K.LOCO.qtl2

for (chromosome in 1:20){
  K.LOCO.192.chr.r[[chromosome]]<- K.LOCO.qtl2[[chromosome]][mice.with.data.r,mice.with.data.r]
}
for (chromosome in 1:20){
  K.LOCO.192.chr.p[[chromosome]] <- K.LOCO.qtl2[[chromosome]][mice.with.data.p,mice.with.data.p]
}  
```

### Establish additive and interactive covaraints

This is important for the scan function to run. It will factor in the effects of the covariants (sex, diet, and tag) on the QTL scans

``` r
# Set up additive and interactive covariates for RNA/protein, only mice with data.
addcovar.r <- model.matrix(~ Sex + Diet, data=covariates.rna.192[mice.with.data.r, ])
intcovar.sex.r <- model.matrix(~ Sex, data=covariates.rna.192[mice.with.data.r, ])
intcovar.diet.r <- model.matrix(~ Diet, data=covariates.rna.192[mice.with.data.r, ])

addcovar.p <- model.matrix(~ Sex + Diet + Tag, data=covariates.protein.192[mice.with.data.p, ])
intcovar.sex.p <- model.matrix(~ Sex, data=covariates.protein.192[mice.with.data.p, ])
intcovar.diet.p <- model.matrix(~ Diet, data=covariates.protein.192[mice.with.data.p, ])
```

### Running the scan

Here the qtl2 function 'scan1' will be used to creae LOD scores for the psuedomarkers. This will allow us to create plots with this data later on. The 'scan1' function also factors in the effects of kinship and covariants on the data.

``` r
# Perform QTL/LOD scans for RNA/protein in conditions of additive, interactive diet, and interactive sex covariates.
    
# RNA
# Scan / RNA / Additive covariate: Sex & Diet
lod.rna.addcovariatesonly <- scan1(genoprobs=probs[mice.with.data.r,],
                                         kinship=K.LOCO.192.chr.r,
                                         pheno=expr.rna.192[mice.with.data.r,index.rna],
                                         addcovar=addcovar.r[,-1],
                                         cores=cores.num, reml=TRUE)
    
      
# Scan / RNA / Interactive covariate: Diet
lod.rna.dietint <- scan1(genoprobs=probs[mice.with.data.r,],
                                kinship=K.LOCO.192.chr.r,
                                pheno=expr.rna.192[mice.with.data.r,index.rna],
                                addcovar=addcovar.r[,-1],
                                intcovar=intcovar.diet.r[,-1],
                                cores=cores.num, reml=TRUE)
    
# Scan / RNA / Interactive covariate: Sex
lod.rna.sexint <- scan1(genoprobs=probs[mice.with.data.r,],
                              kinship=K.LOCO.192.chr.r,
                              pheno=expr.rna.192[mice.with.data.r,index.rna],
                              addcovar=addcovar.r[,-1],
                              intcovar=intcovar.sex.r[,-1],
                              cores=cores.num, reml=TRUE)
```

``` r
      # PROTEIN
      # Scan / Protein / Additive covariate: Sex & Diet
      lod.protein.addcovariatesonly <- scan1(genoprobs=probs[mice.with.data.p,], 
                                              kinship=K.LOCO.192.chr.p,
                                              pheno=expr.protein.192[mice.with.data.p,index.protein],
                                              addcovar=addcovar.p[,-1],
                                              cores=cores.num, reml=TRUE)
      
    
      
      # Scan / Protein / Interactive covariate: Diet
      lod.protein.dietint <- scan1(genoprobs=probs[mice.with.data.p,], 
                                    kinship=K.LOCO.192.chr.p,
                                    pheno=expr.protein.192[mice.with.data.p,index.protein],
                                    addcovar=addcovar.p[,-1],
                                    intcovar=intcovar.diet.p[,-1],
                                    cores=cores.num, reml=TRUE)
      
      
      # Scan / Protein / Interactive covariate: Sex
      lod.protein.sexint <- scan1(genoprobs=probs[mice.with.data.p,], 
                                  kinship=K.LOCO.192.chr.p, 
                                  pheno=expr.protein.192[mice.with.data.p,index.protein], 
                                  addcovar=addcovar.p[,-1], 
                                  intcovar=intcovar.sex.p[,-1], 
                                  cores=cores.num, reml=TRUE)
```

Name and save the file
----------------------

Here we will name the file with identifying information. The working directory will be set so that the scan files will be saved in the 'Output' data folder.

``` r
# Systematically name the file.
file_name <- paste(annotations.protein.192[index.protein,]$Associated.Gene.Name, "_", annotations.protein.192[index.protein,]$Ensembl.Protein.ID, ".DO192Liver.DietSexInts.RData", sep="")

# Set the working directory (FOR THIS NOTEBOOK ONLY)
setwd("/Users/c-allanb/Desktop/Benjamin_Allan-Rahill/D0192Liver/Data/Output")
```

    ## Error in setwd("/Users/c-allanb/Desktop/Benjamin_Allan-Rahill/D0192Liver/Data/Output"): cannot change working directory

``` r
# Save the scan data to the file.
save(index.rna, index.protein, lod.rna.addcovariatesonly, lod.rna.dietint, lod.rna.sexint, lod.protein.addcovariatesonly, lod.protein.dietint, lod.protein.sexint, file = file_name)
```

Plot
----

### Save plot files to a specific folder

This script currently creates a seperate folder for each gene. This can be easily rewritten to put the plots into 'RNA' and 'Protein' folders.

We will name the plot file with identifying information from the RData file. we will then create a new file for the plots of that gene.

``` r
plot_file_name <- paste(gene_name, "_", annotations.protein.192[index.protein,]$Ensembl.Protein.ID, ".DO192Liver.DietSexInts", "_plot", sep="")

setwd(paste(path, "Data/Output/Plots", sep= ""))

# Make a new folder for the gene and its plots
dir.create(paste(path, "Data/Output/Plots/", gene_name, sep ="" ))
```

### RNA

This chunk set the working directory to the folder of the gene. It then greates a PNG with 3 facets, one for each scan. The y-axis are a set scale (0, 30) and each plot has three threshold lines.

-   Red = 4
-   Yellow = 7
-   Black = 10

``` r
# Set the working directory to put the plots in
setwd(paste(path, "Data/Output/Plots/", gene_name, sep ="" ))
```

    ## Error in setwd(paste(path, "Data/Output/Plots/", gene_name, sep = "")): cannot change working directory

``` r
# Set the device to PNG so that R writes .png files when it creates the plots 
png(paste(gene_name, "_RNA", sep =""), height = 780, width = 780 )

# Set the number of rows to 3 and columns to 1
par(mfrow= c(3,1))

# Create 3 plots with 3 threshold lines

plot(x = lod.rna.addcovariatesonly, map = map, main = paste(gene_name, "~Sex + Diet RNA"), ylim = c(0,30))
abline(h = 4, col = "red", lwd = 2)
abline(h = 7, col = "yellow", lwd = 2)
abline(h = 10, col = "black", lwd = 2)

plot(x = lod.rna.dietint, map = map, main = paste(gene_name, "Diet Int RNA"), ylim = c(0,30))
abline(h = 4, col = "red", lwd = 2)
abline(h = 7, col = "yellow", lwd = 2)

plot(x = lod.rna.sexint, map = map, main = paste(gene_name, "Sex Int RNA"), ylim = c(0,30))
abline(h = 4, col = "red", lwd = 2)
abline(h = 7, col = "yellow", lwd = 2)
abline(h = 10, col = "black", lwd = 2)

dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Protein

This chunk set the working directory to the folder of the gene. It then greates a PNG with 3 facets, one for each scan. The y-axis are a set scale (0, 40) and each plot has three threshold lines.

``` r
# Set the working directory to put the plots in
setwd(paste(path, "Data/Output/Plots/", gene_name, sep ="" ))
```

    ## Error in setwd(paste(path, "Data/Output/Plots/", gene_name, sep = "")): cannot change working directory

``` r
# Set the device to PNG so that R writes .png files when it creates the plots 
png(paste(path, "Data/Output/Plots/", gene_name, "_Protein", sep ="" ))

# Set the number of rows to 3 and columns to 1
par(mfrow= c(3,1))

# Create 3 plots with 3 threshold lines with set y-axis scale 

plot(x = lod.protein.addcovariatesonly, map = map, main = paste(gene_name, "~Sex + Diet Protein"), ylim = c(0,40))
abline(h = 4, col = "red", lwd = 2)
abline(h = 7, col = "yellow", lwd = 2)
abline(h = 10, col = "black", lwd = 2)

plot(x = lod.protein.dietint, map = map, main = paste(gene_name, "Diet Protein"), ylim = c(0,40))
abline(h = 4, col = "red", lwd = 2)
abline(h = 7, col = "yellow", lwd = 2)
abline(h = 10, col = "black", lwd = 2)

plot(x = lod.protein.sexint, map = map, main = paste(gene_name, "Sex Protein"), ylim = c(0,40))
abline(h = 4, col = "red", lwd = 2)
abline(h = 7, col = "yellow", lwd = 2)
abline(h = 10, col = "black", lwd = 2)

dev.off()
```

    ## quartz_off_screen 
    ##                 2
