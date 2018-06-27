### Scan, Collate, and Analyze
### Function definitions
### June 13th, 2018
### Doug Perkins

################################################################################################################################################
# ARGUMENTS
# make.data.generic

# scan.generic
# rna=T, protein=F - scans for eQTLs only
# rna=F, protein=T - scans for pQTLs only
# rna=T, protein=T, rna.with.proteins.only=T - scans for eQTLs that affect genes that also have expression data.
# rna=T, protein=T, rna.with.proteins.only=F - nothing yet. probably not worth making.

################################################################################################################################################

# Install packages
# This is almost all from the Manhattan/AE plot script.
# I'm not sure that it's all necessary for the current version of this script.
#install.packages("devtools")
#library(devtools)
#install_github("dmgatti/DOQTL")
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Mmusculus.UCSC.mm10")
#install_github("simecek/intermediate")
#install_github("kbroman/qtlcharts")
#install.packages("tidyverse")

# Load packages
library(qtl2)
#library(qtlcharts)
#library(DOQTL)
#library(intermediate)
library(dplyr)
library(ggplot2)
#The following aren't necessary yet - not sure which script used this.
#library(qtl2db)

options(stringsAsFactors = F)

################################################################################################################################################
# First rename the columns that are common to all annotations matrices, etc.
# Then if/else renaming depending on whether or not the column you're renaming is actually in the dataset.
# Ensuring that chromosome 20 is always 20 and not X? It seems to be X even in the data that's formatted for qtl2.
# I might need to do something with the covariates dataframes. I might need to know more about them to know what to name the generic columns, or I can just infer/name them based on how they're used in the other scripts.


make.data.generic <- function(rna.annotations, rna.expression, rna.covariates, 
                              rna.ENSMUSG.col, rna.ENSMUST.col, rna.ENSMUSP.col,
                              rna.gene.col, rna.chromosome.col, rna.start.mbp.col, rna.end.mbp.col, rna.gene.biotype.col,
                              
                              protein.annotations, protein.expression, protein.covariates, 
                              protein.ENSMUSG.col, protein.ENSMUST.col, protein.ENSMUSP.col,
                              protein.gene.col, protein.chromosome.col, protein.start.mbp.col, protein.end.mbp.col, protein.gene.biotype.col,
                              
                              kinship, founder.probs, 
                              
                              markers,
                              markers.id.col, markers.chromosome.col, markers.mbp.col, markers.cM.col, markers.bp.col){
  
  print(paste("rna.ENSMUSG.col: ", rna.ENSMUSG.col))
  print(paste("rna.gene.col: ", rna.gene.col))
  print(paste("rna.chromosome.col: ", rna.chromosome.col))
  print(paste("rna.start.mbp.col: ", rna.start.mbp.col))
  print(paste("rna.end.mbp.col: ", rna.end.mbp.col))
  print(paste("rna.gene.biotype.col: ", rna.gene.biotype.col))
  
  # Ideas of what the problem is: try rna.annotations %>% rename(), check the column names on the HPC version of the data, try 
  rna.annotations <- rna.annotations %>% select(
                           ENSMUSG = rna.ENSMUSG.col,
                           Gene = rna.gene.col,
                           Chr = rna.chromosome.col,
                           Mbp.Start = rna.start.mbp.col,
                           Mbp.End = rna.end.mbp.col,
                           Gene.Biotype = rna.gene.biotype.col)
  
  protein.annotations <- protein.annotations %>% select(
                               ENSMUSG = protein.ENSMUSG.col,
                               ENSMUSP = protein.ENSMUSP.col,
                               Gene = protein.gene.col,
                               Chr = protein.chromosome.col,
                               Mbp.Start = protein.start.mbp.col,
                               Mbp.End = protein.end.mbp.col,
                               Gene.Biotype = protein.gene.biotype.col)

  markers <- markers %>% select(
                   Marker = markers.id.col,
                   Chr = markers.chromosome.col,
                   Mbp = markers.mbp.col,
                   cM = markers.cM.col,
                   bp = markers.bp.col
                   )
  
  # Change any X or x chromosomes to "20"
  if (any(rna.annotations$Chr=='X')){
  rna.annotations[rna.annotations$Chr=='X',]$Chr <- '20'
  }
  
  if (any(protein.annotations$Chr=='X')){
  protein.annotations[protein.annotations$Chr=='X',]$Chr <- '20'
  }
  
  if (any(markers$Chr=='X')){
  markers[markers$Chr=='X',]$Chr <- '20'
  }
  
  if (any(rna.annotations$Chr=='x')){
    rna.annotations[rna.annotations$Chr=='x',]$Chr <- '20'
  }
  
  if (any(protein.annotations$Chr=='x')){
    protein.annotations[protein.annotations$Chr=='x',]$Chr <- '20'
  }
  
  if (any(markers$Chr=='x')){
    markers[markers$Chr=='x',]$Chr <- '20'
  }
  
  list(rna.annotations, rna.expression, rna.covariates, protein.annotations, protein.expression, protein.covariates, kinship, founder.probs, markers)
  
}

################################################################################################################################################  
# Do not edit these. They load in the new generic data frames.
make.generic.data.structures <- function(){
  rna.annotations <- generic.data[[1]]
  rna.expression <- generic.data[[2]]
  rna.covariates <- generic.data[[3]]
  protein.annotations <- generic.data[[4]]
  protein.expression <- generic.data[[5]]
  protein.covariates <- generic.data[[6]]
  kinship <- generic.data[[7]]
  founder.probs <- generic.data[[8]]
  markers <- generic.data[[9]]
  
  assign("rna.annotations", rna.annotations, envir = .GlobalEnv)
  assign("rna.expression", rna.expression, envir = .GlobalEnv)
  assign("rna.covariates", rna.covariates, envir = .GlobalEnv)
  assign("protein.annotations", protein.annotations, envir = .GlobalEnv)
  assign("protein.expression", protein.expression, envir = .GlobalEnv)
  assign("protein.covariates", protein.covariates, envir = .GlobalEnv)
  assign("kinship", kinship, envir = .GlobalEnv)
  assign("founder.probs", founder.probs, envir = .GlobalEnv)
  assign("markers", markers, envir = .GlobalEnv)
  
}
################################################################################################################################################  

# Maybe make this two smaller functions, rna.scan.generic and protein.scan.generic? Then feed it the nrows of each annotations file as a default num.samples.
scan.generic <- function(rna=T, protein=F, phos=F,
                         rna.beg=1, rna.end=1,
                         protein.beg=1, protein.end=1,
                         rna.with.proteins.only=F){
  
  # Set the number of processor cores to be used. Use 1 on the HPC.
  cores.num = 1
  
  if (test==T){
    for (j in rna.beg:rna.end) {
      # Get the gene's index in the rna data.
      index.rna <- j  
      
      # Get the chromosome that the gene is on
      gene_chr <- rna.annotations[index.rna,]$Chr  
      
      # Determine the SNP closest to the center of the gene.
      closestSNP <- markers[markers$Chr==gene_chr,][which.min(abs(mean(rna.annotations[index.rna,]$Mbp.Start, rna.annotations[index.rna,]$Mbp.End) - markers[which(markers$Chr==gene_chr),]$bp)),]
      
      # Only subset mice with protein expression data. I don't remember what issue this prevented.
      # mice.with.data.r <- !is.na(rna.expression[,index.rna])
      
      # Set up kinship matrices for the mice that have data for your given RNA/protein. Does it need to be for all chromosomes?
      kinship.chr.r <- kinship
      
      #for (chromosome in 1:20){
      #  kinship.chr.r[[chromosome]] <- kinship[[chromosome]][mice.with.data.r,mice.with.data.r]
      #}
      
      # Set up additive and interactive covariates for RNA/protein, only mice with data.
      addcovar.r <- model.matrix(~ Sex + Diet, data=rna.covariates[, ])
      intcovar.sex.r <- model.matrix(~ Sex, data=rna.covariates[, ])
      intcovar.diet.r <- model.matrix(~ Diet, data=rna.covariates[, ])
      
      # Perform QTL/LOD scans for RNA/protein in conditions of additive, interactive diet, and interactive sex covariates.
      cat("Scanning ", j, " out of ", nrow(rna.annotations), " ",rna.annotations$Gene[j], "\n")
      # RNA
      # Scan / RNA / Additive covariate: Sex & Diet
      lod.rna.addcovariatesonly <- scan1(genoprobs=probs[,],
                                         kinship=kinship.chr.r,
                                         pheno=rna.expression[,index.rna],
                                         addcovar=addcovar.r[,-1],
                                         cores=cores.num, reml=TRUE)
      
      
      # Scan / RNA / Interactive covariate: Diet
      lod.rna.dietint <- scan1(genoprobs=probs[,],
                               kinship=kinship.chr.r,
                               pheno=rna.expression[,index.rna],
                               addcovar=addcovar.r[,-1],
                               intcovar=intcovar.diet.r[,-1],
                               cores=cores.num, reml=TRUE)
      
      
      # Scan / RNA / Interactive covariate: Sex
      lod.rna.sexint <- scan1(genoprobs=probs[,],
                              kinship=kinship.chr.r,
                              pheno=rna.expression[,index.rna],
                              addcovar=addcovar.r[,-1],
                              intcovar=intcovar.sex.r[,-1],
                              cores=cores.num, reml=TRUE)
      
      # Systematically name the file.
      file_name <- paste(rna.annotations[index.rna,]$Gene, "_", rna.annotations[index.rna,]$ENSMUSG, ".RNA.DO192Liver.DietSexIntScan.RData", sep="")
      
      # Save the scan data to the file.
      setwd(dir.output)
      if (dir.exists("rna_scan")){
        setwd("rna_scan")
      }
      else{
        dir.create(file.path(dir.output, "rna_scan"))
        setwd("rna_scan")
      }
      save(index.rna, closestSNP, lod.rna.addcovariatesonly, lod.rna.dietint, lod.rna.sexint, file = file_name)
      
      # Display progress info by count and gene name.
      cat("Done scanning ", j, " out of ", nrow(rna.annotations), " ",rna.annotations$Gene[j], "\n")
    }
  }
  
    if (rna==T & rna.with.proteins.only==F){
      for (j in rna.beg:rna.end) {
      # Get the gene's index in the rna data.
      index.rna <- j  
      
      # Get the chromosome that the gene is on
      gene_chr <- rna.annotations[index.rna,]$Chr  
      
      # Determine the SNP closest to the center of the gene.
      closestSNP <- markers[markers$Chr==gene_chr,][which.min(abs(mean(rna.annotations[index.rna,]$Mbp.Start, rna.annotations[index.rna,]$Mbp.End) - markers[which(markers$Chr==gene_chr),]$bp)),]
      
      # Only subset mice with protein expression data. I don't remember what issue this prevented.
      mice.with.data.r <- !is.na(rna.expression[,index.rna])
      
      # Set up kinship matrices for the mice that have data for your given RNA/protein. Does it need to be for all chromosomes?
      kinship.chr.r <- kinship

      for (chromosome in 1:20){
      kinship.chr.r[[chromosome]] <- kinship[[chromosome]][mice.with.data.r,mice.with.data.r]
      }
      
      # Set up additive and interactive covariates for RNA/protein, only mice with data.
      addcovar.r <- model.matrix(~ Sex + Diet, data=rna.covariates[mice.with.data.r, ])
      intcovar.sex.r <- model.matrix(~ Sex, data=rna.covariates[mice.with.data.r, ])
      intcovar.diet.r <- model.matrix(~ Diet, data=rna.covariates[mice.with.data.r, ])
      
      # Perform QTL/LOD scans for RNA/protein in conditions of additive, interactive diet, and interactive sex covariates.
      cat("Scanning ", j, " out of ", nrow(rna.annotations), " ",rna.annotations$Gene[j], "\n")
      # RNA
      # Scan / RNA / Additive covariate: Sex & Diet
      lod.rna.addcovariatesonly <- scan1(genoprobs=probs[mice.with.data.r,],
                                         kinship=kinship.chr.r,
                                         pheno=rna.expression[mice.with.data.r,index.rna],
                                         addcovar=addcovar.r[,-1],
                                         cores=cores.num, reml=TRUE)
      
      
      # Scan / RNA / Interactive covariate: Diet
      lod.rna.dietint <- scan1(genoprobs=probs[mice.with.data.r,],
                               kinship=kinship.chr.r,
                               pheno=rna.expression[mice.with.data.r,index.rna],
                               addcovar=addcovar.r[,-1],
                               intcovar=intcovar.diet.r[,-1],
                               cores=cores.num, reml=TRUE)
      
      
      # Scan / RNA / Interactive covariate: Sex
      lod.rna.sexint <- scan1(genoprobs=probs[mice.with.data.r,],
                              kinship=kinship.chr.r,
                              pheno=rna.expression[mice.with.data.r,index.rna],
                              addcovar=addcovar.r[,-1],
                              intcovar=intcovar.sex.r[,-1],
                              cores=cores.num, reml=TRUE)
      
      # Systematically name the file.
      file_name <- paste(rna.annotations[index.rna,]$Gene, "_", rna.annotations[index.rna,]$ENSMUSG, ".RNA.DO192Liver.DietSexIntScan.RData", sep="")
      
      # Save the scan data to the file.
      setwd(dir.output)
      if (dir.exists("rna_scan")){
        setwd("rna_scan")
      }
      else{
        dir.create(file.path(dir.output, "rna_scan"))
        setwd("rna_scan")
      }
      save(index.rna, closestSNP, lod.rna.addcovariatesonly, lod.rna.dietint, lod.rna.sexint, file = file_name)
      
      # Display progress info by count and gene name.
      cat("Done scanning ", j, " out of ", nrow(rna.annotations), " ",rna.annotations$Gene[j], "\n")
      }
    }
    
    if (protein==T & rna.with.proteins.only==F){
      for (j in protein.beg:protein.end) {
      # Get the gene's index in the protein data.
      index.protein <- j

      # Get the chromosome that the gene is on
      gene_chr <- protein.annotations[index.protein,]$Chr  
      
      # Determine the SNP closest to the center of the gene.
      closestSNP <- markers[markers$Chr==gene_chr,][which.min(abs(mean(protein.annotations[index.protein,]$Mbp.Start, protein.annotations[index.protein,]$Mbp.End) - markers[which(markers$Chr==gene_chr),]$bp)),]
      
      # Only subset mice with protein expression data. I don't remember what issue this prevented.
      mice.with.data.p <- !is.na(protein.expression[,index.protein])
      
      # Set up kinship matrices for the mice that have data for your given RNA/protein. Does it need to be for all chromosomes?
      kinship.chr.p <- kinship
      for (chromosome in 1:20){
        kinship.chr.p[[chromosome]] <- kinship[[chromosome]][mice.with.data.p,mice.with.data.p]
      }  
      
      addcovar.p <- model.matrix(~ Sex + Diet + Tag, data=protein.covariates[mice.with.data.p, ])
      intcovar.sex.p <- model.matrix(~ Sex, data=protein.covariates[mice.with.data.p, ])
      intcovar.diet.p <- model.matrix(~ Diet, data=protein.covariates[mice.with.data.p, ])
      
      # Perform QTL/LOD scans for RNA/protein in conditions of additive, interactive diet, and interactive sex covariates.
      cat("Scanning ", j, " out of ", nrow(protein.annotations), " ",protein.annotations$Gene[j], "\n")
      
      # PROTEIN
      # Scan / Protein / Additive covariate: Sex & Diet
      lod.protein.addcovariatesonly <- scan1(genoprobs=probs[mice.with.data.p,], 
                                             kinship=kinship.chr.p,
                                             pheno=protein.expression[mice.with.data.p,index.protein],
                                             addcovar=addcovar.p[,-1],
                                             cores=cores.num, reml=TRUE)
      
      # Scan / Protein / Interactive covariate: Diet
      lod.protein.dietint <- scan1(genoprobs=probs[mice.with.data.p,], 
                                   kinship=kinship.chr.p,
                                   pheno=protein.expression[mice.with.data.p,index.protein],
                                   addcovar=addcovar.p[,-1],
                                   intcovar=intcovar.diet.p[,-1],
                                   cores=cores.num, reml=TRUE)
      
      # Scan / Protein / Interactive covariate: Sex
      lod.protein.sexint <- scan1(genoprobs=probs[mice.with.data.p,], 
                                  kinship=kinship.chr.p, 
                                  pheno=protein.expression[mice.with.data.p,index.protein], 
                                  addcovar=addcovar.p[,-1], 
                                  intcovar=intcovar.sex.p[,-1], 
                                  cores=cores.num, reml=TRUE)
      
      # Systematically name the file.
      file_name <- paste(protein.annotations[index.protein,]$Gene, "_", protein.annotations[index.protein,]$ENSMUSP, ".Protein.DO192Liver.DietSexIntScan.RData", sep="")
      
      # Save the scan data to the file.
      setwd(dir.output)
      if (dir.exists("protein_scan")){
        setwd("protein_scan")
      }
      else{
        dir.create(file.path(dir.output, "protein_scan"))
        setwd("protein_scan")
      }
      save(index.protein, closestSNP, lod.protein.addcovariatesonly, lod.protein.dietint, lod.protein.sexint, file = file_name)
      
      # Display progress info by count and gene name.
      cat("Done scanning ", j, " out of ", nrow(protein.annotations), " ", protein.annotations$Gene[j], "\n")
      }
    }
    
    if (rna==T & protein==T & rna.with.proteins.only==T){
      for (j in protein.beg:protein.end) {
    # Get the gene's index in the protein and rna data.
    index.protein <- j
    index.rna <- which(rna.annotations$ENSMUSG == protein.annotations[index.protein,]$ENSMUSG)  
    
    # Get the chromosome that the gene is on
    gene_chr <- protein.annotations[index.protein,]$Chr  
    
    # Determine the SNP closest to the center of the gene.
    closestSNP <- markers[markers$Chr==gene_chr,][which.min(abs(mean(protein.annotations[index.protein,]$Mbp.Start, protein.annotations[index.protein,]$Mbp.End) - markers[which(markers$Chr==gene_chr),]$bp)),]
    
    # Only subset mice with protein expression data. I don't remember what issue this prevented.
    mice.with.data.r <- !is.na(rna.expression[,index.rna])
    mice.with.data.p <- !is.na(protein.expression[,index.protein])
    
    # Set up kinship matrices for the mice that have data for your given RNA/protein. Does it need to be for all chromosomes?
    kinship.chr.r <- kinship
    kinship.chr.p <- kinship
    for (chromosome in 1:20){
      kinship.chr.r[[chromosome]] <- kinship[[chromosome]][mice.with.data.r,mice.with.data.r]
    }
    for (chromosome in 1:20){
      kinship.chr.p[[chromosome]] <- kinship[[chromosome]][mice.with.data.p,mice.with.data.p]
    }  
    
    # Set up additive and interactive covariates for RNA/protein, only mice with data.
    addcovar.r <- model.matrix(~ Sex + Diet, data=rna.covariates[mice.with.data.r, ])
    intcovar.sex.r <- model.matrix(~ Sex, data=rna.covariates[mice.with.data.r, ])
    intcovar.diet.r <- model.matrix(~ Diet, data=rna.covariates[mice.with.data.r, ])
    
    addcovar.p <- model.matrix(~ Sex + Diet + Tag, data=protein.covariates[mice.with.data.p, ])
    intcovar.sex.p <- model.matrix(~ Sex, data=protein.covariates[mice.with.data.p, ])
    intcovar.diet.p <- model.matrix(~ Diet, data=protein.covariates[mice.with.data.p, ])
    
    # Perform QTL/LOD scans for RNA/protein in conditions of additive, interactive diet, and interactive sex covariates.
    cat("Scanning ", j, " out of ", nrow(protein.annotations), " ",protein.annotations$Gene[j], "\n")
    # RNA
    # Scan / RNA / Additive covariate: Sex & Diet
    lod.rna.addcovariatesonly <- scan1(genoprobs=probs[mice.with.data.r,],
                                       kinship=kinship.chr.r,
                                       pheno=rna.expression[mice.with.data.r,index.rna],
                                       addcovar=addcovar.r[,-1],
                                       cores=cores.num, reml=TRUE)
    
    
    # Scan / RNA / Interactive covariate: Diet
    lod.rna.dietint <- scan1(genoprobs=probs[mice.with.data.r,],
                             kinship=kinship.chr.r,
                             pheno=rna.expression[mice.with.data.r,index.rna],
                             addcovar=addcovar.r[,-1],
                             intcovar=intcovar.diet.r[,-1],
                             cores=cores.num, reml=TRUE)
    
    
    # Scan / RNA / Interactive covariate: Sex
    lod.rna.sexint <- scan1(genoprobs=probs[mice.with.data.r,],
                            kinship=kinship.chr.r,
                            pheno=rna.expression[mice.with.data.r,index.rna],
                            addcovar=addcovar.r[,-1],
                            intcovar=intcovar.sex.r[,-1],
                            cores=cores.num, reml=TRUE)
    
    # PROTEIN
    # Scan / Protein / Additive covariate: Sex & Diet
    lod.protein.addcovariatesonly <- scan1(genoprobs=probs[mice.with.data.p,], 
                                           kinship=kinship.chr.p,
                                           pheno=protein.expression[mice.with.data.p,index.protein],
                                           addcovar=addcovar.p[,-1],
                                           cores=cores.num, reml=TRUE)
    
    
    
    # Scan / Protein / Interactive covariate: Diet
    lod.protein.dietint <- scan1(genoprobs=probs[mice.with.data.p,], 
                                 kinship=kinship.chr.p,
                                 pheno=protein.expression[mice.with.data.p,index.protein],
                                 addcovar=addcovar.p[,-1],
                                 intcovar=intcovar.diet.p[,-1],
                                 cores=cores.num, reml=TRUE)
    
    
    # Scan / Protein / Interactive covariate: Sex
    lod.protein.sexint <- scan1(genoprobs=probs[mice.with.data.p,], 
                                kinship=kinship.chr.p, 
                                pheno=protein.expression[mice.with.data.p,index.protein], 
                                addcovar=addcovar.p[,-1], 
                                intcovar=intcovar.sex.p[,-1], 
                                cores=cores.num, reml=TRUE)
    
    # Systematically name the file.
    file_name <- paste(protein.annotations[index.protein,]$Gene, "_", protein.annotations[index.protein,]$ENSMUSP, "RNA.Protein.DO192Liver.DietSexIntScan.RData", sep="")
    
    # Save the scan data to the file.
    setwd(dir.output)
    if (dir.exists("rna_protein_scan")){
      setwd("rna_protein_scan")
    }
    else{
      dir.create(file.path(dir.output, "rna_protein_scan"))
      setwd("rna_protein_scan")
    }    
    save(index.rna, index.protein, closestSNP, lod.rna.addcovariatesonly, lod.rna.dietint, lod.rna.sexint, lod.protein.addcovariatesonly, lod.protein.dietint, lod.protein.sexint, file = file_name)
    
    # Display progress info by count and gene name.
    cat("Done scanning ", j, " out of ", nrow(protein.annotations), " ",protein.annotations$Gene[j], "\n")
      }
  }
  
}

################################################################################################################################################

collate.generic <- function(rna=T, protein=F, phos=F,
                            rna.diet.int=T, rna.sex.int=T, protein.diet.int=T, protein.sex.int=T,
                            lod.thresholds=c(4,5,7), window.size.mbp=5, window.type="center", peaks.only=T){

  setwd(dir.output)
  
  if (rna==T & protein==F){
    if (dir.exists("rna_scan")){
      setwd("rna_scan")
    }
    else{
      dir.create(dir.output.rna.scan)
      setwd("rna_scan")
    }
    dir.output.rna.scan <- file.path(dir.output, "rna_scan")
    dir.files <- dir.output.rna.scan
    filenames <- dir(path = dir.output.rna.scan, pattern = "..DietSexIntScan.RData")
    setwd(dir.output.rna.scan)
  }
  if (rna==F & protein==T){
    if (dir.exists("protein_scan")){
      setwd("protein_scan")
    }
    else{
      dir.create(dir.output.protein.scan)
      setwd("protein_scan")
    }
    dir.output.protein.scan <- file.path(dir.output, "protein_scan")
    dir.files <- dir.output.protein.scan
    filenames <- dir(path = dir.output.protein.scan, pattern = "..DietSexIntScan.RData")
    setwd(dir.output.protein.scan)
  }
  if (rna==T & protein==T){
    if (dir.exists("rna_protein_scan")){
      setwd("rna_protein_scan")
    }
    else{
      dir.create(dir.output.rna.protein.scan)
      setwd("rna_protein_scan")
    }
    dir.output.rna.protein.scan <- file.path(dir.output, "rna_protein_scan")
    dir.files <- dir.output.rna.protein.scan
    filenames <- dir(path = dir.output.rna.protein.scan, pattern = "..DietSexIntScan.RData")
    setwd(dir.output.rna.protein.scan)
  }

for (threshold in lod.thresholds){
  # Initialize empty table.
  print(paste("Threshold: ", threshold))
  tmp = rep(0,0)
  
  if (rna == T & protein == F){
    if (rna.diet.int == T){
      results.rna.diet.final = data.frame(ENSMUSG=tmp, Gene=tmp, Chr=tmp,
                                          Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                          SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                          r.LOD.Add=tmp,r.LOD.Diet.Int=tmp,r.LOD.Sex.Int=tmp,r.LOD.Diet.Diff=tmp,r.LOD.Sex.Diff=tmp)
    }
  
    if (rna.sex.int == T){  
      results.rna.sex.final = data.frame(ENSMUSG=tmp, Gene=tmp,Chr=tmp,
                                         Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                         SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                         r.LOD.Add=tmp,r.LOD.Diet.Int=tmp,r.LOD.Sex.Int=tmp,r.LOD.Diet.Diff=tmp,r.LOD.Sex.Diff=tmp)
    }
  }
  
  else if (rna == F & protein == T){
    if (protein.diet.int == T){
      results.protein.diet.final = data.frame(ENSMUSP=tmp, ENSMUSG=tmp, Gene=tmp,Chr=tmp,
                                              Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                              SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                              p.LOD.Add=tmp,p.LOD.Diet.Int=tmp,p.LOD.Sex.Int=tmp,p.LOD.Diet.Diff=tmp,p.LOD.Sex.Diff=tmp)
    }
    if (protein.sex.int == T){
      results.protein.sex.final = data.frame(ENSMUSP=tmp, ENSMUSG=tmp, Gene=tmp,Chr=tmp,
                                             Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                             SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                             p.LOD.Add=tmp,p.LOD.Diet.Int=tmp,p.LOD.Sex.Int=tmp,p.LOD.Diet.Diff=tmp,p.LOD.Sex.Diff=tmp)
    }
  }
  
  else if (rna == T & protein == T){
    if (rna.diet.int == T){
      results.rna.diet.final = data.frame(ENSMUSG=tmp, Gene=tmp, Chr=tmp,
                                          Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                          SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                          r.LOD.Add=tmp,r.LOD.Diet.Int=tmp,r.LOD.Sex.Int=tmp,r.LOD.Diet.Diff=tmp,r.LOD.Sex.Diff=tmp,
                                          p.LOD.Add=tmp,p.LOD.Diet.Int=tmp,p.LOD.Sex.Int=tmp,p.LOD.Diet.Diff=tmp,p.LOD.Sex.Diff=tmp)
                                          
    }
    if (rna.sex.int == T){
      results.rna.sex.final = data.frame(ENSMUSG=tmp, Gene=tmp,Chr=tmp,
                                         Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                         SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                         r.LOD.Add=tmp,r.LOD.Diet.Int=tmp,r.LOD.Sex.Int=tmp,r.LOD.Diet.Diff=tmp,r.LOD.Sex.Diff=tmp,
                                         p.LOD.Add=tmp,p.LOD.Diet.Int=tmp,p.LOD.Sex.Int=tmp,p.LOD.Diet.Diff=tmp,p.LOD.Sex.Diff=tmp)
    }
    if (protein.diet.int == T){
      results.protein.diet.final = data.frame(ENSMUSP=tmp, ENSMUSG=tmp, Gene=tmp,Chr=tmp,
                                              Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                              SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                              p.LOD.Add=tmp,p.LOD.Diet.Int=tmp,p.LOD.Sex.Int=tmp,p.LOD.Diet.Diff=tmp,p.LOD.Sex.Diff=tmp,
                                              r.LOD.Add=tmp,r.LOD.Diet.Int=tmp,r.LOD.Sex.Int=tmp,r.LOD.Diet.Diff=tmp,r.LOD.Sex.Diff=tmp)
    }
    if (protein.sex.int == T){
      results.protein.sex.final = data.frame(ENSMUSP=tmp, ENSMUSG=tmp, Gene=tmp,Chr=tmp,
                                             Mbp.Start=tmp, Mbp.End=tmp, n.Expressed=tmp,
                                             SNP.ID=tmp,SNP.Chr=tmp,SNP.Mbp=tmp,SNP.cM=tmp, QTL.Type=tmp, QTL.Distance=tmp,
                                             p.LOD.Add=tmp,p.LOD.Diet.Int=tmp,p.LOD.Sex.Int=tmp,p.LOD.Diet.Diff=tmp,p.LOD.Sex.Diff=tmp,
                                             r.LOD.Add=tmp,r.LOD.Diet.Int=tmp,r.LOD.Sex.Int=tmp,r.LOD.Diet.Diff=tmp,r.LOD.Sex.Diff=tmp)
    }
  }
  
  
  for (f in 1:length(filenames)){
    setwd(dir.files)
    load(filenames[f])
    map.peakfinder =  map_df_to_list(map = markers, marker_column = "Marker", chr_column = "Chr", pos_column = "bp")

    if (rna==T){
    num.expressed.r <- length(which(!is.na(rna.expression[,index.rna])))
    }
    if (protein==T){
    num.expressed.p <- length(which(!is.na(protein.expression[,index.protein])))
    }

    # Display progress by count and gene name to the user
    print(paste("Processing", f,"of",length(filenames), rna.annotations$Gene[index.rna], sep=" "))
    
    if (rna==T & protein==F){
      
      if (rna.diet.int==T){
        print("eQTL/Diet Table Processing")
        print("find_peaks started")
        qtl2.peakfinder <- find_peaks(lod.rna.dietint - lod.rna.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)
        print("find_peaks completed")

            for (qtl.marker in rownames(lod.rna.dietint)){
              print(paste("qtl.marker: ", qtl.marker))
              
              if ((lod.rna.dietint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){
                print(paste("Marker", qtl.marker, "above threshold"), sep="")
                qtl.marker.data <- markers[markers$Marker==qtl.marker,]
                
                # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
                if (window.type=="center"){
                  gene.center <- abs(rna.annotations$Mbp.Start[index.rna] - rna.annotations$Mbp.End[index.rna]) / 2
                  qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
                }
                else if (window.type=="edge"){
                  gene.left.edge <- rna.annotations$Mbp.Start[index.rna]
                  gene.right.edge <- rna.annotations$Mbp.End[index.rna]
                  if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                    closer.edge <- gene.left.edge
                  }
                  else {
                    closer.edge <- gene.right.edge
                  }
                  qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
                }
                
                # Qualify QTL as either local or distant depending on the selected window size.
                if (qtl.distance > window.size.mbp/2){
                  qtl.type <- "distant"
                }
                else {
                  qtl.type <- "local"
                }
                
                this.lod.add.r <- round(lod.rna.addcovariatesonly[qtl.marker,], digits=2)
                this.lod.diet.int.r <- round(lod.rna.dietint[qtl.marker,], digits=2)
                this.lod.sex.int.r <- round(lod.rna.sexint[qtl.marker,], digits=2)
                this.lod.diet.diff.r <- round(this.lod.diet.int.r - this.lod.add.r, digits=2)
                this.lod.sex.diff.r <- round(this.lod.sex.int.r - this.lod.add.r, digits=2)
                
                print(paste("Adding marker to table: ", qtl.marker))
                
                if (peaks.only==F){
                  results.rna.diet.final[nrow(results.rna.diet.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                                               rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                                               qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                               this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r)
                }
  
                if (peaks.only==T){
                  for (peak in qtl2.peakfinder$pos){
                    if (peak == qtl.marker.data$Mbp*1e6){
                      results.rna.diet.final[nrow(results.rna.diet.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                                                   rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                                                   qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                   this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r)
                    }
                  }              
                }
              }
            }
          }
        
      if (rna.sex.int==T){
        print("eQTL/Sex Table Processing")
        qtl2.peakfinder <- find_peaks(lod.rna.sexint - lod.rna.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)

           for (qtl.marker in rownames(lod.rna.sexint)){
              if ((lod.rna.sexint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){
                print(paste("Marker", qtl.marker, "above threshold"))
                qtl.marker.data <- markers[markers$Marker==qtl.marker,]
                
                # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
                if (window.type=="center"){
                  gene.center <- abs(rna.annotations$Mbp.Start[index.rna] - rna.annotations$Mbp.End[index.rna]) / 2
                  qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
                }
                else if (window.type=="edge"){
                  gene.left.edge <- rna.annotations$Mbp.Start[index.rna]
                  gene.right.edge <- rna.annotations$Mbp.End[index.rna]
                  if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                    closer.edge <- gene.left.edge
                  }
                  else {
                    closer.edge <- gene.right.edge
                  }
                  qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
                }
                
                # Qualify QTL as either local or distant depending on the selected window size.
                if (qtl.distance > window.size.mbp/2){
                  qtl.type <- "distant"
                }
                else {
                  qtl.type <- "local"
                }
                
                this.lod.add.r <- lod.rna.addcovariatesonly[rownames(lod.rna.addcovariatesonly)==qtl.marker]
                this.lod.diet.int.r <- lod.rna.dietint[rownames(lod.rna.dietint)==qtl.marker]
                this.lod.sex.int.r <- lod.rna.sexint[rownames(lod.rna.sexint)==qtl.marker]
                this.lod.diet.diff.r <- this.lod.diet.int.r - this.lod.add.r
                this.lod.sex.diff.r <- this.lod.sex.int.r - this.lod.add.r
                
                print(paste("Adding marker to table: ", qtl.marker))
                
                if (peaks.only == F){
                  results.rna.sex.final[nrow(results.rna.sex.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                                             rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                                             qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                             this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r)
                }
                
                if (peaks.only == T){
                  for (peak in qtl2.peakfinder$pos){
                    if (peak == qtl.marker.data$Mbp*1e6){
                      results.rna.sex.final[nrow(results.rna.sex.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                                                 rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                                                 qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                 this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r)
                    }
                  }
                }
              } 
            }
          }
        }


    if (rna==F & protein==T){
      
      if (protein.diet.int==T){
        print("pQTL/Diet Table Processing")
        qtl2.peakfinder <- find_peaks(lod.protein.dietint - lod.protein.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)
          
          for (qtl.marker in rownames(lod.protein.dietint)){
            if ((lod.protein.dietint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){

              qtl.marker.data <- markers[markers$Marker==qtl.marker,]
              
              # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
              if (window.type=="center"){
                gene.center <- abs(protein.annotations$Mbp.Start[index.protein] - protein.annotations$Mbp.End[index.protein]) / 2
                qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
              }
              else if (window.type=="edge"){
                gene.left.edge <- protein.annotations$Mbp.Start[index.protein]
                gene.right.edge <- protein.annotations$Mbp.End[index.protein]
                if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                  closer.edge <- gene.left.edge
                }
                else {
                  closer.edge <- gene.right.edge
                }
                qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
              }
              
              # Qualify QTL as either local or distant depending on the selected window size.
              if (qtl.distance > window.size.mbp/2){
                qtl.type <- "distant"
              }
              else {
                qtl.type <- "local"
              }
              
              this.lod.add.p <- lod.protein.addcovariatesonly[rownames(lod.protein.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.p <- lod.protein.dietint[rownames(lod.protein.dietint)==qtl.marker]
              this.lod.sex.int.p <- lod.protein.sexint[rownames(lod.protein.sexint)==qtl.marker]
              this.lod.diet.diff.p <- this.lod.diet.int.p - this.lod.add.p
              this.lod.sex.diff.p <- this.lod.sex.int.p - this.lod.add.p
              
              print(paste("Adding marker to table: ", qtl.marker))
              
              if (peaks.only==F){
              results.protein.diet.final[nrow(results.protein.diet.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                  protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                  qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                  this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
              }
              
              if (peaks.only==T){
                for (peak in qtl2.peakfinder$pos){
                  if (peak == qtl.marker.data$Mbp*1e6){
                    results.protein.diet.final[nrow(results.protein.diet.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                                                         protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                                                         qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                         this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
                  }
                }
              }
            }
          }
        }
      
      if (protein.sex.int==T){
        print("pQTL/Sex Table Processing")
        qtl2.peakfinder <- find_peaks(lod.protein.sexint - lod.protein.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)
          
          for (qtl.marker in rownames(lod.protein.sexint)){
            if ((lod.protein.sexint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){
              
              qtl.marker.data <- markers[markers$Marker==qtl.marker,]
              
              # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
              if (window.type=="center"){
                gene.center <- abs(protein.annotations$Mbp.Start[index.protein] - protein.annotations$Mbp.End[index.protein]) / 2
                qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
              }
              else if (window.type=="edge"){
                gene.left.edge <- protein.annotations$Mbp.Start[index.protein]
                gene.right.edge <- protein.annotations$Mbp.End[index.protein]
                if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                  closer.edge <- gene.left.edge
                }
                else {
                  closer.edge <- gene.right.edge
                }
                qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
              }
              
              # Qualify QTL as either local or distant depending on the selected window size.
              if (qtl.distance > window.size.mbp/2){
                qtl.type <- "distant"
              }
              else {
                qtl.type <- "local"
              }
              
              this.lod.add.p <- lod.protein.addcovariatesonly[rownames(lod.protein.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.p <- lod.protein.dietint[rownames(lod.protein.dietint)==qtl.marker]
              this.lod.sex.int.p <- lod.protein.sexint[rownames(lod.protein.sexint)==qtl.marker]
              this.lod.diet.diff.p <- this.lod.diet.int.p - this.lod.add.p
              this.lod.sex.diff.p <- this.lod.sex.int.p - this.lod.add.p
              
              print(paste("Adding marker to table: ", qtl.marker))
              
              if (peaks.only==F){
              results.protein.sex.final[nrow(results.protein.sex.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                 protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                 qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                 this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
              }
              
              if (peaks.only==T){
                for (peak in qtl2.peakfinder$pos){
                  if (peak == qtl.marker.data$Mbp*1e6){
                    results.protein.sex.final[nrow(results.protein.sex.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                                                       protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                                                       qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                       this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
                  }
                }
              }
            }
          }
        }
      }
    
    
    if (rna==T & protein==T){  
      
      if (rna.diet.int==T){
        print("eQTL/Diet Table Processing")
        qtl2.peakfinder <- find_peaks(lod.rna.dietint - lod.rna.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)

          for (qtl.marker in rownames(lod.rna.dietint)){
            if ((lod.rna.dietint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){

              qtl.marker.data <- markers[markers$Marker==qtl.marker,]
              
              # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
              if (window.type=="center"){
                gene.center <- abs(rna.annotations$Mbp.Start[index.rna] - rna.annotations$Mbp.End[index.rna]) / 2
                qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
              }
              else if (window.type=="edge"){
                gene.left.edge <- rna.annotations$Mbp.Start[index.rna]
                gene.right.edge <- rna.annotations$Mbp.End[index.rna]
                if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                  closer.edge <- gene.left.edge
                }
                else {
                  closer.edge <- gene.right.edge
                }
                qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
              }
              
              # Qualify QTL as either local or distant depending on the selected window size.
              if (qtl.distance > window.size.mbp/2){
                qtl.type <- "distant"
              }
              else {
                qtl.type <- "local"
              }
              
              this.lod.add.r <- lod.rna.addcovariatesonly[rownames(lod.rna.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.r <- lod.rna.dietint[rownames(lod.rna.dietint)==qtl.marker]
              this.lod.sex.int.r <- lod.rna.sexint[rownames(lod.rna.sexint)==qtl.marker]
              this.lod.diet.diff.r <- this.lod.diet.int.r - this.lod.add.r
              this.lod.sex.diff.r <- this.lod.sex.int.r - this.lod.add.r
              
              this.lod.add.p <- lod.protein.addcovariatesonly[rownames(lod.protein.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.p <- lod.protein.dietint[rownames(lod.protein.dietint)==qtl.marker]
              this.lod.sex.int.p <- lod.protein.sexint[rownames(lod.protein.sexint)==qtl.marker]
              this.lod.diet.diff.p <- this.lod.diet.int.p - this.lod.add.p
              this.lod.sex.diff.p <- this.lod.sex.int.p - this.lod.add.p
              
              print(paste("Adding marker to table: ", qtl.marker))
              
              if (peaks.only==F){
                results.rna.diet.final[nrow(results.rna.diet.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r,
                                                this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
              }
              
              if (peaks.only==T){
                for (peak in qtl2.peakfinder$pos){
                  if (peak == qtl.marker.data$Mbp*1e6){
                    results.rna.diet.final[nrow(results.rna.diet.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                                                 rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                                                 qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                 this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r,
                                                                                 this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
                  }
                }
              }
            }
          }
        }
      
      if (rna.sex.int==T){
        print("eQTL/Sex Table Processing")
        qtl2.peakfinder <- find_peaks(lod.rna.sexint - lod.rna.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)
          
          for (qtl.marker in rownames(lod.rna.sexint)){
            if ((lod.rna.sexint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){

              qtl.marker.data <- markers[markers$Marker==qtl.marker,]
              
              # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
              if (window.type=="center"){
                gene.center <- abs(rna.annotations$Mbp.Start[index.rna] - rna.annotations$Mbp.End[index.rna]) / 2
                qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
              }
              else if (window.type=="edge"){
                gene.left.edge <- rna.annotations$Mbp.Start[index.rna]
                gene.right.edge <- rna.annotations$Mbp.End[index.rna]
                if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                  closer.edge <- gene.left.edge
                }
                else {
                  closer.edge <- gene.right.edge
                }
                qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
              }
              
              # Qualify QTL as either local or distant depending on the selected window size.
              if (qtl.distance > window.size.mbp/2){
                qtl.type <- "distant"
              }
              else {
                qtl.type <- "local"
              }
              
              this.lod.add.r <- lod.rna.addcovariatesonly[rownames(lod.rna.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.r <- lod.rna.dietint[rownames(lod.rna.dietint)==qtl.marker]
              this.lod.sex.int.r <- lod.rna.sexint[rownames(lod.rna.sexint)==qtl.marker]
              this.lod.diet.diff.r <- this.lod.diet.int.r - this.lod.add.r
              this.lod.sex.diff.r <- this.lod.sex.int.r - this.lod.add.r
              
              this.lod.add.p <- lod.protein.addcovariatesonly[rownames(lod.protein.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.p <- lod.protein.dietint[rownames(lod.protein.dietint)==qtl.marker]
              this.lod.sex.int.p <- lod.protein.sexint[rownames(lod.protein.sexint)==qtl.marker]
              this.lod.diet.diff.p <- this.lod.diet.int.p - this.lod.add.p
              this.lod.sex.diff.p <- this.lod.sex.int.p - this.lod.add.p
              
              print(paste("Adding marker to table: ", qtl.marker))
              
              if (peaks.only==F){
                results.rna.sex.final[nrow(results.rna.sex.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                               rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                               qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                               this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r,
                                               this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
              }
              
              if (peaks.only==T){
                for (peak in qtl2.peakfinder$pos){
                  if (peak == qtl.marker.data$Mbp*1e6){
                    results.rna.sex.final[nrow(results.rna.sex.final)+1,] <- c(rna.annotations$ENSMUSG[index.rna], rna.annotations$Gene[index.rna], rna.annotations$Chr[index.rna],
                                                                               rna.annotations$Mbp.Start[index.rna], rna.annotations$Mbp.End[index.rna], num.expressed.r,
                                                                               qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                               this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r,
                                                                               this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p)
                  }
                }
              }
            } 
          }
        }
      
      if (protein.diet.int==T){
        print("pQTL/Diet Table Processing")
        qtl2.peakfinder <- find_peaks(lod.protein.dietint - lod.protein.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)
          
          for (qtl.marker in rownames(lod.protein.dietint)){
            if ((lod.protein.dietint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){

              qtl.marker.data <- markers[markers$Marker==qtl.marker,]
              
              # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
              if (window.type=="center"){
                gene.center <- abs(protein.annotations$Mbp.Start[index.protein] - protein.annotations$Mbp.End[index.protein]) / 2
                qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
              }
              else if (window.type=="edge"){
                gene.left.edge <- protein.annotations$Mbp.Start[index.protein]
                gene.right.edge <- protein.annotations$Mbp.End[index.protein]
                if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                  closer.edge <- gene.left.edge
                }
                else {
                  closer.edge <- gene.right.edge
                }
                qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
              }
              
              # Qualify QTL as either local or distant depending on the selected window size.
              if (qtl.distance > window.size.mbp/2){
                qtl.type <- "distant"
              }
              else {
                qtl.type <- "local"
              }
              
              this.lod.add.r <- lod.rna.addcovariatesonly[rownames(lod.rna.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.r <- lod.rna.dietint[rownames(lod.rna.dietint)==qtl.marker]
              this.lod.sex.int.r <- lod.rna.sexint[rownames(lod.rna.sexint)==qtl.marker]
              this.lod.diet.diff.r <- this.lod.diet.int.r - this.lod.add.r
              this.lod.sex.diff.r <- this.lod.sex.int.r - this.lod.add.r
              
              this.lod.add.p <- lod.protein.addcovariatesonly[rownames(lod.protein.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.p <- lod.protein.dietint[rownames(lod.protein.dietint)==qtl.marker]
              this.lod.sex.int.p <- lod.protein.sexint[rownames(lod.protein.sexint)==qtl.marker]
              this.lod.diet.diff.p <- this.lod.diet.int.p - this.lod.add.p
              this.lod.sex.diff.p <- this.lod.sex.int.p - this.lod.add.p
              
              print(paste("Adding marker to table: ", qtl.marker))
              
              if (peaks.only==F){
                results.protein.diet.final[nrow(results.protein.diet.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                    protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                    qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                    this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p,
                                                    this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r)
              }
              
              if (peaks.only==T){
                for (peak in qtl2.peakfinder$pos){
                  if (peak == qtl.marker.data$Mbp*1e6){
                    results.protein.diet.final[nrow(results.protein.diet.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                                                         protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                                                         qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                         this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p,
                                                                                         this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r)
                  }
                }
              }
            }
          }
        }
      
      if (protein.sex.int==T){
        print("pQTL/Sex Table Processing")
        qtl2.peakfinder <- find_peaks(lod.protein.sexint - lod.protein.addcovariatesonly, map=map.peakfinder, threshold = threshold, drop = 1.5, peakdrop = 1.8)
          
          for (qtl.marker in rownames(lod.protein.sexint)){
            if ((lod.protein.sexint[qtl.marker,] - lod.rna.addcovariatesonly[qtl.marker,]) >= threshold){

              qtl.marker.data <- markers[markers$Marker==qtl.marker,]
              
              # Determine the distance from the gene to the QTL, either by the center of the gene or the closer edge of the gene.
              if (window.type=="center"){
                gene.center <- abs(protein.annotations$Mbp.Start[index.protein] - protein.annotations$Mbp.End[index.protein]) / 2
                qtl.distance <- abs(qtl.marker.data$Mbp - gene.center)
              }
              else if (window.type=="edge"){
                gene.left.edge <- protein.annotations$Mbp.Start[index.protein]
                gene.right.edge <- protein.annotations$Mbp.End[index.protein]
                if (abs(qtl.marker.data$Mbp - gene.left.edge) < abs(qtl.marker.data$Mbp - gene.right.edge)) {
                  closer.edge <- gene.left.edge
                }
                else {
                  closer.edge <- gene.right.edge
                }
                qtl.distance <- abs(qtl.marker.data$Mbp - closer.edge)
              }
              
              # Qualify QTL as either local or distant depending on the selected window size.
              if (qtl.distance > window.size.mbp/2){
                qtl.type <- "distant"
              }
              else {
                qtl.type <- "local"
              }
              
              this.lod.add.r <- lod.rna.addcovariatesonly[rownames(lod.rna.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.r <- lod.rna.dietint[rownames(lod.rna.dietint)==qtl.marker]
              this.lod.sex.int.r <- lod.rna.sexint[rownames(lod.rna.sexint)==qtl.marker]
              this.lod.diet.diff.r <- this.lod.diet.int.r - this.lod.add.r
              this.lod.sex.diff.r <- this.lod.sex.int.r - this.lod.add.r
              
              this.lod.add.p <- lod.protein.addcovariatesonly[rownames(lod.protein.addcovariatesonly)==qtl.marker]
              this.lod.diet.int.p <- lod.protein.dietint[rownames(lod.protein.dietint)==qtl.marker]
              this.lod.sex.int.p <- lod.protein.sexint[rownames(lod.protein.sexint)==qtl.marker]
              this.lod.diet.diff.p <- this.lod.diet.int.p - this.lod.add.p
              this.lod.sex.diff.p <- this.lod.sex.int.p - this.lod.add.p
              
              print(paste("Adding marker to table: ", qtl.marker))
              
              if (peaks.only==F){
                results.protein.sex.final[nrow(results.protein.sex.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                   protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                   qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                   this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p,
                                                   this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r) 
              }
              
              if (peaks.only==T){
                for (peak in qtl2.peakfinder$pos){
                  if (peak == qtl.marker.data$Mbp*1e6){
                    results.protein.sex.final[nrow(results.protein.sex.final)+1,] <- c(protein.annotations$ENSMUSP[index.protein], protein.annotations$ENSMUSG[index.protein], protein.annotations$Gene[index.protein], protein.annotations$Chr[index.protein],
                                                                                       protein.annotations$Mbp.Start[index.protein], protein.annotations$Mbp.End[index.protein], num.expressed.p,
                                                                                       qtl.marker.data$Marker, qtl.marker.data$Chr, qtl.marker.data$Mbp, qtl.marker.data$cM, qtl.type, qtl.distance,
                                                                                       this.lod.add.p, this.lod.diet.int.p, this.lod.sex.int.p, this.lod.diet.diff.p, this.lod.sex.diff.p,
                                                                                       this.lod.add.r, this.lod.diet.int.r, this.lod.sex.int.r, this.lod.diet.diff.r, this.lod.sex.diff.r) 
                  }
                }
              }
            }
          }
        }
      }
    }
  if (rna == T){
    if (rna.diet.int==T){
      setwd(dir.output)
      if (dir.exists("rna_collate")){
        setwd("rna_collate")
      }
      else{
        dir.create("rna_collate")
        setwd("rna_collate")
      }
      write.table(results.rna.diet.final, file.path(file.path(dir.output, "rna_collate"), paste(dataset.name, "_LOD", threshold, "_Window", window.size.mbp, "Mbp_", window.type, "_DietInt_eQTL.tsv", sep = "")), sep="\t")
    }
    if (rna.sex.int==T){
      setwd(dir.output)
      if (dir.exists("rna_collate")){
        setwd("rna_collate")
      }
      else{
        dir.create("rna_collate")
        setwd("rna_collate")
      }
      write.table(results.rna.sex.final, file.path(file.path(dir.output, "rna_collate"), paste(dataset.name, "_LOD", threshold, "_Window", window.size.mbp, "Mbp_", window.type, "_SexInt_eQTL.tsv", sep = "")), sep="\t")
    }
  }
  if (protein==T){
    if (protein.diet.int==T){
      setwd(dir.output)
      if (dir.exists("protein_collate")){
        setwd("protein_collate")
      }
      else{
        dir.create("protein_collate")
        setwd("protein_collate")
      }
      write.table(results.protein.diet.final, file.path(file.path(dir.output, "protein_collate"), paste(dataset.name, "_LOD", threshold, "_Window", window.size.mbp, "Mbp_", window.type, "_DietInt_pQTL.tsv", sep = "")), sep="\t")
    }
    if (protein.sex.int==T){
      setwd(dir.output)
      if (dir.exists("protein_collate")){
        setwd("protein_collate")
      }
      else{
        dir.create("protein_collate")
        setwd("protein_collate")
      }
      write.table(results.protein.sex.final, file.path(file.path(dir.output, "protein_collate"), paste(dataset.name, "_LOD", threshold, "_Window", window.size.mbp, "Mbp_", window.type, "_SexInt_pQTL.tsv", sep = "")), sep="\t")
    }
  }
}
}

################################################################################################################################################

manhattan.generic <- function(){
  
}

################################################################################################################################################

ome.map.generic <- function(){
  
}

################################################################################################################################################

moving.window.generic <- function(){
  
}

################################################################################################################################################

coefficient.generic <- function(){
  
}

################################################################################################################################################

mediation.generic <- function(){
  
}


