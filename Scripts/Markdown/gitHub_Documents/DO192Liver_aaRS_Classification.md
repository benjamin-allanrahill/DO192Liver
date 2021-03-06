Local/Distant Classification
================
Ben Allan-Rahill
6/28/2018

-   [Local/Distant Classification](#localdistant-classification)
    -   [Libraries](#libraries)
    -   [Classification](#classification)
    -   [CLASSIFICATION FOR DIFFERENCE PEAKS](#classification-for-difference-peaks)

Local/Distant Classification
============================

### Libraries

``` r
library(tidyverse)
```

Classification
--------------

This is for peaks that are not thresholded by interactive difference.

``` r
peak_data_path <- "/Users/c-allanb/Desktop/Benjamin_Allan-Rahill/DO192Liver/Data/Output/peak_data"
gene.names <- read.csv("/Users/c-allanb/Desktop/DO192Liver/Data/Input/aa-tRNA_Sequence_Names.csv", header = FALSE)

# Proteins

# Set working directory to the 'proteins' folder
setwd(paste(peak_data_path, "/Proteins", sep=""))
getwd()

# list files from current working directory 
file_names.prot <- list.files()

# Establish a global counter variable (g) to be used later to iterate through gene names 
g <-  1

for (file.p in file_names.prot){
  
  # Use the list of Gene names to name the files
  gene_names <-  gene.names
  gene_file_name <-  paste(gene_names[[1]][g]) # Use the counter to locate the gene in the list
  g <- g + 1 # Add to the counter so the next for loop uses the next gene 
  
  gene_name <- gene_file_name 
  
  # Read in the .tsv file
  peak_data.prot <- read.table(file.p, header = TRUE, sep = "\t")
  
  # Calculate center of gene 
  gene_center <- mean(c(peak_data.prot$Gene_End, peak_data.prot$Gene_Start))
  qtl_distance <- abs(peak_data.prot$pos - gene_center)
  
  # Add new column of classification 
    
  # Determine the QTL Type based off distance and chromisome
  qtl_type <- ifelse(peak_data.prot$chr != peak_data.prot$Chromosome | qtl_distance >= 2.5, "Distant", "Local")

  # Add to df with QTL type and distance 
  peak_data.prot <- mutate(peak_data.prot, QTL_Distance = qtl_distance, QTL_Type = qtl_type)

  # save file
  write.table(peak_data.prot, file = paste(gene_name, "_peaks_comb_Prot.tsv", sep =""), sep = "\t", row.names = FALSE)
  
}

# RNA

# Set working directory to the 'proteins' folder
setwd(paste(peak_data_path, "/RNA", sep=""))
getwd()

# list files from current working directory 
file_names.rna <- list.files()

# Establish a global counter variable (g) to be used later to iterate through gene names 
g <-  1

for (file.rna in file_names.rna){
  
  # Use the list of Gene names to name the files
  gene_names <-  gene.names
  gene_file_name <-  paste(gene_names[[1]][g]) # Use the counter to locate the gene in the list
  g <- g + 1 # Add to the counter so the next for loop uses the next gene 
  
  gene_name <- gene_file_name 
  
  # Read in the .tsv file
  peak_data.rna <- read.table(file.rna, header = TRUE, sep = "\t")
  
  # Calculate center of gene 
  gene_center <- mean(c(peak_data.rna$Gene_End, peak_data.rna$Gene_Start))
  qtl_distance <- abs(peak_data.rna$pos - gene_center)
  
 
   # Add new column of classification #
    
  # Determine the QTL Type based off distance and chromisome
  qtl_type <- ifelse(peak_data.rna$chr != peak_data.rna$Chromosome | qtl_distance >= 2.5, "Distant", "Local")

  # Add to df with QTL type and distance 
  peak_data.rna <- mutate(peak_data.rna, QTL_Distance = qtl_distance, QTL_Type = qtl_type)
  
  
  write.table(peak_data.rna, file = paste(gene_name, "_peaks_comb_RNA.tsv", sep =""), sep = "\t", row.names = FALSE)
}
```

CLASSIFICATION FOR DIFFERENCE PEAKS
-----------------------------------

``` r
# Read in large tsv of all the peaks 
diff_tsv.p <- as_tibble(read.table("/Users/c-allanb/Desktop/DO192Liver/Data/Output/Scans/pQTL_coef_data", sep = "\t", header = TRUE)) %>%
  filter(lod >= 7)
```

    ## Warning: package 'bindrcpp' was built under R version 3.4.4

``` r
diff_tsv.r <- as_tibble(read.table("/Users/c-allanb/Desktop/DO192Liver/Data/Output/Scans/eQTL_coef_data", sep = "\t", header = TRUE)) %>%
  filter(lod >= 7)
```

``` r
# Calculate distance and type 

qtl_dist_classification <- function(gene_end, gene_start, qtl_pos, qtl_chr, gene_chr){
  
  gene_center <- mean(c(gene_end, gene_start))
  qtl_distance <- abs(qtl_pos - gene_center)

  # Determine the QTL Type based off distance and chromisome
  qtl_type <- ifelse(qtl_chr != gene_chr | qtl_distance >= 2.5, "Distant", "Local")
  
  return(qtl_type)
}
```

``` r
# Rowise read in the arguments to the function to classify the QTL

diff_tsv.p <- diff_tsv.p %>%
  rowwise()%>%
  mutate(QTL_Type = qtl_dist_classification(Gene_End, Gene_Start, pos, chr, Chromosome))

diff_tsv.r <- diff_tsv.r %>%
  rowwise()%>%
  mutate(QTL_Type = qtl_dist_classification(Gene_End, Gene_Start, pos, chr, Chromosome))
```

``` r
write.table(diff_tsv.p, file = "/Users/c-allanb/Desktop/DO192Liver/Data/Output/peak_data/pQTL_coef_data", sep = "\t")
write.table(diff_tsv.r, file = "/Users/c-allanb/Desktop/DO192Liver/Data/Output/peak_data/eQTL_coef_data", sep = "\t")
```
