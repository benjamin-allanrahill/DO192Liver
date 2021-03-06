DO192Liver aaRS Peaks and Founder Strain Plotting
================
Ben Allan-Rahill
6/27/2018

-   [Purpose](#purpose)
    -   [Import Libraries and load data & Gene Names](#import-libraries-and-load-data-gene-names)
-   [Multiple Genes](#multiple-genes)
-   [For Loop](#for-loop)
    -   [Peaks Diff Function](#peaks-diff-function)
    -   [Covar Assign Function](#covar-assign-function)
    -   [Marker Position Function](#marker-position-function)
    -   [For Loop for finding the difference of peaks for each gene](#for-loop-for-finding-the-difference-of-peaks-for-each-gene)
-   [COEF SCANS and PlOTS](#coef-scans-and-plots)
    -   [Define functions](#define-functions)
    -   [COEF \#\#s](#coef-s)
    -   [Save the data frame](#save-the-data-frame)

Purpose
-------

In this script we will find the peaks of each QTL scan so we can locate the marker location of that peak. We will use this loacation to find the founder allele effects at that location. This will allow us to relate the QTLs to the founders, allowing us to understand the tRNA charging pathways of the 8 founder strains.

### Import Libraries and load data & Gene Names

``` r
library(qtl2)
library(dplyr)
library(MALDIquant)

# load for annotations 
load("/Users/c-allanb/Desktop/DO192Liver/Data/Input/DO192LiverData_Formattedfor_rQTL2.Rdata")

#load to be used in for loop (naming files)
gene.names <- read.csv("/Users/c-allanb/Desktop/DO192Liver/Data/Input/aa-tRNA_Sequence_Names.csv", header = FALSE)
```

Multiple Genes
--------------

This part of the script is to be used with many scan files. I used 46 files and the whole loop took less than a minute. If more files are used, it can be run on the HPC.

Here we will establish the files to iterate over and the working directory to send the finished files to.

``` r
# Establish a set path to use for this code chunk 
scans_path <- '/Users/c-allanb/Desktop/DO192Liver/Data/Input/'

# Generate list of scan file names (list of 44)
scan_file_names <- list.files(paste(scans_path,'Scans', sep = ""))

# Set the working directory to '/peak_data' (only need to run once)
dir.create(paste(scans_path, "peak_data", sep= ""))
```

    ## Warning in dir.create(paste(scans_path, "peak_data", sep = "")): '/Users/c-
    ## allanb/Desktop/DO192Liver/Data/Input/peak_data' already exists

``` r
setwd(paste(scans_path,'peak_data', sep = ""))

# Create two directories in the '/Plots' folder. One for RNA and one for proteins. (ONLY RUN ONCE)
dir.create(paste(scans_path, "peak_data/", "Proteins",  sep= ""))
```

    ## Warning in dir.create(paste(scans_path, "peak_data/", "Proteins", sep =
    ## "")): '/Users/c-allanb/Desktop/DO192Liver/Data/Input/peak_data/Proteins'
    ## already exists

``` r
dir.create(paste(scans_path, "peak_data/", "RNA",  sep= ""))
```

    ## Warning in dir.create(paste(scans_path, "peak_data/", "RNA", sep = "")): '/
    ## Users/c-allanb/Desktop/DO192Liver/Data/Input/peak_data/RNA' already exists

For Loop
--------

This for loop will take all the scan files and collate the correct data (i.e. gene name, LOD score, position). It will then save the collated data in .tsv files.

FOR CLARIFICATION:

The `marker intervals` are just collated tables combining the LOD scores of a scan (additive/interactive) to the marker and base pair positions using the `Markers.64k`. This is for the purpose of relating LOD scores to their marker position to be used in the later founder coefficient part.

### Peaks Diff Function

This function is called later in the for loop and is takes in the marker table and wether it is diet or sex interactive LOD scores and it takes the additive LOD peak and it's confidence intervals (lo & hi)

``` r
# Establish the function that will find the difference between the additive and interactive peaks
peaks_diff <- function(marker_interval_type, marker_interval, lod_peak, lo, hi) {
  
  #Test
  print(paste(marker_interval_type))
  
  # Calculate the difference for Diet
  if (marker_interval_type == "Diet"){
    
    # find the max LOD score in the base pair interval (found from the confidence interval which was retrieved from find peaks)
    max_lod_diet <- max(marker_interval[which(marker_interval$bp==lo):which(marker_interval$bp==hi),]$pheno1)
    
    print(max_lod_diet)
    
    #Make generic
    max_lod <- max_lod_diet

    # Subtract the max LOD from that interval from the additive peak
    subtr_peak <- lod_peak - max_lod

    return(subtr_peak)
  
  }
  
  #Calculate the difference for Sex
  if (marker_interval_type == "Sex") {
    
    # find the max LOD score in the base pair interval (found from the confidence interval which was retrieved from find peaks)
    max_lod_sex <- max(marker_interval[which(marker_interval$bp==lo):which(marker_interval$bp==hi),]$pheno1)
   
    max_lod <- max_lod_sex
    
    # Subtract the max LOD from that interval from the additive peak
    subtr_peak <- lod_peak - max_lod

    return(subtr_peak)

  }
  
}
```

### Covar Assign Function

This function takes in the diet and sex difference (calculated in the `peaks diff` function). It checks to see if they are greater than 4. If so, it assigns the covaraite value as the either "Diet Int" or "Sex Int".

``` r
# Very basic function to return the type of QTL it is
covar_assign <-  function(diet, sex){
  if ( abs(diet) >= 4) { 
    return("Diet Int")
  }

  if ( abs(sex) >= 4){
    
    return("Sex Int")
  }

  else {
    return("Additive")
  }
}
```

### Marker Position Function

This function takes: \* The additive QTL position \* The additive QTL covariate \* The peak confidence interval beginning and end \* All 6 marker data frames \* Wether it is and eQTL or a pQTL

The function takes these arguments and find the marker position closest to the QTL. If the QTL is interactive, it calculates the base pair postion of the max interactive LOD score in the confidence interval and then gets the marker position and that locus. It then ret

``` r
marker_pos <- function(position, chromosome, covariate, int_beg, int_end, marker_int.a.p, marker_int.a.r, marker_int.s.p, marker_int.d.p, marker_int.s.r, marker_int.d.r, kind){

  if (kind == "pQTL"){
    if (covariate == "Additive"){
      
      marker_int.a.p <- marker_int.a.p%>%
        filter(chr == chromosome)
      
      marker_pos.add <- marker_int.a.p$marker[match.closest(position, marker_int.a.p$bp)]
     
      return(marker_pos.add)
    }
    else{
      
      if (covariate == "Diet Int"){
        
        # subset the marker interval based off the chromosome and the index positions from above
        marker_int.d.p <- marker_int.d.p %>%
          filter(chr == chromosome)
          
        # establish the index from the Markers using the interval beginning and end 
        diet_int_beg <- which(marker_int.d.p$bp==int_beg)
        diet_int_end <- which(marker_int.d.p$bp==int_end)
        
        #subset the marker df into just the rows within the interval
        marker_int.d.p <- slice(marker_int.d.p, diet_int_beg:diet_int_end)
       
        # retrieve the base pair position where the lod score is the highest
        bp_pos.diet <- marker_int.d.p$bp[which.max(marker_int.d.p$pheno1)]

        # use `match.closest` to get the closest marker to the max lod score
        marker_pos.diet <- marker_int.d.p$marker[match.closest(bp_pos.diet, marker_int.d.p$bp)]
    
        return(marker_pos.diet)
    
      }
    
      if (covariate == "Sex Int") {
        
        # subset the marker interval based off the chromosome and the index positions from above
        marker_int.s.p <- marker_int.s.p %>%
          filter(chr == chromosome)
          
        # establish the index from the Markers using the interval beginning and end 
        sex_int_beg <- which(marker_int.s.p$bp==int_beg)
        sex_int_end <- which(marker_int.s.p$bp==int_end)
        
        #subset the marker df into just the rows within the interval
        marker_int.s.p <- slice(marker_int.s.p, sex_int_beg:sex_int_end)

        # retrieve the base pair position where the lod score is the highest
        bp_pos.sex <- marker_int.s.p$bp[which.max(marker_int.s.p$pheno1)]

        # use `match.closest` to get the closest marker to the max lod score
        marker_pos.sex <- marker_int.s.p$marker[match.closest(bp_pos.sex, marker_int.s.p$bp)]
        print(paste("Sex MARKER.p", marker_pos.sex))
        
        # return the marker to the mutate function 
        return(marker_pos.sex)
      }
    }
  }
  

  ### RNA ###
  
  if (kind == "eQTL"){
    if (covariate == "Additive"){
      marker_int.a.r <- marker_int.a.r %>%
        filter(chr == chromosome)

      marker_pos.add <- marker_int.a.r$marker[match.closest(position, marker_int.a.r$bp)]
      
      return(marker_pos.add)
    }
    else{
      
      if (covariate == "Diet Int"){
        # subset the marker interval based off the chromosome and the index positions from above
        marker_int.d.r <- marker_int.d.r %>%
          filter(chr == chromosome)
          
        # establish the index from the Markers using the interval beginning and end 
        diet_int_beg <- which(marker_int.d.r$bp==int_beg)
        diet_int_end <- which(marker_int.d.r$bp==int_end)
        
        #subset the marker df into just the rows within the interval
        marker_int.d.r <- slice(marker_int.d.r, diet_int_beg:diet_int_end)
       
        # retrieve the base pair position where the lod score is the highest
        bp_pos.diet <- marker_int.d.r$bp[which.max(marker_int.d.r$pheno1)]

        # use `match.closest` to get the closest marker to the max lod score
        marker_pos.diet <- marker_int.d.r$marker[match.closest(bp_pos.diet, marker_int.d.r$bp)]
    
        return(marker_pos.diet)
    
      }
    
      if (covariate == "Sex Int") {
        
        # subset the marker interval based off the chromosome and the index positions from above
        marker_int.s.r <- marker_int.s.r %>%
          filter(chr == chromosome)
          
        # establish the index from the Markers using the interval beginning and end 
        sex_int_beg <- which(marker_int.s.r$bp==int_beg)
        sex_int_end <- which(marker_int.s.r$bp==int_end)
        
        #subset the marker df into just the rows within the interval
        marker_int.s.r <- slice(marker_int.s.r, sex_int_beg:sex_int_end)

        # retrieve the base pair position where the lod score is the highest
        bp_pos.sex <- marker_int.s.r$bp[which.max(marker_int.s.r$pheno1)]

        # use `match.closest` to get the closest marker to the max lod score
        marker_pos.sex <- marker_int.s.r$marker[match.closest(bp_pos.sex, marker_int.s.r$bp)]
        
        return(marker_pos.sex)
      }
    }
  }
  

}
```

### For Loop for finding the difference of peaks for each gene

``` r
# establish important global variables for the `for loop`

#establish empty data frame to add to inside of the loop 
peaks_diff_comb_df <- data.frame()

# set g = 1
g <- 1
```

``` r
# Iterate through files in the list of file names 
for (f in scan_file_names) {
  
  if (g <= 44) {
    
    # Use the list of Gene names to name the files
    gene_names <-  gene.names
    gene_file_name <-  paste(gene_names[[1]][g]) # Use the counter to locate the gene in the list
    g <- g + 1 # Add to the counter so the next for loop uses the next gene 
    
    gene_name <- gene_file_name 
    
    # Load the file to be able to access the scan data 
    load(paste(scans_path, "/Scans/", f, sep = ""))
    
        
    # join LOD scores with basepair pos and markers
    # TODO: Change name from "marker interval" to something more explanatory.
    marker_interval_add.p <- bind_cols(data.frame(lod.protein.addcovariatesonly),markers.64K) 
    
    marker_interval_diet.p <- bind_cols(data.frame(lod.protein.dietint),markers.64K) 
    
    marker_interval_sex.p <- bind_cols(data.frame(lod.protein.sexint),markers.64K) 
    
    marker_interval_add.r <- bind_cols(data.frame(lod.rna.addcovariatesonly),markers.64K) 
    
    marker_interval_diet.r <- bind_cols(data.frame(lod.rna.dietint),markers.64K) 
    
    marker_interval_sex.r <- bind_cols(data.frame(lod.rna.sexint),markers.64K) 
                                       ## FIND PEAKS ##
    
    #########################
    ### PROTEIN
    #########################
    
    # Run these find peaks function for each LOD scan 
    peaks_additive.p <- find_peaks(lod.protein.addcovariatesonly, 
                                   map = map, 
                                   threshold = 5, 
                                   drop = 1.5, 
                                   peakdrop = 1.8)
    
    # run the `peaks diff` function in the mutate call using rowwise so it runs for each peak
    peaks_diff.p <- peaks_additive.p %>%
      rowwise() %>%
      mutate(diff_diet = peaks_diff("Diet", marker_interval_diet.p, lod, ci_lo, ci_hi), 
         diff_sex = peaks_diff("Sex", marker_interval_sex.p, lod, ci_lo, ci_hi)) 
     
    # assign the covariates based off their diet and sex difference  
    peaks_diff.p <- peaks_diff.p %>%
      rowwise() %>%
      mutate(Covariates = covar_assign(diff_diet, diff_sex))
    
    # if there is a diet or sex significant peak, duplicate the row w/ the covariate at "Additive"
    dup_int_peaks_diff.p <- peaks_diff.p %>%
      filter(Covariates == "Diet Int" | Covariates == "Sex Int") %>%
      mutate(Covariates = rep("Additive"))
    
    # Bind both dfs back together
    peaks_diff.p <- bind_rows(peaks_diff.p, dup_int_peaks_diff.p) %>%
      arrange(chr)

    # add "pQTL" into 'Type' column; aslo add the gene name into a column
    peaks_diff.p <- peaks_diff.p %>%
      mutate(Type = rep("pQTL"), Gene = rep(gene_name)) %>%
      select(Type, Gene, chr, pos, ci_lo, ci_hi, lod, diff_diet, diff_sex, Covariates)
    
    ##PROTEIN DATA##
    
    # retrieve gene specific data from the annotations data frame. This is important for classifying the QTL as distant or local 
    gene_specific_data.p <- filter(annotations.protein.192, annotations.protein.192$Associated.Gene.Name == gene_name) %>%
      select(
        Ensembl_Gene_ID = Ensembl.Gene.ID, 
        Ensembl_Protein_ID = Ensembl.Protein.ID, 
        Gene_Name = Associated.Gene.Name, 
        Chromosome = Chromosome.Name, 
        Gene_Start = Gene.Start..bp., 
        Gene_End = Gene.End..bp.) 
    
    # Add rows so that the combine colls will work 
      if (nrow(peaks_diff.p) > 1){
        gene_specific_data.p <- add_row(gene_specific_data.p, 
          Ensembl_Gene_ID = rep(gene_specific_data.p[1,1], times = nrow(peaks_diff.p)-1), 
          Ensembl_Protein_ID = gene_specific_data.p[1,2], 
          Gene_Name = gene_specific_data.p[1,3], 
          Chromosome = gene_specific_data.p[1,4], 
          Gene_Start = gene_specific_data.p[1,5], 
          Gene_End = gene_specific_data.p[1,6])
      }
    
    # Combine QTLs and their gene specific data 
    peaks_diff_comp.p <- bind_cols(gene_specific_data.p, peaks_diff.p)
    
    
     #########################
    ### RNA
    #########################
    
    ### Repeat process above for RNA
    
    # Run these find peaks function for each LOD scan 
    peaks_additive.rna <- find_peaks(lod.rna.addcovariatesonly, 
                                     map = map, 
                                     threshold = 5, 
                                     drop = 1.5, 
                                     peakdrop = 1.8)
    
    peaks_diff.r <- peaks_additive.rna %>%
      rowwise() %>%
      mutate(diff_diet = peaks_diff("Diet", marker_interval_diet.r, lod, ci_lo, ci_hi), 
         diff_sex = peaks_diff("Sex", marker_interval_sex.r, lod, ci_lo, ci_hi)) 
    
    
    peaks_diff.r <- peaks_diff.r %>%
      rowwise() %>%
      mutate(Covariates = covar_assign(diff_diet, diff_sex))
    
    dup_int_peaks_diff.r <- peaks_diff.r %>%
      filter(Covariates == "Diet Int" | Covariates == "Sex Int") %>%
      mutate(Covariates = rep("Additive"))
    
    peaks_diff.r <- bind_rows(peaks_diff.r, dup_int_peaks_diff.r) %>%
      arrange(chr)    
    
    peaks_diff.r <- peaks_diff.r %>%
      mutate(Type = rep("eQTL"), Gene = rep(gene_name)) %>%
      select(Type, Gene, chr, pos, ci_lo, ci_hi, lod, diff_diet, diff_sex, Covariates)
    
    ##RNA## 
    gene_specific_data.rna <-  filter(annotations.rna.192, annotations.rna.192$Gene == gene_name) %>%
      select(Ensembl_Gene_ID = EnsemblID, 
             Gene_Name = Gene, 
             Chromosome = Chr, 
             Gene_Start = Start.Mbp, 
             Gene_End = End.Mbp)
    
    # Add rows so that the combine colls will work 
    if (nrow(peaks_diff.r) > 1){
      gene_specific_data.rna <- add_row(gene_specific_data.rna, 
        Ensembl_Gene_ID = rep(gene_specific_data.rna[1,1], times = nrow(peaks_diff.r) -1), 
        Gene_Name = gene_specific_data.rna[1,2], 
        Chromosome = gene_specific_data.rna[1,3], 
        Gene_Start = gene_specific_data.rna[1,4], 
        Gene_End = gene_specific_data.rna[1,5])      
    }

      
    peaks_diff_comp.r <- bind_cols(gene_specific_data.rna, peaks_diff.r)
    
    
    # Add to data frame with the combined pQTL and eQTL data 
    peaks_diff_comb_df <- bind_rows(peaks_diff_comb_df, peaks_diff_comp.p, peaks_diff_comp.r)
    
    
    # Determine the marker positon of each QTL using the `marker pos` function
    peaks_diff_marker_df <- peaks_diff_comb_df %>%
      rowwise() %>%
      mutate(QTL_Marker_pos = marker_pos(pos, as.numeric(chr), Covariates, ci_lo, ci_hi, marker_interval_add.p,  marker_interval_add.r, marker_interval_sex.p, marker_interval_diet.p, marker_interval_sex.r, marker_interval_diet.r, Type))

    
    }
  } 
```

``` r
# save the table 
write.table(peaks_diff_marker_df, file = "/Users/c-allanb/Desktop/DO192Liver/Data/Output/peak_data/diff_data", row.names = FALSE, sep = "\t")
```

COEF SCANS and PlOTS
--------------------

This part of the script takes in the distinct peaks (one per gene per chroomosome) and determines the founder coefficients at that locus. These are then returned to the same dataframe

``` r
file_path <- "/Users/c-allanb/Desktop/DO192Liver/Data/Output/peak_data/diff_data"

#read in saved .tsv
peaks_diff_data <- read.table(file_path, header = TRUE, sep = "\t")

# collate table to get rid of dupes
distinct_peaks_data <-  distinct(peaks_diff_data, chr, Gene, Type, Covariates, .keep_all = TRUE)
```

``` r
# reestablish  paths 
scans_path <- '/Users/c-allanb/Desktop/DO192Liver/Data/Input/Scans/'

path <-  "/Users/c-allanb/Desktop/DO192Liver/"
```

### Define functions

#### Protein

``` r
coef_scan_plot.prot <-  function(chromosome, type, diff_diet, diff_sex, gene_name) {

  # Get the gene's index in the protein and rna data.
  index.protein <- which(annotations.protein.192$Associated.Gene.Name == gene_name)
  # print(index.protein)
  
  #collate the data to make sure there are no NA from the expression data
  mice.with.data.p <- !is.na(expr.protein.192[,index.protein])
  
  # Load Scan data to go with gene 
  load(paste(scans_path, list.files(path = scans_path, paste("^", gene_name, "_", sep= "")), sep = ""))

  #kinship 
  K.LOCO.192.chr.p <- K.LOCO.qtl2

  # Calculate the K.LOCO for each chromosome
  for (kloco.chr in 1:20){
    K.LOCO.192.chr.p[[kloco.chr]] <- K.LOCO.qtl2[[kloco.chr]][mice.with.data.p,mice.with.data.p]
  }   

  # Set up additive and interactive covariates for protein, only mice with data.
  addcovar.p <- model.matrix(~ Sex + Diet + Tag, data=covariates.protein.192[mice.with.data.p, ])
  intcovar.sex.p <- model.matrix(~ Sex, data=covariates.protein.192[mice.with.data.p, ])
  intcovar.diet.p <- model.matrix(~ Diet, data=covariates.protein.192[mice.with.data.p, ])

  if (chromosome != "20" && diff.d <=4 && diff.s <=4 ) {
    # Scan // Save 
    
    # Run the coef scan using the LOD scan loaded above
    coef_scan_add.prot <-  scan1coef(genoprobs = probs[mice.with.data.p, chromosome],
                                kinship=K.LOCO.192.chr.p[[chromosome]],
                                pheno=expr.protein.192[mice.with.data.p,index.protein],
                                addcovar=addcovar.p[,-1],
                                reml=TRUE )
    #  print(gene_name)
    
    # save the scan with a specific name 
    save(coef_scan_add.prot, file = paste(path, "Data/Output/Scans/pQTL/", gene_name, "_chr", chromosome, "_Add_coef_scan.Rdata", sep = ""))
    
    # Plot
    
    # Set the device to png to write a png file of the plot
    #   Use specified dimensions to make the plot easier to read
    png(filename= paste(path, "Data/Output/Plots/Founder_effects/", type, "/", gene_name, "_chr", chromosome, "_add_", type, sep =""), height = 1018, width = 1280 )
    
    # Create the plot with the CC founder strain colors and using the additive LOD scan (will appear on bottom of plot)
    coef_plot_add.p <- plot_coefCC(coef_scan_add.prot, map, scan1_output = lod.protein.addcovariatesonly)
    dev.off()
    #cat("Done with ", gene_name, "_chr", chromosome, "\n")
  }
  
  if (chromosome != "20" && diff_diet >=4) {
    # Scan // Save both additive and int scans for peak
    coef_scan_add.prot <-  scan1coef(genoprobs = probs[mice.with.data.p, chromosome],
                                kinship=K.LOCO.192.chr.p[[chromosome]],
                                pheno=expr.protein.192[mice.with.data.p,index.protein],
                                addcovar=addcovar.p[,-1],
                                reml=TRUE )
    # Save
    save(coef_scan_add.prot, file = paste(path, "Data/Output/Scans/pQTL/", gene_name, "_chr", chromosome, "_Add_coef_scan.Rdata", sep = ""))    
    
    # Run scan coef with interactive covariate matrix specified 
    coef_scan_diet_int.prot <- scan1coef(genoprobs=probs[mice.with.data.p, chromosome],
                                    kinship=K.LOCO.192.chr.p[[chromosome]],
                                    pheno=expr.protein.192[mice.with.data.p,index.protein],
                                    addcovar=addcovar.p[,-1],
                                    intcovar=intcovar.diet.p[,-1],
                                    reml=TRUE)
    # Save
    save(coef_scan_diet_int.prot, file = paste(path, "Data/Output/Scans/pQTL/", gene_name, "_chr", chromosome, "_Diet_coef_scan.Rdata", sep = ""))

    # Plot
    png(filename= paste(path, "Data/Output/Plots/Founder_effects/", type, "/", gene_name, "_chr", chromosome, "_diff_diet" , sep =""), height = 1018, width = 1280 )
    coef_plot_diet_int.p <- plot_coefCC(coef_scan_diet_int.prot, map, scan1_output = lod.protein.dietint)

    dev.off()
    #cat("Done with ", gene_name, "_chr", chromosome, "\n")
    
  }
  
  if (chromosome != "20" && diff_sex >=4){
    # Scan / Protein / Interactive covariate: Sex
    coef_scan_add.prot <-  scan1coef(genoprobs = probs[mice.with.data.p, chromosome],
                                kinship=K.LOCO.192.chr.p[[chromosome]],
                                pheno=expr.protein.192[mice.with.data.p,index.protein],
                                addcovar=addcovar.p[,-1],
                                reml=TRUE )
    save(coef_scan_add.prot, file = paste(path, "Data/Output/Scans/pQTL/", gene_name, "_chr", chromosome, "_Add_coef_scan.Rdata", sep = ""))    
    
    coef_scan_sex_int.prot <- scan1coef(genoprobs=probs[mice.with.data.p, chromosome],
                            kinship=K.LOCO.192.chr.p[[chromosome]],
                            pheno=expr.protein.192[mice.with.data.p,index.protein],
                            addcovar=addcovar.p[,-1],
                            intcovar=intcovar.sex.p[,-1],
                            reml=TRUE)
 
    save(coef_scan_sex_int.prot, file = paste(path, "Data/Output/Scans/pQTL/", gene_name, "_chr", chromosome, "_Sex_coef_scan.Rdata", sep = ""))
    
    #Plot
    par(mfrow = c(2,1))
    png(filename= paste(path, "Data/Output/Plots/Founder_effects/", type, "/", gene_name, "_chr", chromosome, "_diff_sex", sep =""), height = 1018, width = 1280 )
    coef_plot_sex_int.p <- plot_coefCC(coef_scan_sex_int.prot, map, scan1_output = lod.protein.sexint)
    dev.off()
    #cat("Done with ", gene_name, "_chr", chromosome, "\n")
    #
  }
  
}
```

#### RNA

``` r
# Same as function above with parameters set for eQTLS 
coef_scan_plot.rna <-  function(chromosome, type, diff_diet, diff_sex, gene_name) {
  #print(gene_name)
  
  # Get the gene's index
  index.rna <- which(annotations.rna.192$Gene == gene_name)
  #print(index.rna)
  mice.with.data.r <- !is.na(expr.rna.192[,index.rna])
  
  # kinship 
  K.LOCO.192.chr.r <- K.LOCO.qtl2
  

  for (kloco.chr in 1:20){
    K.LOCO.192.chr.r[[kloco.chr]]<- K.LOCO.qtl2[[kloco.chr]][mice.with.data.r,mice.with.data.r]
  }
 
  # Set up additive and interactive covariates for RNA, only mice with data.
  addcovar.r <- model.matrix(~ Sex + Diet, data=covariates.rna.192[mice.with.data.r, ])
  intcovar.sex.r <- model.matrix(~ Sex, data=covariates.rna.192[mice.with.data.r, ])
  intcovar.diet.r <- model.matrix(~ Diet, data=covariates.rna.192[mice.with.data.r, ])
    
  # Load scan data 
  load(paste(scans_path, list.files(path = scans_path, paste("^", gene_name, "_", sep= "")), sep = ""))
  
  if (chromosome != "20" && diff_diet <= 4 && diff_sex <=4) {
    coef_scan_add.rna <-  scan1coef(genoprobs = probs[mice.with.data.r, chromosome],
                                kinship=K.LOCO.192.chr.r[chromosome],
                                pheno=expr.rna.192[mice.with.data.r,index.rna],
                                addcovar=addcovar.r[,-1],
                                reml=TRUE )
    #print(gene_name)
    
    save(coef_scan_add.rna, file = paste(path, "Data/Output/Scans/eQTL/", gene_name, "_chr", chromosome, "_Add_coef_scan.Rdata", sep = ""))
    
    # Plot
    png(filename= paste(path, "Data/Output/Plots/Founder_effects/", type, "/", gene_name, "_chr", chromosome, "_add_", type, sep =""), height = 1018, width = 1280 )
    coef_plot_add.r <- plot_coefCC(coef_scan_add.rna, map, scan1_output = lod.rna.addcovariatesonly)
    dev.off()
    #cat("Done with ", gene_name, "_chr", chromosome, "\n")
  }
  
  if (chromosome != "20" && diff_diet >=4) {
    # Scan
    coef_scan_add.rna <-  scan1coef(genoprobs = probs[mice.with.data.r, chromosome],
                                kinship=K.LOCO.192.chr.r[chromosome],
                                pheno=expr.rna.192[mice.with.data.r,index.rna],
                                addcovar=addcovar.r[,-1],
                                reml=TRUE )
    save(coef_scan_add.rna, file = paste(path, "Data/Output/Scans/eQTL/", gene_name, "_chr", chromosome, "_Add_coef_scan.Rdata", sep = ""))
    
    
    coef_scan_diet_int.rna <- scan1coef(genoprobs=probs[mice.with.data.r, chromosome],
                                    kinship=K.LOCO.192.chr.r[chromosome],
                                    pheno=expr.rna.192[mice.with.data.r,index.rna],
                                    addcovar=addcovar.r[,-1],
                                    intcovar=intcovar.diet.r[,-1],
                                    reml=TRUE)
    save(coef_scan_diet_int.rna, file = paste(path, "Data/Output/Scans/eQTL/", gene_name, "_chr", chromosome, "_Diet_coef_scan.Rdata", sep = ""))

    # Plot
    png(filename= paste(path, "Data/Output/Plots/Founder_effects/", type, "/", gene_name, "_chr", chromosome, "_diff_diet_", type, sep =""), height = 1018, width = 1280 )
    coef_plot_diet_int.r <- plot_coefCC(coef_scan_diet_int.rna, map, scan1_output = lod.rna.dietint )
    dev.off()
    #cat("Done with ", gene_name, "_chr", chromosome, "\n")
    
  }
  
  if (chromosome != "20" && diff_sex >= 4){
    # Scan / Protein / Interactive covariate: Sex
    coef_scan_add.rna <-  scan1coef(genoprobs = probs[mice.with.data.r, chromosome],
                                kinship=K.LOCO.192.chr.r[chromosome],
                                pheno=expr.rna.192[mice.with.data.r,index.rna],
                                addcovar=addcovar.r[,-1],
                                reml=TRUE )
    save(coef_scan_add.rna, file = paste(path, "Data/Output/Scans/eQTL/", gene_name, "_chr", chromosome, "_Add_coef_scan.Rdata", sep = ""))

    coef_scan_sex_int.rna <- scan1coef(genoprobs=probs[mice.with.data.r, chromosome],
                            kinship=K.LOCO.192.chr.r[chromosome],
                            pheno=expr.rna.192[mice.with.data.r,index.rna],
                            addcovar=addcovar.r[,-1],
                            intcovar=intcovar.sex.r[,-1],
                            reml=TRUE)
    save(coef_scan_sex_int.rna, file = paste(path, "Data/Output/Scans/eQTL/", gene_name, "_chr", chromosome, "_Sex_coef_scan.Rdata", sep = ""))
    

    #Plot
    png(filename= paste(path, "Data/Output/Plots/Founder_effects/", type, "/", gene_name, "_chr", chromosome, "__diff_sex_", type, sep =""), height = 1018, width = 1280 )
    coef_plot_sex_int.r <- plot_coefCC(coef_scan_sex_int.rna, map, scan1_output = lod.rna.sexint)
    dev.off()
    #cat("Done with ", gene_name, "_chr", chromosome, "\n")
   
  }
  
}
```

``` r
#### Define functions below before runnning ####


# For each row in  the distinct peaks df, determine id it is a pQTL or an eQTL and then run the coef scan funtion 
for (i in 1:nrow(distinct_peaks_data)) {
  chr <- distinct_peaks_data$chr[i]
  diff.d <- abs(distinct_peaks_data$diff_diet[i])
  diff.s <- abs(distinct_peaks_data$diff_sex[i])
  gene <- distinct_peaks_data$Gene[i]
  
  # if the type of QTL is pQTL run the protein coef scan plot function
  if (distinct_peaks_data$Type[i] == "pQTL") {
    # print("pQTL is running")
    
    # pass the QTL chromosome, the type, the diet difference, the sex difference, and the gene name to the function
    coef_scan_plot.prot(chr, "Protein", diff.d , diff.s, gene)
      }
      
  # if the type of QTL is eQTL run the protein coef scan plot function
  if (distinct_peaks_data$Type[i] == "eQTL") {
    # print("eQTL is running")
    
    # pass the QTL chromosome, the type, the diet difference, the sex difference, and the gene name to the function
    coef_scan_plot.rna(chr, "RNA", diff.d, diff.s, gene)
      }
      
}
```

### COEF \#\#s

``` r
# Read in scan files 


coef_scan_path <- "/Users/c-allanb/Desktop/DO192Liver/Data/Output/Scans"

# create lists of the file names
coef_scan_files.prot <- list.files(paste(coef_scan_path, "/pQTL", sep = ""))
coef_scan_files.rna <- list.files(paste(coef_scan_path, "/eQTL", sep = ""))

# read in peak data 
file_path <- "/Users/c-allanb/Desktop/DO192Liver/Data/Output/peak_data/diff_data"
peaks_diff_data <- read.table(file_path, header = TRUE, sep = "\t")

# filter the data for QTLs not on chr 20 
peaks_diff_data <- peaks_diff_data %>%
  filter(chr != 20 )
```

    ## Warning: package 'bindrcpp' was built under R version 3.4.4

``` r
# Create empty df to store coefs in
coef_df <- data.frame(matrix(nrow = nrow(peaks_diff_data), ncol = 8)) %>%
  rename(A = X1, B = X2, C= X3, D = X4, E = X5, F = X6, G = X7, H = X8) %>%
  mutate(A = 0, B = 0, C = 0, D = 0, E = 0, F = 0, G = 0, H = 0)

# bind columns
peaks_coef_num_df <- bind_cols(peaks_diff_data, coef_df)

# add columns to peak data with founder coefs
peaks_coef_num_df <- peaks_coef_num_df %>%
  rowwise() %>%
  mutate(
    A = find_peaks_coef("A", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    B = find_peaks_coef("B", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    C = find_peaks_coef("C", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    D = find_peaks_coef("D", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    E = find_peaks_coef("E", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    "F" = find_peaks_coef("F", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    G = find_peaks_coef("G", Gene_Name, Type, chr, QTL_Marker_pos, Covariates),
    H = find_peaks_coef("H", Gene_Name, Type, chr, QTL_Marker_pos, Covariates)
  )
```

``` r
# Create function to find the value of the coef for the position where the QTL is 

find_peaks_coef <- function(founder, name, type, chrom, marker, cov){
  
  # load scan needed 
  #print("test")
  #print(as.character(name))
  #print(as.character(type))
  #print(chrom)
  #print(marker)
  #print(cov)
  coef_scan_path <- paste(coef_scan_path, "/", as.character(type), "/", as.character(name), "_chr", chrom, "_Add_coef_scan.Rdata", sep = "" )
  coef_scan <- get(load(coef_scan_path))
  #print("TEST")

  # for each founder retrieve the coef from the coef scan 
  if (founder == "A") {
    a_founder_coef <- coef_scan[as.character(marker), 1] 
    return(as.numeric(a_founder_coef))
  }
  if (founder == "B") {
    b_founder_coef <- coef_scan[as.character(marker), 2]
    return(b_founder_coef)
  }  
  if (founder == "C") {
    c_founder_coef <- coef_scan[as.character(marker), 3]
    return(c_founder_coef)
  }
  if (founder == "D") {
    d_founder_coef <- coef_scan[as.character(marker), 4]
    return(d_founder_coef)
  }
  if (founder == "E") {
    e_founder_coef <- coef_scan[as.character(marker), 5]
    return(e_founder_coef)
  }
  if (founder == "F") {
    f_founder_coef <- coef_scan[as.character(marker), 6]
    return(f_founder_coef)
  }
  if (founder == "G") {
    g_founder_coef <- coef_scan[as.character(marker), 7]
    return(g_founder_coef)
  }
  if (founder == "H") {
    h_founder_coef <- coef_scan[as.character(marker), 8]
    return(h_founder_coef)
  }  
}
```

### Save the data frame

``` r
# Write the table to an output file 
write.table(peaks_coef_num_df, file = paste(coef_scan_path, "/coef_data", sep = ""), sep = "\t", row.names = FALSE)
```
