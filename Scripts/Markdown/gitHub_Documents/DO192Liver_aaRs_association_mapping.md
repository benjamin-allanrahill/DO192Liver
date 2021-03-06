DO192Liver aaRS SNP association
================
Ben Allan-Rahill

-   [DO192 Liver aaRS SNP Asssociation](#do192-liver-aars-snp-asssociation)
    -   [Load libraries](#load-libraries)
-   [Plot](#plot)

DO192 Liver aaRS SNP Asssociation
---------------------------------

``` r
# Download files 
download.file("https://ndownloader.figshare.com/files/9746485", "cc_variants.sqlite")
download.file("https://ndownloader.figshare.com/files/9746452", "mouse_genes_mgi.sqlite")
download.file("https://ndownloader.figshare.com/files/9746452", "mouse_mgi.sqlite")
```

### Load libraries

``` r
library(qtl2)
library(tidyverse)
```

``` r
# Load DO192 Liver Global Data 
load("/Users/c-allanb/Desktop/DO192Liver/Data/Input/DO192LiverData_Formattedfor_rQTL2.Rdata")
```

``` r
# Create querying function
path <- "/Users/c-allanb/Desktop/DO192Liver/Data/Input/snp_association"
query_variants <- create_variant_query_func(paste(path,"/cc_variants.sqlite", sep = ""))
query_genes <- create_gene_query_func(paste(path,"/mouse_mgi.sqlite", sep = ""))
```

``` r
# Import peak data and filter based off LOD = 7 
coef_scan_data.p <- as_tibble(read.table("/Users/c-allanb/Desktop/DO192Liver/Data/Output/Scans/pQTL_coef_data", sep = "\t", header = TRUE)) %>%
  filter(lod >= 7)
```

    ## Warning: package 'bindrcpp' was built under R version 3.4.4

``` r
coef_scan_data.r <- as_tibble(read.table("/Users/c-allanb/Desktop/DO192Liver/Data/Output/Scans/eQTL_coef_data", sep = "\t", header = TRUE)) %>%
  filter(lod >= 7)
```

``` r
# Specifiy the gene that you want to run SNP association on 
gene_name <- "Rars"
```

``` r
# Get the confidence interval for the LOD peak at that gene 
interval_beg.p <- coef_scan_data.p$ci_lo[coef_scan_data.p$Gene == gene_name]

interval_end.p <- coef_scan_data.p$ci_hi[coef_scan_data.p$Gene == gene_name]

interval_beg.r <- coef_scan_data.r$ci_lo[coef_scan_data.r$Gene == gene_name]

interval_end.r <- coef_scan_data.r$ci_hi[coef_scan_data.r$Gene == gene_name]
```

``` r
# Establish the chromosome of the varaint you want to look at 
chr.p <- coef_scan_data.p$chr[coef_scan_data.p$Gene == gene_name]
chr.r <- coef_scan_data.r$chr[coef_scan_data.r$Gene == gene_name]
```

``` r
# Get the varaints (SNPs, indels, SV) in the peak interval 
variants.p <- query_variants(chr.p, interval_beg.p, interval_end.p)
variants.r <- query_variants(chr.r, interval_beg.r, interval_end.r)
```

``` r
##PROTEIN##
# Get the gene's index
index.p <- which(annotations.protein.192$Associated.Gene.Name == gene_name)
print(index.p)
```

    ## [1] 4925

``` r
mice.with.data.p <- !is.na(expr.protein.192[,index.p])

# kinship 
K.LOCO.192.chr.p <- K.LOCO.qtl2


for (kloco.chr in 1:20){
  K.LOCO.192.chr.p[[chr.p]]<- K.LOCO.qtl2[[chr.p]][mice.with.data.p,mice.with.data.p]
}

# Set up additive and interactive covariates for RNA, only mice with data.
addcovar.p <- model.matrix(~ Sex + Diet, data=covariates.protein.192[mice.with.data.p, ])
intcovar.sex.p <- model.matrix(~ Sex, data=covariates.protein.192[mice.with.data.p, ])
intcovar.diet.p <- model.matrix(~ Diet, data=covariates.protein.192[mice.with.data.p, ])
```

``` r
##RNA##
# Get the gene's index
index.r <- which(annotations.rna.192$Gene == gene_name)
print(index.r)
```

    ## [1] 13118

``` r
mice.with.data.r <- !is.na(expr.rna.192[,index.r])

# kinship 
K.LOCO.192.chr.r <- K.LOCO.qtl2


for (kloco.chr in 1:20){
  K.LOCO.192.chr.r[[chr.r]]<- K.LOCO.qtl2[[chr.r]][mice.with.data.r,mice.with.data.r]
}

# Set up additive and interactive covariates for RNA, only mice with data.
addcovar.r <- model.matrix(~ Sex + Diet, data=covariates.rna.192[mice.with.data.r, ])
intcovar.sex.r <- model.matrix(~ Sex, data=covariates.rna.192[mice.with.data.r, ])
intcovar.diet.r <- model.matrix(~ Diet, data=covariates.rna.192[mice.with.data.r, ])
```

``` r
# get variant LOD scores 
snps.p <- scan1snps(
                    probs, #genome probabilities  
                    map, #physical map
                    pheno=expr.protein.192[mice.with.data.p,index.p], # matrix of phenotypes
                    K.LOCO.qtl2[[paste(chr.p)]], #kinship matrix for the chr
                    addcovar = addcovar.p, # matrix for additive covariates 
                    query_func=query_variants, 
                    chr=chr.p, 
                    start=interval_beg.p, 
                    end=interval_end.p, 
                    keep_all_snps=TRUE
                    )

snps.r <- scan1snps(
                    probs, #genome probabilities  
                    map, #physical map
                    pheno=expr.rna.192[mice.with.data.r,index.r], # matrix of phenotypes
                    K.LOCO.qtl2[[paste(chr.r)]], #kinship matrix for the chr
                    addcovar = addcovar.r, # matrix for additive covariates 
                    query_func=query_variants, 
                    chr=chr.r, 
                    start=interval_beg.r, 
                    end=interval_end.r, 
                    keep_all_snps=TRUE
                    )
```

Plot
----

``` r
##PRTOETIN##
genes.p <- query_genes(chr.p, interval_beg.p, interval_end.p)
par(mar=c(4.1, 4.1, 0.6, 0.6))

png(filename = paste("/Users/c-allanb/Desktop/DO192Liver/Data/Output/Plots/Association_plots/DO192Liver_aaRS_assocaition_plots_",gene_name, "chr", chr.p, "pQTL", sep = "" ), height = 1280, width = 1920)

plot(snps.p$lod, snps.p$snpinfo, drop_hilit=1.5, genes=genes.p, main = paste("SNP_Association", gene_name, "chr", chr.p, "pQTL", sep = "_"))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
plot(snps.p$lod, snps.p$snpinfo, drop_hilit=1.5, genes=genes.p, main = paste("SNP_Association", gene_name, "chr", chr.p, "pQTL", sep = "_"))
```

![](DO192Liver_aaRs_association_mapping_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
##RNA##
genes.r <- query_genes(chr.p, interval_beg.p, interval_end.p)
par(mar=c(4.1, 4.1, 0.6, 0.6))

png(filename = paste("/Users/c-allanb/Desktop/DO192Liver/Data/Output/Plots/Association_plots/DO192Liver_aaRS_assocaition_plots_",gene_name, "chr", chr.p, "eQTL", sep = "" ), height = 1280, width = 1920)

plot(snps.r$lod, snps.r$snpinfo, drop_hilit=1.5, genes=genes.r, main = paste("SNP_Association", gene_name, "chr", chr.r, "eQTL", sep = "_"))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
plot(snps.r$lod, snps.r$snpinfo, drop_hilit=1.5, genes=genes.r, main = paste("SNP_Association", gene_name, "chr", chr.r, "eQTL", sep = "_"))
```

![](DO192Liver_aaRs_association_mapping_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
# find top SNPs 
top.p <- top_snps(snps.p$lod, snps.p$snpinfo)
print(top.p[5,c(1, 8:15, 20)], row.names=FALSE)
```

    ##       snp_id A_J C57BL_6J 129S1_SvImJ NOD_ShiLtJ NZO_HlLtJ CAST_EiJ
    ##  rs251435921   1        1           1          1         1        2
    ##  PWK_PhJ WSB_EiJ      lod
    ##        1       1 7.670241

``` r
top.r <- top_snps(snps.r$lod, snps.r$snpinfo)
print(top.r[5,c(1, 8:15, 20)], row.names=FALSE)
```

    ##       snp_id A_J C57BL_6J 129S1_SvImJ NOD_ShiLtJ NZO_HlLtJ CAST_EiJ
    ##  rs236393446   1        1           1          1         1        2
    ##  PWK_PhJ WSB_EiJ      lod
    ##        1       1 8.165269
