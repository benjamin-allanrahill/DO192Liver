### Scan, Collate, and Analyze
### June 15th, 2018
### Doug Perkins

#### OVERVIEW ####
# This script allows you to describe your Diversity Outbred dataset and then have it automatically
# scanned, collated, and analyzed, according to whichever parameters and analyses you select.
# Scans will be output in folders called rna_scan, protein_scan, and rna_protein_scan in the output directory you choose. 
# Collated tables will be output in folders called rna_collate, etc. with subfolders for each interactive covariate.

#### Choose your input, output, and scripts directories. ####
dir.input <- "/projects/munger-lab/DO192Liver/new/scan/qtlscan/data/input"
dir.output <- "/projects/munger-lab/DO192Liver/new/scan/qtlscan/data/output"
dir.scripts.r <- "/projects/munger-lab/DO192Liver/new/scan/qtlscan/scripts/r"
dir.scripts.python <- "/projects/munger-lab/DO192Liver/new/scan/qtlscan/scripts/python"

#### Provide some variable names to be used in output filenames. ####
dataset.name <- "DO192Liver"

# Install packages, load packages, and load the functions.
setwd(dir.scripts.r)
source("pipeline_functions_partial_2018_06_26.R")

options(stringsAsFactors = F)
setwd(dir.input)
load("DO192LiverData_Formattedfor_rQTL2.Rdata")

#### Prepare generic data frames ####
# Tell the functions how to interpret your data by pointing them to each type of data (ex: rna.expression) or column name (ex: rna.chromosome.col).
# DO192Liver
generic.data <- make.data.generic(
                  rna.annotations = annotations.rna.192, 
                  rna.expression = expr.rna.192, 
                  rna.covariates = covariates.rna.192, 
                  rna.ENSMUSG.col = "EnsemblID",
                  rna.gene.col = "Gene",
                  rna.chromosome.col = "Chr",
                  rna.start.mbp.col = "Start.Mbp",
                  rna.end.mbp.col = "End.Mbp",
                  rna.gene.biotype.col = "Gene.Biotype",
                  
                  protein.annotations = annotations.protein.192, 
                  protein.expression = expr.protein.192, 
                  protein.covariates = covariates.protein.192, 
                  protein.ENSMUSG.col = "Ensembl.Gene.ID", 
                  protein.ENSMUST.col = "Ensembl.Transcript.ID", 
                  protein.ENSMUSP.col = "Ensembl.Protein.ID",
                  protein.gene.col = "Associated.Gene.Name",
                  protein.chromosome.col = "Chromosome.Name",
                  protein.start.mbp.col = "Gene.Start..bp.",
                  protein.end.mbp.col = "Gene.End..bp.",
                  protein.gene.biotype.col = "Gene.Biotype",
                  
                  kinship = K.LOCO.qtl2, 
                  founder.probs = probs, 
                  
                  markers = markers.64K,
                  markers.id.col = "marker",
                  markers.chromosome.col = "chr",
                  markers.mbp.col = "bp",
                  markers.cM.col = "cM",
                  markers.bp.col = "pos")

# Load in the new generic data frames.
make.generic.data.structures()
rm(annotations.rna.192); rm(annotations.protein.192); rm(covariates.rna.192); rm(covariates.protein.192); rm(expr.rna.192); rm(expr.protein.192); rm(markers.64K); rm(X.protein); rm(X.rna); rm(K.LOCO.qtl2)

#### Scan, Collate, and Analyze ####

#scan.generic(rna=T, protein=F, rna.with.proteins.only=F, rna.beg=rna.start.num, rna.end=rna.end.num, protein.beg=protein.start.num, protein.end=protein.end.num)
collate.generic(rna=T, protein=T, rna.diet.int=T, rna.sex.int=T, protein.diet.int=T, protein.sex.int=T, lod.thresholds=c(4,5,6,7), window.type="center", peaks.only = T)

#classification()
#manhattan_generic()
#ome_map_generic()
#moving_window_generic()
#coefficient_generic()
#mediation_generic()