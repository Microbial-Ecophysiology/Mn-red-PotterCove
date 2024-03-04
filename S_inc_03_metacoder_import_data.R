# This script imports sequence data produced by the dada2 pipeline for use with metacoder into R

# packages ####
library(tidyverse)
library(taxa)
library(metacoder)


# import data ####
setwd("WORKINGDIRECTORY")

## path to scripts to load
sloc <- "SCRIPTLOCATION/Scripts_to_source/"

## project name to name files
prj <- "PROJECT"


# ASV table
asv_counts <- read.table(paste0("seq_table_for_metacoder_", prj, ".txt"))

# taxonomy table
tax_table <- read.table(paste0("tax_table_for_metacoder_", prj, ".txt"), sep = "\t")

# metadata
mdata <- read_tsv("MDATAFILE.txt")
sample.names <- mdata$Name.Mapping.File


# create taxmap object ####
obj <- parse_dada2(seq_table = asv_counts,
                   tax_table = tax_table)
## remove samples which are not in sample.names
obj$data$asv_table <- obj$data$asv_table[c("taxon_id", "asv_id", sample.names)]

ini_reads <- sum(obj$data$asv_table[, sample.names])

# process original table ####
## remove all none 'Bacteria' sequences
obj <- filter_taxa(obj,
                   taxon_names == "Bacteria",
                   subtaxa = T)

sum(obj$data$asv_table[, sample.names])
sum(obj$data$asv_table[, sample.names])/ini_reads


## remove unwanted taxa Chloroplasts and Mitochondria
obj <- filter_taxa(obj, taxon_names == "Mitochondria", invert = TRUE, 
                   subtaxa = TRUE, reassign_obs = FALSE)
sum(obj$data$asv_table[, sample.names])
sum(obj$data$asv_table[, sample.names])/ini_reads


obj <- filter_taxa(obj, taxon_names == "Chloroplast", invert = TRUE, 
                   subtaxa = TRUE, reassign_obs = FALSE)
sum(obj$data$asv_table[, sample.names])
sum(obj$data$asv_table[, sample.names])/ini_reads

obj



## remove dubletons and singletons
obj$data$asv_table <- zero_low_counts(obj, "asv_table",
                                       min_count = 3,
                                       use_total = T,
                                       other_cols = T)


## check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)

# remove empty ASVs
obj <- filter_obs(obj, "asv_table",
                  ! no_reads,
                  drop_taxa = T)
obj

sum(obj$data$asv_table[, sample.names]) 
sum(obj$data$asv_table[, sample.names])/ini_reads


print(taxonomy_table(obj))


# calculate further tables ####
## Calculate relative abundance
# relative abundance per ASV
obj$data$rel_abd <- calc_obs_props(obj, "asv_table", other_cols = T)
# relative abundance per taxon
obj$data$rel_tax_abd <- calc_taxon_abund(obj, "rel_abd")
print(obj)

# save taxmap object for other scripts
save(obj, file = paste0("taxmap_", prj, ".RData"))




## export asv table ####
export_name <- prj

# function to export OTU table from taxmap
as.otu.table <- function(object, filename) {
  o1 <- as.data.frame(object$data$save)
  o2 <- as.data.frame(classifications(object))
  o2$taxon_id <- rownames(o2)
  o3 <- dplyr::left_join(o1, o2)
  write.table(o3, filename, sep = "\t", row.names = FALSE)
  otu_tax <<- o3
}

# name object table as object 'save' for exporting
obj$data$save <- obj$data$asv_table
as.otu.table(obj, paste0("asv_count_table_curated_", export_name, ".txt"))

obj$data$save <- obj$data$rel_abd
as.otu.table(obj, paste0("rel_abd_table_", export_name, ".txt"))

obj$data$save <- obj$data$rel_tax_abd
as.otu.table(obj, paste0("rel_tax_abd_table_", export_name, ".txt"))
