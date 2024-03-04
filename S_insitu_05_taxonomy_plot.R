# This script imports taxmap objects created in script metacoder_import_data.R produced through dada2 pipeline

## produce stacked barplots for overview
# follow stackoverflow question https://stackoverflow.com/questions/62627480/how-to-creat-a-bar-graph-of-microbiota-data-with-one-color-for-higher-taxonomic

# packages ####
library(tidyverse)
library(taxa)
library(metacoder) # version 0.3.6
library(phyloseq)
library(hues)
library(cowplot)
library(devEMF)
library(patchwork)
library(viridis)

# import data ####
setwd("WORKINGDIRECTORY")

## path to scripts to load
sloc <- "SCRIPTLOCATION/Scripts_to_source/"

## project name to name files
prj <- "PROJECT"

## taxmap object from metacoder_import_data R script, samples with insufficient coverage removed
load(paste0("taxmap_", prj, ".RData"))

## metadata
mdata_all <- read_tsv("MDATAFILE.txt")



# create phyloseq object ####
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "asv_id", sample_data = mdata, sample_id_col = "Sequencing-ID-Bac")

ntaxa(objphy)                  # how many taxa (including subtaxa)
nsamples(objphy)               # how many samples
sample_names(objphy)           # check sample names
sample_variables(objphy)       # sample variables from mdata
otu_table(objphy)              # ASV count table
tax_table(objphy)              # taxonomy for ASVs
taxa_names(objphy)             # ASV IDs

rank_names(objphy)             # check rank names (need to be Kingdom, Phylum etc.)

## calculate relative abundance
opr <- transform_sample_counts(objphy, function(x) x/sum(x))

## melt phyloseq object into dataframe
phylo_melt <- psmelt(opr)

### save object phylo_melt used to perform all following operations and remove not needed objects
save(phylo_melt, file = paste0(prj, "_phylomelt.RData"))
rm(obj, objphy, opr)
load(paste0(prj, "_phylomelt.RData"))

# rename taxa depending on max abundance ####
## find correct threshold to include selected taxa
phylo_melt %>% 
  filter(Family == "Arcobacteraceae") %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  filter(Abundance == max(Abundance)) %>% 
  pull(Abundance)
# 0.0005513779

phylo_melt %>% 
  filter(Genus == "Desulfuromusa") %>% 
  group_by(Sample) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  filter(Abundance == max(Abundance)) %>% 
  pull(Abundance)
# 0.001743044

## use function 'sort_abundant_taxa'
source(paste0(sloc, "filter_taxa_above_threshold.R"))

## taxa above 5% in at least one samples
# ta5 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 5, Abundance = "Abundance", Sample = "Sample")
# 3 genera, 9 families, 10 orders, 10 classes, 7 phyla

## taxa above 2% in at least one samples
# ta2 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 2, Abundance = "Abundance", Sample = "Sample")
# 13 genera, 18 families, 19 orders, 15 classes, 13 phyla

## taxa above 4% in at least one samples - same threshold as complete in situ plot
# ta4 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 4, Abundance = "Abundance", Sample = "Sample")
# # 21 genera, 30 families, 27 orders, 18 classes, 12 phyla

## taxa above 0.1% in at least one samples
# ta0.1 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.1, Abundance = "Abundance", Sample = "Sample")
# # 112 genera, 106 families, 97 orders, 68 classes, 39 phyla

## taxa above 0.05% in at least one samples - needed for Mn-MS plot
# ta0.05 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 0.05, Abundance = "Abundance", Sample = "Sample")
# # 148 genera, 137 families, 130 orders, 73 classes, 40 phyla


# rm(ta2, ta5, ta0.1)

# MS supl figure overview class 4% ####
## decide on rank level to plot: Genus and sorted by Class
unique(ta4$ASV_table_taxa_abv_4_perc$Class_mod) # 18 taxa to sort by -> use this
unique(ta4$ASV_table_taxa_abv_4_perc$Order_mod) # 29 taxa to sort by
unique(ta4$ASV_table_taxa_abv_4_perc$Phylum_mod) # 8 taxa

ta4_cl <- select(ta4$ASV_table_taxa_abv_4_perc, OTU, Sample, Phylum_mod, Class_mod, Abundance, Station, Depth..cm., Core) %>% 
  mutate(Phylum_mod = factor(Phylum_mod, levels = c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Desulfobacterota", 
                                                    "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "other_Bacteria_<4%")),
         Class_mod=factor(Class_mod, levels=c("other_p_Acidobacteriota_<4%",
                                              "Acidimicrobiia", "other_p_Actinobacteriota_<4%",
                                              "Bacteroidia", "other_p_Bacteroidota_<4%",
                                              "Desulfobacteria", "Desulfobulbia", "Desulfuromonadia", "other_p_Desulfobacterota_<4%",
                                              "Phycisphaerae", "Planctomycetes", "other_p_Planctomycetota_<4%",
                                              "Alphaproteobacteria", "Gammaproteobacteria", "other_p_Proteobacteria_<4%",
                                              "Verrucomicrobiae", "other_p_Verrucomicrobiota_<4%",
                                              "other_Bacteria_<4%"),
                          labels = c("p_Acidobacteriota",
                                     "Acidimicrobiia", "other_p_Actinobacteriota_<4%",
                                     "Bacteroidia", "other_p_Bacteroidota_<4%",
                                     "Desulfobacteria", "Desulfobulbia", "Desulfuromonadia", "other_p_Desulfobacterota_<4%",
                                     "Phycisphaerae", "Planctomycetes", "other_p_Planctomycetota_<4%",
                                     "Alphaproteobacteria", "Gammaproteobacteria", "other_p_Proteobacteria_<4%",
                                     "Verrucomicrobiae", "other_p_Verrucomicrobiota_<4%",
                                     "other_Bacteria_<4%"))
  ) %>%
  group_by(Class_mod)    # reorder classes so classes of same phylum are together



## sum ASV abundance by taxon ####
## need to decide which rank to plot and to group by which other rank
# here plot genus, group by class
ta4_cl_s <- ta4_cl %>% 
  group_by(Sample, Phylum_mod, Class_mod, Depth..cm.) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup()

save(ta4_cl, file = paste0(prj, "_ta4_cl_class_4perc.RData"))
save(ta4_cl_s, file = paste0(prj, "_ta4_cl_summed_class_4perc.RData"))
load(paste0(prj, "_ta4_cl_class_4perc.RData"))


## plot taxa ####
### get automatic colors ####
# function to define multiple colour pallete based on subgroups from here
# https://stackoverflow.com/questions/49818271/stacked-barplot-with-colour-gradients-for-each-bar?rq=1
ColourPalleteMulti <- function(df, group, subgroup){
  
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 40)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}

## create colour gradient with color pallet function
colours <- ColourPalleteMulti(ta4_cl_s, "Phylum_mod", "Class_mod")
### import hand edited color table
colours_edit <- scan(paste0(prj, "_class_color-pallet-edit.txt"), what = "")


### plot class < 4% ####
plot_class <- ggplot(ta4_cl_s, aes(x = Depth..cm., y = Abundance*100, 
                                   fill = Class_mod, color = Class_mod)) + 
  geom_bar(position = position_stack(reverse = T), stat = "identity", 
           width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  coord_flip() +
  scale_x_reverse() +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_fill_manual("", values = colours_edit) +
  scale_color_manual("", values = colours_edit) +
  guides(fill = guide_legend(ncol = 1), color = "none") +
  labs(x = "Depth (cm)", y = "Relative abundance (%)", 
       title = "All classes > 4%") +
  theme_cowplot() +
  theme(legend.position = "right",
        text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_line(linewidth = 0.5), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 9))

ggsave("plots/PC_insitu_STA01.02_class_4perc.pdf", plot_class,
       width = 18, height = 12, units = "cm")
ggsave("plots/PC_insitu_STA01.02_class_4perc.jpg", plot_class,
       width = 18, height = 12, units = "cm")
ggsave("plots/PC_insitu_STA01.02_class_4perc.emf", plot_class,
       width = 18, height = 12, units = "cm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


# MS supl figure detail enriched taxa ####
### only taxa of interest for Mn-MS paper:
# Desulfuromonadia, Campylobacteria
# use 0.05% table
data_Mn <- subset(ta0.05$ASV_table_taxa_abv_0.05_perc, 
                  ta0.05$ASV_table_taxa_abv_0.05_perc$Class %in% 
                    c("Campylobacteria", "Desulfuromonadia"))
unique(data_Mn$Genus_mod) # 13
unique(data_Mn$Family_mod) # 10

## rename taxa and sum ####
data_Mn <- subset(ta0.05$ASV_table_taxa_abv_0.05_perc, 
                  ta0.05$ASV_table_taxa_abv_0.05_perc$Class %in% 
                    c("Campylobacteria", "Desulfuromonadia"))  %>% 
  select(OTU, Sample, Family_mod, Genus_mod, Class, Abundance, Depth..cm.) %>% 
  mutate(Family_mod = factor(Family_mod, levels = c("Arcobacteraceae", "Sulfurimonadaceae", "Sulfurovaceae", "other_o_Campylobacterales_<0.05%",
                                                    "other_o_Bradymonadales_<0.05%", "Desulfuromonadaceae", "Geopsychrobacteraceae", "Sva1033",
                                                    "other_o_Desulfuromonadales_<0.05%", "other_c_Desulfuromonadia_<0.05%")),
         Genus_mod = factor(Genus_mod, levels=c("other_f_Arcobacteraceae_<0.05%", "Sulfurimonas", "Sulfurovum", 
                                                "other_o_Campylobacterales_<0.05%",
                                                
                                                "other_o_Bradymonadales_<0.05%", 
                                                "Desulfuromonas", "other_f_Desulfuromonadaceae_<0.05%",
                                                "Desulfuromusa", "Geopsychrobacter", "other_f_Geopsychrobacteraceae_<0.05%",
                                                "other_f_Sva1033_<0.05%",
                                                "other_o_Desulfuromonadales_<0.05%", "other_c_Desulfuromonadia_<0.05%" ),
                          labels = c("f_Arcobacteraceae", "Sulfurimonas", "Sulfurovum", 
                                     "other_o_Campylobacterales_<0.05%",
                                     
                                     "o_Bradymonadales", 
                                     "Desulfuromonas", "other_f_Desulfuromonadaceae_<0.05%",
                                     "Desulfuromusa", "Geopsychrobacter", "other_f_Geopsychrobacteraceae_<0.05%",
                                     "f_Sva1033",
                                     "other_o_Desulfuromonadales_<0.05%", "other_c_Desulfuromonadia_<0.05%" ))
  ) %>%
  group_by(Genus_mod)    # reorder genera so genera of same family are together



data_Mn_sum <- data_Mn %>% 
  group_by(Sample, Family_mod, Genus_mod, Class, Depth..cm.) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup()

## separate for Desulfuromonadia and Campylobacteria ####
data_Desu <- subset(data_Mn_sum, data_Mn_sum$Class == "Desulfuromonadia")
data_Camp <- subset(data_Mn_sum, data_Mn_sum$Class == "Campylobacteria")

## plot ####

### Desulfuromonadia
# function to define multiple colour pallete based on subgroups from here
# https://stackoverflow.com/questions/49818271/stacked-barplot-with-colour-gradients-for-each-bar?rq=1
ColourPalleteMulti <- function(df, group, subgroup){
  
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 100)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 40)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}

## create colour gradient with color pallet function
colours_Desu <- ColourPalleteMulti(data_Desu, "Family_mod", "Genus_mod")
write(colours_Desu, "plots/MS-Mn_PC_insitu-reseq_Desul_color-pallet-auto.txt")
### import hand edited color table
colours_Desu_edit <- scan(paste0(prj, "_color-pallet-edit.txt"), what = "")

plot_Desu <- ggplot(data_Desu, aes(x = Depth..cm., y = Abundance*100, 
                                   fill = Genus_mod, color = Genus_mod)) + 
  geom_bar(position = position_stack(reverse = T), stat = "identity", 
           width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  coord_flip() +
  scale_x_reverse() +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_fill_manual("", values = colours_Desu_edit) +
  scale_color_manual("", values = colours_Desu_edit) +
  guides(fill = guide_legend(ncol = 1), color = "none") +
  labs(x = "Depth (cm)", y = "Relative abundance (%)", 
       title = "Desulfuromonadia") +
  theme_cowplot() +
  theme(legend.position = "right",
        text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_line(linewidth = 0.5), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 9))


### Campylobacteria
plot_Camp <- ggplot(data_Camp, aes(x = Depth..cm., y = Abundance*100, 
                                   fill = Genus_mod, color = Genus_mod)) + 
  geom_bar(position = position_stack(reverse = T), stat = "identity", 
           width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  coord_flip() +
  scale_x_reverse() +
  scale_y_continuous(position = "right", expand = c(0,0)) +
  scale_fill_viridis("", discrete = T) +
  scale_color_viridis("", discrete = T) +
  guides(fill = guide_legend(ncol = 1), color = "none") +
  labs(x = "Depth (cm)", y = "Relative abundance (%)", 
       title = "Campylobacteria") +
  theme_cowplot() +
  theme(legend.position = "right",
        text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_line(linewidth = 0.5), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 9))


# combine supl Fig in situ community ####
supl_fig <- (plot_class + guides(fill = guide_legend(ncol = 3)) + theme(legend.position = "bottom")) / 
  ((plot_Desu + guides(fill = guide_legend(ncol = 2)) + 
      theme(legend.position = "bottom", legend.box.just = "top", legend.box.margin = margin(t = 0, l = -1, unit = "cm"))) + 
     (plot_Camp + theme(legend.position = "bottom", legend.box.just = "top", legend.box.margin = margin(t = -0.5, l = 1, unit = "cm")))
   + plot_layout(ncol = 2, widths = c(1, 0.7))) &
  plot_annotation(tag_levels = "A", title = "Bacterial 16S rRNA gene community STA01.02", 
                  theme = theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold")))

ggsave("plots/MS_supl_FigSX_insitu_STA01.02.jpg", supl_fig, width = 180, height = 200, units = "mm")
ggsave("plots/MS_supl_FigSX_insitu_STA01.02.pdf", supl_fig, width = 180, height = 200, units = "mm")
ggsave("plots/MS_supl_FigSX_insitu_STA01.02.emf", supl_fig, width = 180, height = 200, units = "mm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
