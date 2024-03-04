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
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "asv_id", sample_data = mdata, sample_id_col = "Name.Mapping.File")

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
## use function 'sort_abundant_taxa'
source(paste0(sloc, "filter_taxa_above_threshold.R"))

## taxa above 5% in at least one samples
ta5 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 5, Abundance = "Abundance")
# 6 genera, 12 families, 9 orders, 7 classes, 6 phyla

## taxa above 2% in at least one samples
# ta2 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 2, Abundance = "Abundance")
# 12 genera, 24 families, 18 orders, 16 classes, 12 phyla

# use 5% threshold
# rm(ta2, ta3, ta4)

# sort and make nice names for taxa ####
## other columns needed as factor for plotting
cols_factor <- c("Incubation", "Replicate", "Substrate", "NucleicAcid")

## decide on rank level to plot: Genus and sorted by Class
unique(ta5$ASV_table_taxa_abv_5_perc$Class_mod) # 13 taxa to sort by
unique(ta5$ASV_table_taxa_abv_5_perc$Genus_mod) # 33 taxa to plot

ta5n <- select(ta5$ASV_table_taxa_abv_5_perc, OTU, Sample, Class_mod, Genus_mod, Abundance, Incubation, Substrate, Replicate, Timepoint, NucleicAcid, Genus, Family, Class, Order, Phylum) %>% 
  mutate(Class_mod=factor(Class_mod, levels = c("Bacteroidia", "other_p_Bacteroidota_<5%",
                                                "Campylobacteria",
                                                "Desulfobacteria", "Desulfobulbia", "Desulfuromonadia", "other_p_Desulfobacterota_<5%",
                                                "Planctomycetes", "other_p_Planctomycetota_<5%",
                                                "Gammaproteobacteria", "other_p_Proteobacteria_<5%",
                                                "other_p_Verrucomicrobiota_<5%",
                                                "other_Bacteria_<5%")),
         
         Genus_mod=factor(Genus_mod, levels = c("other_o_Bacteroidales_<5%", 
                                                "other_f_Flavobacteriaceae_<5%", "other_o_Flavobacteriales_<5%",
                                                "other_c_Bacteroidia_<5%", "other_p_Bacteroidota_<5%",
                                                
                                                
                                                "other_f_Arcobacteraceae_<5%", "Sulfurimonas",
                                                "other_o_Campylobacterales_<5%",
                                                
                                                
                                                "other_f_Desulfobacteraceae_<5%",
                                                "Desulfofrigus", "other_f_Desulfolunaceae_<5%",
                                                "other_o_Desulfobacterales_<5%", "other_c_Desulfobacteria_<5%",
                                                
                                                "Desulforhopalus", "other_f_Desulfocapsaceae_<5%", 
                                                "other_o_Desulfobulbales_<5%",
                                                
                                                "Desulfuromonas", "other_f_Desulfuromonadaceae_<5%",
                                                "Desulfuromusa", "other_f_Geopsychrobacteraceae_<5%",
                                                "other_f_Sva1033_<5%",
                                                "other_o_Desulfuromonadales_<5%", "other_c_Desulfuromonadia_<5%",
                                                
                                                "other_p_Desulfobacterota_<5%",
                                                
                                                
                                                "other_f_Pirellulaceae_<5%", "other_c_Planctomycetes_<5%", "other_p_Planctomycetota_<5%",
                                                
                                                "other_f_Unknown Family_<5%",
                                                "Woeseia",
                                                "other_c_Gammaproteobacteria_<5%", "other_p_Proteobacteria_<5%",
                                                
                                                "other_p_Verrucomicrobiota_<5%",
                                                "other_Bacteria_<5%"),
                          labels = c("o_Bacteroidales", 
                                     "f_Flavobacteriaceae", "other_o_Flavobacteriales_<5%",
                                     "other_c_Bacteroidia_<5%", "other_p_Bacteroidota_<5%",
                                     
                                     
                                     "f_Arcobacteraceae", "Sulfurimonas",
                                     "other_o_Campylobacterales_<5%",
                                     
                                     
                                     "f_Desulfobacteraceae",
                                     "Desulfofrigus", "other_f_Desulfolunaceae_<5%",
                                     "other_o_Desulfobacterales_<5%", "other_c_Desulfobacteria_<5%",
                                     
                                     "Desulforhopalus", "other_f_Desulfocapsaceae_<5%", 
                                     "other_o_Desulfobulbales_<5%",
                                     
                                     "Desulfuromonas", "other_f_Desulfuromonadaceae_<5%",
                                     "Desulfuromusa", "other_f_Geopsychrobacteraceae_<5%",
                                     "f_Sva1033",
                                     "other_o_Desulfuromonadales_<5%", "other_c_Desulfuromonadia_<5%",
                                     
                                     "other_p_Desulfobacterota_<5%",
                                     
                                     
                                     "f_Pirellulaceae", "other_c_Planctomycetes_<5%", "other_p_Planctomycetota_<5%",
                                     
                                     "Unknown_Family_Gammapr.",
                                     "Woeseia",
                                     "other_c_Gammaproteobacteria_<5%", "other_p_Proteobacteria_<5%",
                                     
                                     "p_Verrucomicrobiota_<5%",
                                     "other_Bacteria_<5%")),
         ASV.taxa = factor(paste0(OTU, "\n", Genus)),  # new column with ASV ID and genus name
         across(all_of(cols_factor), as.factor)) %>%   # reorder families so families of same phylum are together
  group_by(Genus_mod)


# sum ASV abundance by taxon ####
## need to decide which rank to plot and to group by which other rank
 # here plot genus, group by class
ta5s <- ta5n %>% 
  group_by(Sample, Class_mod, Genus_mod, Incubation, Substrate, Replicate, Timepoint, NucleicAcid) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup()

save(ta5n, ta5s, file = paste0(prj, "_ta5n_ta5s_rel.abd.ASV.taxa.long_ab5perc_.RData"))

load(paste0(prj, "_ta5n_ta5s_rel.abd.ASV.taxa.long_ab5perc_.RData"))

# subset RNA or DNA ####
ta5n_R <- subset(ta5n, ta5n$NucleicAcid == "RNA")
ta5n_D <- subset(ta5n, ta5n$NucleicAcid == "DNA")
ta5s_R <- subset(ta5s, ta5s$NucleicAcid == "RNA")
ta5s_D <- subset(ta5s, ta5s$NucleicAcid == "DNA")

# plot RNA level ####
## get automatic colors ####
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
# colours <- ColourPalleteMulti(ta5s, "Class_mod", "Genus_mod")
### export generated color table
# write(colours, paste0(prj, "_color-pallet-automatical.txt"))
### import hand edited color table
colours_edit <- scan(paste0(prj, "_color-pallet-edit.txt"), what = "")


## plot basic taxa ####
ggplot(ta5s_R, aes(x = interaction(Replicate, Timepoint), y = Abundance, fill = Genus_mod, color = Genus_mod)) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual("", values = colours_edit) +
  scale_color_manual("", values = colours_edit) +
  guides(fill = guide_legend(ncol = 5), color = guide_legend(ncol = 5)) +
  labs(y = "Relative abundance") +
  facet_wrap(~Substrate, ncol = 5) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        axis.line = element_line(linewidth = 0), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5))


### sort facets differently
plot_basic <- ggplot(ta5s_R, aes(x = Replicate, y = Abundance*100, fill = Genus_mod, color = Genus_mod)) + # for supplementary plot
# plot_basic <- ggplot(ta2s_R, aes(x = Replicate, y = Abundance, fill = Genus_mod, color = Genus_mod)) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual("", values = colours_edit) +
  scale_color_manual("", values = colours_edit) +
  guides(fill = guide_legend(ncol = 3), color = guide_legend(ncol = 3)) +
  labs(y = "Relative abundance (%)", x = "Time (days)") + # for supplementary plot
  # labs(y = "Relative abundance", x = "Time (days)") +
  facet_wrap(~Timepoint, nrow = 1, strip.position = "bottom") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_line(linewidth = 0.5), 
        strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 9))

#### subset by treatment
plot_birnessite.acetate <- plot_basic %+% subset(ta5s_R, ta5s_R$Substrate == "birnessite.acetate") +
  labs(title = "Birnessite + acetate")

plot_acetate <- plot_basic %+% subset(ta5s_R, ta5s_R$Substrate == "acetate") +
  labs(title = "Acetate")

plot_day0 <- plot_basic %+% subset(ta5s_R, ta5s_R$Timepoint == 0) +
  labs(title = "Slurry")


## MS supl Fig taxonomy RNA ####
MS_supl_plot_RNA <- plot_day0 + plot_birnessite.acetate + plot_acetate + 
  plot_layout(guides = "collect", widths = c(1, 4, 4)) + 
  plot_annotation(tag_levels = "A", title = "Bacterial community RNA level",
                  theme = theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))) &
  theme(legend.position = "bottom")

ggsave("plots/MS_supl_FigSX_bac_barplot_RNA.jpg", plot = MS_supl_plot_RNA, width = 180, height = 140, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_RNA.pdf", plot = MS_supl_plot_RNA, width = 180, height = 140, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_RNA.emf", plot = MS_supl_plot_RNA, width = 180, height = 140, units = "mm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


#
# plot DNA level ####
ggplot(ta5s_D, aes(x = interaction(Replicate, Timepoint), y = Abundance, fill = Genus_mod, color = Genus_mod)) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual("", values = colours_edit) +
  scale_color_manual("", values = colours_edit) +
  guides(fill = guide_legend(ncol = 5), color = guide_legend(ncol = 5)) +
  labs(y = "Relative abundance") +
  facet_wrap(~Substrate, ncol = 5) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        axis.line = element_line(linewidth = 0), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.text = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5))

plot_D_basic <- ggplot(ta5s_D, aes(x = Replicate, y = Abundance*100, fill = Genus_mod, color = Genus_mod)) + # for supplementary figure
# plot_D_basic <- ggplot(ta2s_D, aes(x = Replicate, y = Abundance, fill = Genus_mod, color = Genus_mod)) + 
  geom_bar(position = position_stack(), stat = "identity", width = 0.8, linewidth = 0.01) +  # Stacked 100% barplot
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual("", values = colours_edit) +
  scale_color_manual("", values = colours_edit) +
  # guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) +
  guides(fill = guide_legend(ncol = 3), color = guide_legend(ncol = 3)) + # for supplementary figure
  # labs(y = "Relative abundance", x = "Time (days)") +
  labs(y = "Relative abundance (%)", x = "Time (days)") + # for supplementary figure
  facet_wrap(~Timepoint, nrow = 1, strip.position = "bottom") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_line(linewidth = 0.5), 
        strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 9))

#### subset by treatment
plot_D_birnessite.acetate <- plot_D_basic %+% subset(ta5s_D, ta5s_D$Substrate == "birnessite.acetate") +
  labs(title = "Birnessite + acetate")

plot_D_acetate <- plot_D_basic %+% subset(ta5s_D, ta5s_D$Substrate == "acetate") +
  labs(title = "Acetate")

plot_D_day0 <- plot_D_basic %+% subset(ta5s_D, ta5s_D$Timepoint == 0) +
  labs(title = "Slurry")


## MS supl Fig taxonomy DNA ####
MS_supl_plot_DNA <- plot_D_day0 + plot_D_birnessite.acetate + plot_D_acetate + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A", title = "Bacterial community DNA level",
                  theme = theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))) &
  theme(legend.position = "bottom")

ggsave("plots/MS_supl_FigSX_bac_barplot_DNA.jpg", plot = MS_supl_plot_DNA, width = 180, height = 140, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_DNA.pdf", plot = MS_supl_plot_DNA, width = 180, height = 140, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_DNA.emf", plot = MS_supl_plot_DNA, width = 180, height = 140, units = "mm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})



# plot MS Fig 3 ####
## lineplot, only Mn + acetate, acetate treatment
 # Desulfuromusa, Desulfuromonas, Sva1033, uncl. Arcobacteraceae

## subset data ####
ta5n_MSFig03 <- ta5n %>% 
  filter(Substrate %in% c("acetate", "birnessite.acetate", "DIC") &
           (Genus %in% c("Desulfuromonas", "Desulfuromusa") |
              (Family == "Arcobacteraceae" & is.na(Genus)) |
              Family == "Sva1033") &
           Timepoint != 30)

ta5s_MSFig03 <- ta5n_MSFig03 %>% 
  group_by(Sample, Genus_mod, Substrate, Replicate, Timepoint, NucleicAcid) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(Genus_mod = factor(Genus_mod, levels = c("Desulfuromusa", "f_Sva1033", "Desulfuromonas", "f_Arcobacteraceae"),
                            labels = c("italic(Desulfuromusa)", "Sva1033", "italic(Desulfuromonas)", "italic(Arcobacteraceae)")))

## add day 0 data for all treatments from DIC day 0 mean
ta5s_MSFig03_d0 <- rbind(ta5s_MSFig03,
                         ta5s_MSFig03 %>% filter(Substrate == "DIC" & Timepoint == 0) %>% mutate(Substrate = "acetate"),
                         ta5s_MSFig03 %>% filter(Substrate == "DIC" & Timepoint == 0) %>% mutate(Substrate = "birnessite.acetate")) %>% 
  filter(Substrate != "DIC")

## calculate mean per treatment
calc_sub_ta5s_MSFig03 <- ta5s_MSFig03_d0 %>% 
  select(-Replicate) %>% 
  group_by(Substrate, Timepoint, Genus_mod, NucleicAcid) %>% 
  summarise(mean = mean(Abundance), .groups = "keep")

## combine
ta5s_MSFig03_wmean <- full_join(ta5s_MSFig03_d0, calc_sub_ta5s_MSFig03, relationship = "many-to-many") %>% 
  mutate(Substrate = factor(Substrate, levels = c("birnessite.acetate", "acetate"),
                            labels = c("Birnessite + acetate", "Acetate")))


## plot ####
man_col <- c("#ae017e", "#31688EFF")
plot_Fig03 <- ggplot(ta5s_MSFig03_wmean, aes(x = Timepoint, color = Substrate, linetype = Substrate)) +
  geom_line(aes(y = mean*100)) +
  geom_point(aes(y = Abundance*100, shape = Replicate), size = 1.5) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,2,6,10,14,20), limits = c(-0.1,21)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.2,40.2)) +
  scale_color_manual(values = man_col) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_shape_discrete(na.translate = F) +
  facet_grid(NucleicAcid~Genus_mod, labeller = label_parsed) +
  labs(x = "Incubation time (days)", y = "Relative abundance (%)") +
  guides(shape = guide_legend(title.position = "left", order = 2), 
         color = guide_legend(order = 1), linetype = guide_legend(order = 1)) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 9, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 9.5),
        panel.spacing.y = unit(1, "lines"),
        legend.position = "bottom", legend.box.margin = margin(t = -10),
        legend.key.width = unit(0.8, "cm"))
plot_Fig03
ggsave("plots/MS_Fig03_Inc_tax_new.jpg", plot_Fig03, width = 180, height = 105, units = "mm")
ggsave("plots/MS_Fig03_Inc_tax_new.pdf", plot_Fig03, width = 180, height = 105, units = "mm")
ggsave("plots/MS_Fig03_Inc_tax_new.emf", plot_Fig03, width = 180, height = 105, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


# line plot supl MS sulfate reducers ####
## lineplot, only Mn + acetate, acetate treatment DNA and RNA level
# uncl. Desulfocapsaceae, Desulforhopalus, SEEP-SRB4, Desulfopila, uncl. Desulfobacteraceae, Desulfoconvexum, Desulfofrigus

## subset data ####
### test different subsets which ones make most sense
ta5n_MSFigX_SR <- ta5n %>% 
  filter(Substrate %in% c("acetate", "birnessite.acetate", "DIC") &
           (Genus %in% c("Desulforhopalus", "SEEP-SRB4", "Desulfopila",
                         "Desulfofrigus", 
                         "Desulfoconvexum") |
              (Family == "Desulfobacteraceae" & is.na(Genus)) |
              (Family == "Desulfocapsaceae" & is.na(Genus))) &
           Timepoint != 30)

ta5s_MSFigX_SR <- ta5n_MSFigX_SR %>% 
  group_by(Sample, Family, Genus, Substrate, Replicate, Timepoint, NucleicAcid) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(Fam_Genus = factor(paste0(Family, "_", Genus),
                            levels = c("Desulfocapsaceae_NA", "Desulfocapsaceae_Desulforhopalus", 
                                       "Desulfocapsaceae_SEEP-SRB4", "Desulfocapsaceae_Desulfopila", 
                                       "Desulfolunaceae_Desulfofrigus", 
                                       "Desulfobacteraceae_NA",  "Desulfobacteraceae_Desulfoconvexum"),
                            labels = c("unclassified\nDesulfocapsaceae", "Desulforhopalus\n(Desulfocapsaceae)", 
                                       "SEEP-SRB4\n(Desulfocapsaceae)", "Desulfopila\n(Desulfocapsaceae)", 
                                       "Desulfofrigus\n(Desulfolunaceae)", 
                                       "unclassified\nDesulfobacteraceae",  "Desulfoconvexum\n(Desulfobacteraceae)")))

## add day 0 data for all treatments from DIC day 0 mean
ta5s_MSFigX_SR_d0 <- rbind(ta5s_MSFigX_SR,
                           ta5s_MSFigX_SR %>% filter(Substrate == "DIC" & Timepoint == 0) %>% mutate(Substrate = "acetate"),
                           ta5s_MSFigX_SR %>% filter(Substrate == "DIC" & Timepoint == 0) %>% mutate(Substrate = "birnessite.acetate")) %>% 
  filter(Substrate != "DIC")


## calculate mean per treatment
calc_sub_ta5s_MSFigX_SR <- ta5s_MSFigX_SR_d0 %>% 
  select(-Replicate) %>% 
  group_by(Substrate, Timepoint, Fam_Genus, NucleicAcid) %>% 
  summarise(mean = mean(Abundance), .groups = "keep")

## combine
ta5s_MSFigX_SR_wmean <- full_join(ta5s_MSFigX_SR_d0, calc_sub_ta5s_MSFigX_SR, relationship = "many-to-many") %>% 
  mutate(Substrate = factor(Substrate, levels = c("birnessite.acetate", "acetate"),
                            labels = c("Birnessite + acetate", "Acetate")))


## plot ####
man_col <- c("#ae017e", "#31688EFF")
plot_MSFigX_SR <- ggplot(ta5s_MSFigX_SR_wmean, aes(x = Timepoint, color = Substrate, linetype = Substrate)) +
  geom_line(aes(y = mean*100)) +
  geom_point(aes(y = Abundance*100, shape = Replicate), size = 1.5) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), limits = c(-0.2,21)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.2,19)) +
  scale_color_manual(values = man_col) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_shape_discrete(na.translate = F) +
  facet_grid(NucleicAcid~Fam_Genus) +
  labs(x = "Incubation time (days)", y = "Relative abundance (%)") +
  guides(shape = guide_legend(order = 2), 
         color = guide_legend(order = 1), linetype = guide_legend(order = 1)) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 8, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0.5, "lines"),
        strip.background = element_blank(), strip.text = element_text(size = 6.9),
        legend.position = "bottom", legend.box.margin = margin(t = -8),
        legend.key.width = unit(0.8, "cm"))
plot_MSFigX_SR
ggsave("plots/MS_suplFigX_Inc_tax_SR.jpg", plot_MSFigX_SR, width = 180, height = 80, units = "mm")
ggsave("plots/MS_suplFigX_Inc_tax_SR.pdf", plot_MSFigX_SR, width = 180, height = 80, units = "mm")
ggsave("plots/MS_suplFigX_Inc_tax_SR.emf", plot_MSFigX_SR, width = 180, height = 80, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
