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
mdata <- read_tsv("MDATAFILE.txt")


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
# 5 genera, 12 families, 9 orders, 10 classes, 9 phyla

## taxa above 2% in at least one samples
# ta2 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 2, Abundance = "Abundance")
# 16 genera, 28 families, 23 orders, 19 classes, 15 phyla

# use 5% threshold
# rm(ta2, ta3, ta4)

# sort and make nice names for taxa ####
## other columns needed as factor for plotting
cols_factor <- c("Incubation", "Replicate", "Substrate", "NucleicAcid", "Year")

## decide on rank level to plot: Genus and sorted by Class
ta5_class <- c(unique(ta5$ASV_table_taxa_abv_5_perc$Class_mod)) # 19 taxa to sort by
ta5_genus <- c(unique(ta5$ASV_table_taxa_abv_5_perc$Genus_mod)) # 40 taxa to plot
### reordering phyla and class entries using function from https://github.com/mrdwab/SOfun
source("H:/GitHub/PhD_other-code/moveMe.R")
ta5_genus_reorder <- moveMe(ta5_genus, "other_c_Alphaproteobacteria_<5% after other_p_Planctomycetota_<5%")
ta5_genus_renamed <- ta5_genus_reorder %>% 
  gsub("other_p_Actinobacteriota_<5%", "p_Actinobacteriota", .) %>% 
  gsub("other_o_Bacteroidales_<5%", "o_Bacteroidales", .) %>% 
  gsub("other_f_Flavobacteriaceae_<5%", "f_Flavobacteriaceae", .) %>% 
  gsub("other_f_Arcobacteraceae_<5%", "f_Arcobacteraceae", .) %>% 
  gsub("other_f_Desulfobacteraceae_<5%", "f_Desulfobacteraceae", .) %>% 
  gsub("other_f_Desulfobulbaceae_<5%", "f_Desulfobulbaceae", .) %>% 
  gsub("other_f_Sva1033_<5%", "f_Sva1033", .) %>% 
  gsub("other_c_Clostridia_<5%", "c_Clostridia", .) %>% 
  gsub("other_o_MSBL9_<5%", "o_MSBL9", .) %>% 
  gsub("other_f_Pirellulaceae_<5%", "f_Pirellulaceae", .) %>% 
  gsub("other_c_Alphaproteobacteria_<5%", "c_Alphaproteobacteria", .) %>% 
  gsub("other_f_Unknown Family_<5%", "Unknown_Family_Gammapr.", .) %>% 
  gsub("other_p_Spirochaetota_<5%", "p_Spirochaetota", .) %>% 
  gsub("other_p_Verrucomicrobiota_<5%", "p_Verrucomicrobiota", .)


ta5n <- select(ta5$ASV_table_taxa_abv_5_perc, OTU, Sample, Class_mod, Genus_mod, Abundance, Incubation, Substrate, Replicate, Timepoint, NucleicAcid, Year, Genus, Family, Class, Order, Phylum) %>% 
  mutate(Class_mod=factor(Class_mod, levels = ta5_class),
         Genus_mod=factor(Genus_mod, levels = ta5_genus_reorder, labels = ta5_genus_renamed),
         ASV.taxa = factor(paste0(OTU, "\n", Genus)),  # new column with ASV ID and genus name
         across(all_of(cols_factor), as.factor)) %>%   # reorder families so families of same phylum are together
  group_by(Genus_mod)


# sum ASV abundance by taxon ####
## need to decide which rank to plot and to group by which other rank
 # here plot genus, group by class
ta5s <- ta5n %>% 
  group_by(Sample, Class_mod, Genus_mod, Incubation, Substrate, Replicate, Timepoint, NucleicAcid, Year) %>% 
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

plot_dic <- plot_basic %+% subset(ta5s_R, ta5s_R$Substrate == "DIC" & ta5s_R$Timepoint != 0) +
  labs(title = "DIC (A+B: 2021; C: 2023)")

plot_birnessite <- plot_basic %+% subset(ta5s_R, ta5s_R$Substrate == "birnessite" & ta5s_R$Timepoint != 0) +
  labs(title = "Birnessite + DIC")

plot_day0 <- plot_basic %+% subset(ta5s_R, ta5s_R$Timepoint == 0) +
  labs(title = "Initial slurry") +
  facet_wrap(~Year, drop = T) +
  theme(strip.background = element_rect(color = "black"),
        axis.title.x = element_blank())


## MS supl Fig taxonomy RNA ####
MS_supl_plot_RNA <- plot_birnessite.acetate + plot_acetate + plot_birnessite + plot_dic + plot_day0 + plot_spacer() +  
  plot_layout(guides = "collect", nrow = 2,) + 
  plot_annotation(tag_levels = "A", title = "Bacterial community RNA level",
                  theme = theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))) &
  theme(legend.position = "bottom")


ggsave("plots/MS_supl_FigSX_bac_barplot_RNA.jpg", plot = MS_supl_plot_RNA, width = 180, height = 190, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_RNA.pdf", plot = MS_supl_plot_RNA, width = 180, height = 190, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_RNA.emf", plot = MS_supl_plot_RNA, width = 180, height = 190, units = "mm",
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

plot_D_dic <- plot_D_basic %+% subset(ta5s_D, ta5s_D$Substrate == "DIC" & ta5s_D$Timepoint != 0) +
  labs(title = "DIC (A+B: 2021; C: 2023)")

plot_D_birnessite <- plot_D_basic %+% subset(ta5s_D, ta5s_D$Substrate == "birnessite" & ta5s_D$Timepoint != 0) +
  labs(title = "Birnessite + DIC")

plot_D_day0 <- plot_D_basic %+% subset(ta5s_D, ta5s_D$Timepoint == 0) +
  labs(title = "Initial slurry") +
  facet_wrap(~Year, drop = T) +
  theme(strip.background = element_rect(color = "black"),
        axis.title.x = element_blank())


## MS supl Fig taxonomy DNA ####
MS_supl_plot_DNA <- plot_D_birnessite.acetate + plot_D_acetate + plot_D_birnessite + plot_D_dic + plot_D_day0 + plot_spacer() + 
  plot_layout(guides = "collect", nrow = 2) + 
  plot_annotation(tag_levels = "A", title = "Bacterial community DNA level",
                  theme = theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))) &
  theme(legend.position = "bottom")

ggsave("plots/MS_supl_FigSX_bac_barplot_DNA.jpg", plot = MS_supl_plot_DNA, width = 180, height = 190, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_DNA.pdf", plot = MS_supl_plot_DNA, width = 180, height = 190, units = "mm")
ggsave("plots/MS_supl_FigSX_bac_barplot_DNA.emf", plot = MS_supl_plot_DNA, width = 180, height = 190, units = "mm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})



# plot MS Fig 3 (is really Fig 4 in the end) ####
## lineplot, only Mn + acetate, acetate treatment
 # Desulfuromusa, Desulfuromonas, Sva1033, uncl. Arcobacteraceae

## subset data ####
ta5n_MSFig03 <- ta5n %>% 
  filter(Substrate %in% c("acetate", "birnessite.acetate", "DIC") &
           (Genus %in% c("Desulfuromonas", "Desulfuromusa") |
              (Family == "Arcobacteraceae" & is.na(Genus)) |
              Family == "Sva1033"))

ta5s_MSFig03 <- ta5n_MSFig03 %>% 
  group_by(Sample, Genus_mod, Substrate, Replicate, Timepoint, NucleicAcid, Year) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(Genus_mod = factor(Genus_mod, levels = c("Desulfuromusa", "f_Sva1033", "Desulfuromonas", "f_Arcobacteraceae"),
                            labels = c("italic(Desulfuromusa)", "Sva1033", "italic(Desulfuromonas)", "italic(Arcobacteraceae)")))

## add day 0 data for all treatments from DIC day 0 mean for 2021 treatments and birnessite day 0 to DIC 2023 treatment
ta5s_MSFig03_d0 <- rbind(ta5s_MSFig03,
                         ta5s_MSFig03 %>% filter(Substrate == "DIC" & Timepoint == 0 & Year == "2021") %>% mutate(Substrate = "acetate"),
                         ta5s_MSFig03 %>% filter(Substrate == "DIC" & Timepoint == 0 & Year == "2021") %>% mutate(Substrate = "birnessite.acetate"),
                         ta5s_MSFig03 %>% filter(Substrate == "birnessite" & Timepoint == 0 & Year == "2023") %>% mutate(Substrate = "DIC"))

## calculate mean per treatment
calc_sub_ta5s_MSFig03 <- ta5s_MSFig03_d0 %>% 
  select(-Replicate, -Year) %>% 
  group_by(Substrate, Timepoint, Genus_mod, NucleicAcid) %>% 
  summarise(mean = mean(Abundance), .groups = "keep")

## combine
ta5s_MSFig03_wmean <- full_join(ta5s_MSFig03_d0, calc_sub_ta5s_MSFig03, relationship = "many-to-many") %>% 
  mutate(Substrate = factor(Substrate, levels = c("birnessite.acetate", "birnessite", "acetate",  "DIC"),
                            labels = c("Birnessite + acetate", "Birnessite + DIC", "Acetate", "DIC")))


## plot different scales ####
man_col <- c("#ae017e", "#440154FF", "#31688EFF", "#8FD744FF")

# base plot
fig03_base <- ggplot(ta5s_MSFig03_wmean, aes(x = Timepoint, color = Substrate, linetype = Substrate)) +
  geom_line(aes(y = mean*100)) +
  geom_point(aes(y = Abundance*100, shape = Replicate), size = 1.5) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,5,10,15,20), limits = c(-0.1,21)) +
  scale_color_manual(values = man_col) +
  scale_linetype_manual(values = c("solid", "1343", "22", "13")) +
  scale_shape_discrete(na.translate = F) +
  facet_grid(NucleicAcid~Genus_mod, labeller = label_parsed) +
  labs(x = "Incubation time (days)", y = "Relative abundance (%)") +
  guides(shape = guide_legend(title.position = "left", order = 2), 
         color = guide_legend(order = 1, title = "Treatment"), linetype = guide_legend(order = 1, title = "Treatment")) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 9, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.title.y = element_text(margin = margin(l = 3)),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 9.5),
        panel.spacing.y = unit(1, "lines"),
        legend.position = "bottom", legend.box.margin = margin(t = -10),
        legend.key.width = unit(0.8, "cm"))

# subset plots for different scales
fig3_desulfu_d <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "italic(Desulfuromusa)" &
                                          ta5s_MSFig03_wmean$NucleicAcid == "DNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,17), breaks = c(0,5,10,15)) +
  theme(strip.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig3_sva_d <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "Sva1033" &
                                      ta5s_MSFig03_wmean$NucleicAcid == "DNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,17), breaks = c(0,5,10,15)) +
  theme(strip.text.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig3_desulfo_d <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "italic(Desulfuromonas)" &
                                          ta5s_MSFig03_wmean$NucleicAcid == "DNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.38,63)) +
  theme(strip.text.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

fig3_arco_d <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "italic(Arcobacteraceae)" &
                                       ta5s_MSFig03_wmean$NucleicAcid == "DNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.38,63)) +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


fig3_desulfu_r <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "italic(Desulfuromusa)" &
                                          ta5s_MSFig03_wmean$NucleicAcid == "RNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,17), breaks = c(0,5,10,15)) +
  theme(strip.text.y = element_blank(), strip.text.x = element_blank())

fig3_sva_r <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "Sva1033" &
                                      ta5s_MSFig03_wmean$NucleicAcid == "RNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,17), breaks = c(0,5,10,15)) +
  theme(strip.text.y = element_blank(), axis.title.y = element_blank(), strip.text.x = element_blank())

fig3_desulfo_r <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "italic(Desulfuromonas)" &
                                          ta5s_MSFig03_wmean$NucleicAcid == "RNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.38,63)) +
  theme(strip.text.y = element_blank(), axis.title.y = element_blank(), strip.text.x = element_blank())

fig3_arco_r <- fig03_base %+% subset(ta5s_MSFig03_wmean, ta5s_MSFig03_wmean$Genus_mod == "italic(Arcobacteraceae)" &
                                       ta5s_MSFig03_wmean$NucleicAcid == "RNA") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.38,63)) +
  theme(axis.title.y = element_blank(), strip.text.x = element_blank())



# combine plots
fig3_axis <- fig3_desulfu_d + fig3_sva_d + fig3_desulfo_d + fig3_arco_d + 
  fig3_desulfu_r + fig3_sva_r + fig3_desulfo_r + fig3_arco_r +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2, guides = "collect", axis_titles = "collect_x") &
  theme(legend.position = "bottom", legend.box.margin = margin(t = -5),
        legend.key.width = unit(0.8, "cm"),
        legend.box = "vertical", 
        legend.margin = margin(t = -0.5))
fig3_axis

ggsave("plots/MS_Fig03_Inc_tax_axis.jpg", fig3_axis, width = 180, height = 108, units = "mm", dpi = 600)
ggsave("plots/MS_Fig03_Inc_tax_axis.pdf", fig3_axis, width = 180, height = 108, units = "mm")
ggsave("plots/MS_Fig03_Inc_tax_axis.emf", fig3_axis, width = 180, height = 108, units = "mm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


# line plot like Fig03 distinguished between years for DIC control ####
ta5s_MSFig03 <- ta5n_MSFig03 %>% 
  group_by(Sample, Genus_mod, Substrate, Replicate, Timepoint, NucleicAcid, Year) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(Genus_mod = factor(Genus_mod, levels = c("Desulfuromusa", "f_Sva1033", "Desulfuromonas", "f_Arcobacteraceae"),
                            labels = c("italic(Desulfuromusa)", "Sva1033", "italic(Desulfuromonas)", "italic(Arcobacteraceae)")))

ta5s_MSFigX_sepDIC <- ta5s_MSFig03 %>% 
  mutate(Treatment_sep = case_when(Substrate == "DIC" & Year == "2021" ~ "DIC_2021",
                                   Substrate == "DIC" & Year == "2023" ~ "DIC_2023",
                                   Substrate != "DIC" ~ Substrate))

## add day 0 data for all treatments from DIC day 0 mean for 2021 treatments and birnessite day 0 to DIC 2023 treatment
ta5s_MSFigX_sepDIC_d0 <- rbind(ta5s_MSFigX_sepDIC,
                               ta5s_MSFigX_sepDIC %>% filter(Treatment_sep == "DIC_2021" & Timepoint == 0) %>% mutate(Treatment_sep = "acetate"),
                               ta5s_MSFigX_sepDIC %>% filter(Treatment_sep == "DIC_2021" & Timepoint == 0) %>% mutate(Treatment_sep = "birnessite.acetate"),
                               ta5s_MSFigX_sepDIC %>% filter(Treatment_sep == "birnessite" & Timepoint == 0 & Year == "2023") %>% mutate(Treatment_sep = "DIC_2023"))

## calculate mean per treatment
calc_sub_ta5s_MSFigX_sepDIC <- ta5s_MSFigX_sepDIC_d0 %>% 
  select(-Replicate, -Year) %>% 
  group_by(Substrate, Timepoint, Genus_mod, NucleicAcid, Treatment_sep) %>% 
  summarise(mean = mean(Abundance), .groups = "keep")

## combine
ta5s_MSFigX_sepDIC_wmean <- full_join(ta5s_MSFigX_sepDIC_d0, calc_sub_ta5s_MSFigX_sepDIC, relationship = "many-to-many") %>% 
  mutate(Treatment_sep = factor(Treatment_sep, levels = c("birnessite.acetate", "birnessite", "acetate",  "DIC_2021", "DIC_2023"),
                                labels = c("Birnessite + acetate", "Birnessite + DIC", "Acetate", "DIC_2021", "DIC_2023")))


## plot ####
man_col_sepDIC <- c("#ae017e", "#440154FF", "#31688EFF", "#8FD744FF", "#406f2aff")
plot_sepDIC <- ggplot(ta5s_MSFigX_sepDIC_wmean, aes(x = Timepoint, color = Treatment_sep, linetype = Treatment_sep)) +
  geom_line(aes(y = mean*100)) +
  geom_point(aes(y = Abundance*100, shape = Replicate), size = 1.5) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,2,6,10,14,20), limits = c(-0.1,21)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.2,63)) +
  scale_color_manual(values = man_col_sepDIC) +
  scale_linetype_manual(values = c("solid", "1343", "22", "13", "13")) +
  scale_shape_discrete(na.translate = F) +
  facet_grid(NucleicAcid~Genus_mod, labeller = label_parsed) +
  labs(x = "Incubation time (days)", y = "Relative abundance (%)") +
  guides(shape = guide_legend(title.position = "left", order = 2), 
         color = guide_legend(order = 1, title = "Treatment"), linetype = guide_legend(order = 1, title = "Treatment")) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 9, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 9.5),
        panel.spacing.y = unit(1, "lines"),
        legend.position = "bottom", legend.box.margin = margin(t = -5),
        legend.key.width = unit(0.8, "cm"),
        legend.box = "vertical", legend.spacing.y = unit(0, "cm"), 
        legend.margin = margin(t = -0.5, b = -0.5))
plot_sepDIC


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
              (Family == "Desulfocapsaceae" & is.na(Genus))))

ta5s_MSFigX_SR <- ta5n_MSFigX_SR %>% 
  group_by(Sample, Family, Genus, Substrate, Replicate, Timepoint, NucleicAcid, Year) %>% 
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
                           ta5s_MSFigX_SR %>% filter(Treatment_sep == "DIC_2021" & Timepoint == 0) %>% mutate(Treatment_sep = "acetate"),
                           ta5s_MSFigX_SR %>% filter(Treatment_sep == "DIC_2021" & Timepoint == 0) %>% mutate(Treatment_sep = "birnessite.acetate"),
                           ta5s_MSFigX_SR %>% filter(Treatment_sep == "birnessite" & Timepoint == 0 & Year == "2023") %>% mutate(Treatment_sep = "DIC_2023"))


## calculate mean per treatment
calc_sub_ta5s_MSFigX_SR <- ta5s_MSFigX_SR_d0 %>% 
  select(-Replicate) %>% 
  group_by(Substrate, Timepoint, Fam_Genus, NucleicAcid) %>% 
  summarise(mean = mean(Abundance), .groups = "keep")

## combine
ta5s_MSFigX_SR_wmean <- full_join(ta5s_MSFigX_SR_d0, calc_sub_ta5s_MSFigX_SR, relationship = "many-to-many") %>% 
  mutate(Treatment_sep = factor(Treatment_sep, levels = c("birnessite.acetate", "birnessite", "acetate",  "DIC_2021", "DIC_2023"),
                                labels = c("Birnessite + acetate", "Birnessite + DIC", "Acetate", "DIC_2021", "DIC_2023")))


## plot ####
man_col_sepDIC <- c("#ae017e", "#440154FF", "#31688EFF", "#8FD744FF", "#406f2aff")
plot_MSFigX_SR <- ggplot(ta5s_MSFigX_SR_wmean, aes(x = Timepoint, color = Treatment_sep, linetype = Treatment_sep)) +
  geom_line(aes(y = mean*100)) +
  geom_point(aes(y = Abundance*100, shape = Replicate), size = 1.1) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,10,20), limits = c(-0.2,21)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.2,22)) +
  scale_color_manual(values = man_col_sepDIC) +
  scale_linetype_manual(values = c("solid", "1343", "22", "13", "13")) +
  scale_shape_discrete(na.translate = F) +
  facet_grid(NucleicAcid~Fam_Genus) +
  labs(x = "Incubation time (days)", y = "Relative abundance (%)") +
  guides(shape = guide_legend(order = 2), 
         color = guide_legend(order = 1, title = "Treatment"), linetype = guide_legend(order = 1, title = "Treatment")) +
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
ggsave("plots/MS_suplFigX_Inc_tax_SR.jpg", plot_MSFigX_SR, width = 180, height = 82, units = "mm")
ggsave("plots/MS_suplFigX_Inc_tax_SR.pdf", plot_MSFigX_SR, width = 180, height = 82, units = "mm")
ggsave("plots/MS_suplFigX_Inc_tax_SR.emf", plot_MSFigX_SR, width = 180, height = 82, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
