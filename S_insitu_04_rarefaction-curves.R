# This script imports taxmap objects created in script metacoder_import_data.R produced through dada2 pipeline

# explore data and create rarefaction curves 

# server part ####
## packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(iNEXT)
library(cowplot)

## import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

load(paste0("taxmap_", prj, ".RData"))

mdata <- read_tsv("MDATAFILE.txt")

sample.names <- mdata$Name.Mapping.File

asv_counts_clean <- as.data.frame(obj$data$asv_table[, sample.names])

# rarefaction curves ####
# calculate variety of statistics with iNEXT
stat_all <- iNEXT(asv_counts_clean, q = c(0,2), datatype = "abundance", knots = 500,
                  endpoint = max(colSums(asv_counts_clean)))
save(stat_all, file = paste0("iNext_", prj, ".RData"))


# PC part ####
## load packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(iNEXT)
library(ggrepel)
library(cowplot)

## import data ####
setwd("WORKINGDIRECTORY")

## project name to name files
prj <- "PROJECT"

## taxmap object
load(paste0("taxmap_", prj, ".RData"))

## metadata
mdata <- read_tsv("MDATAFILE.txt")

sample.names <- mdata$Name.Mapping.File

# basic data overview - counts per sample ####
hist(colSums(obj$data$asv_table[, sample.names]))

## total reads and ASVs per sample
mdata <- mdata %>%
  left_join((
    colSums(obj$data$asv_table[, sample.names]) %>%
      as.data.frame() %>%
      rename("totReads" = ".") %>%
      rownames_to_column(var = "Sequencing.ID.Bac")
  )) %>%
  left_join((
    apply(obj$data$asv_table[, sample.names] , 2, function(x) sum(x > 0)) %>%
      as.data.frame() %>%
      rename("totASVs" = ".") %>%
      rownames_to_column(var = "Sequencing.ID.Bac")
  ))

write.table(mdata, "MDATAFILE_totReads-totOTUS.txt", sep = "\t", 
            quote = F, row.names = FALSE)


load(paste0("iNext_", prj, ".RData"))

stat_all$DataInfo
stat_all$AsyEst


# export iNEXT data in separate data frame to plot manually
df <- fortify(stat_all, type = 1)

# for plotting make sure factors are defined as factors
cols_factor <- (c("Method", "Core", "Station", "Section", "replicate..for.graphs.", "Seq.lib.no.", "Depth..cm."))
# add mdata to make nicer plot
dfm <- left_join(df, mdata, by = c("Assemblage" = "Sequencing.ID.Bac")) %>% 
  filter(Core == "STA01.02") %>% 
  mutate(Order.q = factor(Order.q, levels = c("0", "2"), labels = c("Species richness", "Linearized\nsimpson diversity")),
         across(all_of(cols_factor), as.factor),
         label = gsub("Bacteria.", "", Assemblage))

## supl fig rarefaction curve in situ ####
figx <- ggplot(dfm, aes(x = x, y = y, color = Depth..cm.)) +
  geom_line(data = subset(dfm, dfm$Method == "Rarefaction"), 
            aes(group = Assemblage), linetype = "solid") +
  geom_point(data = subset(dfm, dfm$Method == "Observed"), 
             aes(), size = 1.5) +
  scale_x_continuous(expand = c(0.018,0), limits = c(0,NA)) +
  scale_y_continuous(expand = c(0.025,0), limits = c(0,NA)) +
  labs(x = "Number of individuals", y = "Species diversity", linetype = "Method",
       title = NULL) +
  guides(color = guide_legend(title = "Core depth (cm)", nrow = 2)) +
  facet_wrap(.~Order.q, scales = "free_y") +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        strip.text = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.5),
        strip.background = element_rect(colour = "black", fill = "white"),
        legend.position = "bottom")

ggsave("plots/MS_supl_FigX_rarefaction_insitu_STA01.02.jpg", figx, width = 185, height = 100, units = "mm")
ggsave("plots/MS_supl_FigX_rarefaction_insitu_STA01.02.pdf", figx, width = 185, height = 100, units = "mm")
ggsave("plots/MS_supl_FigX_rarefaction_insitu_STA01.02.emf", figx, width = 185, height = 100, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
