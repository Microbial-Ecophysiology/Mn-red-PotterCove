## Potter Cove in situ for manganese reduction manuscript
 # plot dissolved Mn, Fe2+ and sulfide in situ pore water concentrations Station 01
 # cores are in duplicates

# load packages ####
library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(devEMF)
library(patchwork)

# import data ####
setwd("WORKINGDIRECTORY")

df <- read_tsv("INPUTDATA.txt", name_repair = make.names)
df$Core <- as.character(df$Core)

## subset to only STA01 and parameters to plot
df1 <- subset(df, df$Station == "STA01") %>% 
  select("Core_ID", "Station", "Core", "Depth_cm", "Fe2._uM", "Mn_uM", "SO42._mM")


# plot ####
## Mn
plot_Mn <- ggplot(df1, aes(x = Mn_uM, y = Depth_cm, shape = Core_ID, group = Core_ID)) +
  geom_path() +
  geom_point(size = 1.5) +
  scale_y_reverse(name = "Core depth (cm)", expand = c(0,0.1), limits = c(33, 0)) +
  scale_x_continuous(position = "top", name = "Dissolved Mn (?M)", expand = c(0,0.1), limits = c(0,52)) +
  labs(shape = "Replicate core") +
  theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12, color = "black"), 
        axis.title =  element_text(size = 11),
        axis.text = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank(),
        plot.margin = unit(c(0.1,0.4,0.2,0.1), "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = unit(c(-0.3,0,0,0), "cm"))

## Fe2+
plot_Fe <- ggplot(df1, aes(x = Fe2._uM, y = Depth_cm, shape = Core_ID, group = Core_ID)) +
  geom_path() +
  geom_point(size = 1.5) +
  scale_y_reverse(name = "Core depth (cm)", expand = c(0,0.1), limits = c(33, 0)) +
  scale_x_continuous(position = "top", name = expression(Fe^{"2+"}~"(?M)"), expand = c(0,0.1), limits = c(0,110)) +
  labs(shape = "Replicate core") +
  theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12, color = "black"), 
        axis.title =  element_text(size = 11),
        axis.text = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank(),
        plot.margin = unit(c(0.1,0.4,0.2,0.1), "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = unit(c(-0.3,0,0,0), "cm"))

## SO4
plot_SO4 <- ggplot(df1, aes(x = SO42._mM, y = Depth_cm, shape = Core_ID, group = Core_ID)) +
  geom_path() +
  geom_point(size = 1.5) +
  scale_y_reverse(name = "Core depth (cm)", expand = c(0,0.1), limits = c(33, 0)) +
  scale_x_continuous(position = "top", name = expression(SO["4"]^{"2-"}~"(mM)"), expand = c(0,0.1), limits = c(0,30)) +
  labs(shape = "Replicate core") +
  theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 12, color = "black"), 
        axis.title =  element_text(size = 11),
        axis.text = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank(),
        plot.margin = unit(c(0.1,0.4,0.2,0.1), "cm"),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = unit(c(-0.3,0,0,0), "cm"))

## combine plots
plot_comb <- plot_Mn + plot_Fe + plot_SO4 + plot_annotation(tag_levels = "A") + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom",
        text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank())

ggsave("plots/MS_Fig01_insitu_conc.pdf", plot_comb, width = 180, height = 90, units = "mm")
ggsave("plots/MS_Fig01_insitu_conc.jpg", plot_comb, width = 180, height = 90, units = "mm", dpi = 300)
ggsave("plots/MS_Fig01_insitu_conc.emf", plot_comb, width = 180, height = 90, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})
