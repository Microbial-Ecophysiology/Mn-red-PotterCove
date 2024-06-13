## Plotting geochemical data from manganese reduction incubation experiment 

# load packages ####
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(devEMF)
library(multcomp)

# import data ####
setwd("WORKINGDIRECTORY")

data <- read_tsv("INPUTDATA.txt")


# MS plot ####
## Fig 2 - Mn ####
# only take Mn values and calc mean together for both DIC treatments
sub_geo_mod <- data %>% 
  dplyr::select(Treatment, Replicate, Day, Mn.uM) %>% 
  subset(Treatment %in% c("Bir.Ac", "Bir", "Bir.S", "Bir.TS", 
                          "Ac", "C", "C_2") & Day <= 20) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Bir.Ac", "Bir.S", "Bir.TS", "Bir", 
                                                   "Ac", "C", "C_2"),
                            labels = c("Birnessite + acetate",  "Birnessite + sulfur + DIC",
                                       "Birnessite + thiosulfate + DIC",
                                       "Birnessite + DIC", "Acetate", "DIC", "DIC"))) %>% 
  drop_na(Mn.uM)

## calculate mean per treatment
calc_sub_Mn <- sub_geo_mod %>% 
  dplyr::select(-Replicate) %>% 
  group_by(Treatment, Day) %>% 
  summarise(mean = mean(Mn.uM), .groups = "keep")

## combine
Mn_MS <- full_join(sub_geo_mod, calc_sub_Mn)

## plot
man_col <- c("#ae017e", "#bcbddc", "#8c6bb1", "#440154FF","#31688EFF", "#8FD744FF")
fig2 <- ggplot(Mn_MS, aes(x = Day, color = Treatment, linetype = Treatment)) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = Mn.uM, shape = Replicate), size = 1.5) +
  scale_colour_manual(values = man_col) +
  scale_linetype_manual(values = c("solid", "42", "44",  "1343", "22", "13")) +
  scale_x_continuous(expand = c(0.01,0.2), breaks = c(0,5,10,15,20)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 630)) +
  labs(x = "Incubation time (days)", y = "Dissolved Mn (µM)") +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = "bottom", legend.box.margin = margin(l = -0.5, t = -0.3, unit = "cm"),
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"))
fig2
ggsave("plots/MS_Fig02_Inc_conc.jpg", fig2, width = 85, height = 100, units = "mm")
ggsave("plots/MS_Fig02_Inc_conc.pdf", fig2, width = 85, height = 100, units = "mm")
ggsave("plots/MS_Fig02_Inc_conc.emf", fig2, width = 85, height = 100, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})



## supl Fig SX - incubation geochemistry ####
# take Mn, sulfate, iron, sulfide values and calc mean separate for both DIC treatments
sub_geo_supl <- data %>% 
  mutate(SO_uM = SO_mM*1000) %>% 
  dplyr::select(Treatment, Replicate, Day, Mn.uM, Fe.fer.uM, SO_uM, Sulfide_uM) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Bir.Ac", "Bir.S", "Bir.TS", "Bir", 
                                                  "Ac", "C", "C_2"),
                            labels = c("Birnessite + acetate",  "Birnessite + sulfur + DIC",
                                       "Birnessite + thiosulfate + DIC",
                                       "Birnessite + DIC", "Acetate", "DIC", "DIC_2023")))

## calculate mean per treatment
calc_sub_supl <- sub_geo_supl %>% 
  dplyr::select(-Replicate) %>% 
  group_by(Treatment, Day) %>% 
  summarise(across(everything(), ~ mean(., na.rm = T)), .groups = "keep")

## long format
l_calc_sub_supl <- melt(calc_sub_supl, id.vars = c("Treatment", "Day"), value.name = "mean", variable.name = "geo", na.rm = T,
                 measure.vars = c("SO_uM", "Mn.uM", "Fe.fer.uM", "Sulfide_uM"))
l_sub_geo_supl <- melt(sub_geo_supl, id.vars = c("Treatment", "Replicate", "Day"), value.name = "conc", variable.name = "geo", na.rm = T,
                measure.vars = c("SO_uM", "Mn.uM", "Fe.fer.uM", "Sulfide_uM"))

geo_supl_comb <- full_join(l_sub_geo_supl, l_calc_sub_supl) %>% 
  mutate(geo = factor(geo, levels = c("Mn.uM", "Fe.fer.uM", "Sulfide_uM", "SO_uM"),
                      labels = c("Manganese", "Ferrous iron", "Sulfide", "Sulfate")))


## plot
man_col_supl <- c("#ae017e", "#bcbddc", "#8c6bb1", "#440154FF","#31688EFF", "#8FD744FF", "#406f2aff")
fig_supl <- ggplot(geo_supl_comb, aes(x = Day, color = Treatment, linetype = Treatment)) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = conc, shape = Replicate), size = 1.5) +
  scale_colour_manual(values = man_col_supl) +
  scale_linetype_manual(values = c("solid", "42", "44",  "1343", "22", "13", "13")) +
  scale_x_continuous(expand = c(0.01,0.2), breaks = seq(0,20,5)) +
  scale_y_continuous(expand = c(0.04,0), limits = c(0, NA)) +
  labs(x = "Incubation time (days)", y = "Concentration (µM)") +
  guides(color = guide_legend(ncol = 2, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  facet_wrap(~geo, scales = "free_y", nrow = 2) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = NA), strip.text = element_text(size = 10),
        legend.position = "bottom", legend.box.margin = margin(l = -0.5, t = -0.3, unit = "cm"),
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"))
fig_supl

ggsave("plots/MS_supl_FigX_Inc_conc.jpg", fig_supl, width = 185, height = 140, units = "mm")
ggsave("plots/MS_supl_FigX_Inc_conc.pdf", fig_supl, width = 185, height = 140, units = "mm")
ggsave("plots/MS_supl_FigX_Inc_conc.emf", fig_supl, width = 185, height = 140, units = "mm", 
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


# statistics ####
sub_geo_stat <- data %>% 
  dplyr::select(Treatment, Day, Mn.uM) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Bir.Ac", "Bir.S", "Bir.TS", "Bir", 
                                                  "Ac", "C", "C_2"),
                            labels = c("Birnessite + acetate",  "Birnessite + sulfur + DIC",
                                       "Birnessite + thiosulfate + DIC",
                                       "Birnessite + DIC", "Acetate", "DIC", "DIC")),
         Day = factor(Day)) %>% 
  drop_na(Mn.uM)

## function for general linear hypotheses testing
calc_glht_rate <- function(input, day_sub) {
  df <- input %>% 
    subset(Day == day_sub)
  mod <- lm(Mn.uM ~ Treatment, data = df)  # linear model
  res <- glht(mod, mcp(Treatment = "Tukey"))                         # calculate statistic with glht (less sensible to different group sizes)
  s_res <- summary(res)                                            # calculate summary
  cld_s <- cld(s_res)                                              # output different letters for groups which differ significantly
  cld_let <- as.data.frame(cld_s$mcletters$Letters) %>%    # table of only letters
    rownames_to_column(var = "Treatment") %>% 
    rename(Treatment = 1, Letter = 2) %>% mutate(Day = as.numeric(day_sub))
  df <- list("model" = mod, "tested" = res, "results" =  s_res,    # combine produced data to list for tidy output of function  
             "res_display" = cld_s, "letters" = cld_let) 
  return(df)
}

stat0 <- calc_glht_rate(sub_geo_stat, "0")
stat2 <- calc_glht_rate(sub_geo_stat, "2")
stat6 <- calc_glht_rate(sub_geo_stat, "6")
stat10 <- calc_glht_rate(sub_geo_stat, "10")
stat14 <- calc_glht_rate(sub_geo_stat, "14")
stat20 <- calc_glht_rate(sub_geo_stat, "20")
stat_all <- rbind(stat0$letters, stat2$letters, 
                  stat6$letters, 
                  stat10$letters, stat14$letters, stat20$letters)

## combine letters into table
sub_geo_bar <- data %>% 
  dplyr::select(Treatment, Day, Mn.uM, Replicate) %>% 
  subset(Treatment %in% c("Bir.Ac", "Bir", "Bir.S", "Bir.TS", 
                          "Ac", "C", "C_2") & Day %in% c(0,2,6,10,14,20)) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Bir.Ac", "Bir.S", "Bir.TS", "Bir", 
                                                  "Ac", "C", "C_2"),
                            labels = c("Birnessite + acetate",  "Birnessite + sulfur + DIC",
                                       "Birnessite + thiosulfate + DIC",
                                       "Birnessite + DIC", "Acetate", "DIC", "DIC"))) %>% 
  drop_na(Mn.uM) %>% 
  left_join(stat_all, by = c("Day", "Treatment"))


# barplot with statistics ####
## calculate mean per treatment
calc_sub_Mn_stat <- sub_geo_bar %>% 
  dplyr::select(-Replicate) %>% 
  group_by(Treatment, Day) %>% 
  summarise(mean = mean(Mn.uM), .groups = "keep")

## combine
Mn_bar <- full_join(sub_geo_bar, calc_sub_Mn_stat)


man_col <- c("#ae017e", "#bcbddc", "#8c6bb1", "#440154FF","#31688EFF", "#8FD744FF")
Mn_plot_bar <- ggplot(Mn_bar, aes(x = as.factor(Day), fill = Treatment, label = Letter, group = Treatment)) +
  geom_bar(aes(y = Mn.uM), position = "dodge", stat = "summary", fun = "mean", linewidth = 0.5, color = "black") + 
  geom_point(aes(y = Mn.uM, shape = Replicate), color = "black", position = position_dodge(width = 0.9), size = 1.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,640)) +
  scale_fill_manual(values = man_col) +
  geom_text(aes(y = mean), position = position_dodge(width = 0.9), vjust = -2) +
  labs(x = "Incubation time (days)", y = "Dissolved Mn (µM)") +
  guides(fill = guide_legend(ncol = 2, title.position = "top"), 
         shape = guide_legend(ncol = 1, title.position = "top")) +
  theme_bw(base_line_size = 1, base_rect_size = 1) +
  theme(text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 9),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        axis.ticks = element_line(linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = "bottom")


ggsave("plots/MS_supl_FigSX_Mn_stat.png", plot = Mn_plot_bar, width = 185, height = 140, units = "mm")
ggsave("plots/MS_supl_FigSX_Mn_stat.pdf", plot = Mn_plot_bar, width = 185, height = 140, units = "mm")
ggsave("plots/MS_supl_FigSX_Mn_stat.emf", plot = Mn_plot_bar, width = 185, height = 140, units = "mm",
       device = {function(filename, ...) devEMF::emf(file = filename, ...)})


