## Potter Cove previous dataset from Monien et al. 2014 publication
 # Monien, Patrick; Lettmann, Karsten Alexander; Monien, Donata; Asendorf, Sanja; Wölfl, Anne-Cathrin; Lim, Chai Heng; Thal, Janis; Schnetger, Bernhard; Brumsack, Hans-J?rgen (2014): Sediment and pore water geochemistry of sediment cores, Potter Cove, King George Island. PANGAEA, https://doi.org/10.1594/PANGAEA.832335, Supplement to: Monien, P et al. (2014): Redox conditions and trace metal cycling in coastal sediments from the maritime Antarctic. Geochimica et Cosmochimica Acta, 141, 26-44, https://doi.org/10.1016/j.gca.2014.06.003

## combine individual tables
 # plot 

# load packages ####
library(tidyverse)
library(patchwork)
library(viridis)

# import data ####
# save and load workspace
## path to working directory
wd <- "WORKINGDIRECTORY/"
setwd(wd)

## import solid phase data
sp <- grep("\\*/", read_lines("PC-P01b_geochem.tab"), value = F) %>% 
  read_tsv("PC-P01b_geochem.tab", skip = .)

## import porewater data
pq <- grep("\\*/", read_lines("PC-P01a_geochem_pw.tab"), value = F) %>% 
  read_tsv("PC-P01a_geochem_pw.tab", skip = .)


# calulcate wt% into mol/cm^3 ####
## new columns in mol/g dw instead of wt%
molecular_weight_MnO <- 70.9374 # for MnO, Monien data, g/mol = mg/mmol
water_density <- 1 # fresh water, g/cm^3
dry_sed_density <- 2.6 # assumed value, g/cm^3

sp_calc <- sp %>% 
  mutate(MnO_umol_per_g_sed_saltcor = (`MnO [%] (salt corrected, mass percenta...)`*10)/molecular_weight_MnO*1000,
         MnO_umol_per_g_sed = (`MnO [%] (mass percentages, Wave-length...)`*10)/molecular_weight_MnO*1000,
         sediment_wm_perc = 100 - `Water wm [%] (calculated using difference b...)`,
         dry_sed_per_1g_wet_sed_g = 1/100 * sediment_wm_perc,
         water_per_1g_wet_sed_g = 1/100 * `Water wm [%] (calculated using difference b...)`,
         dry_sed_per_1g_wet_sed_cm3 = dry_sed_per_1g_wet_sed_g/dry_sed_density,
         water_per_1g_wet_sed_cm3 = water_per_1g_wet_sed_g/water_density,
         vol_1g_wet_sed_cm3 = dry_sed_per_1g_wet_sed_cm3 + water_per_1g_wet_sed_cm3,
         vol_dry_sed_in_1cm3_wet_sed_cm3 = (dry_sed_per_1g_wet_sed_cm3*1)/vol_1g_wet_sed_cm3,
         dry_sed_per_1cm3_wet_sed_g = vol_dry_sed_in_1cm3_wet_sed_cm3 * dry_sed_density,
         
         MnO_umol_per_cm3 = MnO_umol_per_g_sed * dry_sed_per_1cm3_wet_sed_g
         )

#
# Mn-MS supl plot ####
Mn_mol_MnMS <- sp_calc %>%
  ggplot(aes(y = `Depth sed [m]`*100, x = MnO_umol_per_cm3)) +
  geom_point(size = 1.5) +
  geom_path() +
  scale_y_reverse(limits = c(21, -0.6), expand = c(0,0), name = "Sediment depth (cm)") +
  scale_x_continuous(name = "MnO (µmol/cm³)", limits = c(-0.1,21), breaks = c(0,10,20),
                     position = "top", expand = c(0,0)) +
  labs(title = "PC-P01 Monien et al. 2014") +
  theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(panel.grid.minor.y = element_blank(),
        text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.ticks = element_line(linewidth = 0.5))


Mn_perc_MnMS <- sp_calc %>%
  ggplot(aes(y = `Depth sed [m]`*100, x = `MnO [%] (salt corrected, mass percenta...)`)) +
  geom_point(size = 1.5) +
  geom_path() +
  scale_y_reverse(limits = c(21, -0.6), expand = c(0,0), name = "Sediment depth (cm)") +
  scale_x_continuous(name = "MnO (wt.%)", limits = c(-0.001,0.17), breaks = c(0, 0.1),
                     position = "top", expand = c(0,0)) +
  labs(title = "PC-P01 Monien et al. 2014") +
  theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(panel.grid.minor.y = element_blank(),
        text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.ticks = element_line(linewidth = 0.5))

disMn_MnMS <- pq %>% 
  ggplot(aes(y = `Depth sed [m]`*100, x = `Mn [Âµmol/l]`)) +
  geom_point(size = 1.5) +
  geom_path() +
  scale_y_reverse(limits = c(21, -0.6), expand = c(0,0), name = "Sediment depth (cm)") +
  scale_x_continuous(name = "dissolved Mn (µM)", limits = c(-0.1,110), breaks = c(0,50,100),
                     position = "top", expand = c(0,0)) +
  labs(title = "PC-P01 Monien et al. 2014") +
  theme_bw(base_line_size = 1, base_rect_size = 1) + 
  theme(panel.grid.minor.y = element_blank(),
        text = element_text(size = 10, color = "black"), 
        axis.title =  element_text(size = 8),
        axis.text = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.ticks = element_line(linewidth = 0.5))

MS_comb <- disMn_MnMS + Mn_perc_MnMS + Mn_mol_MnMS & theme(plot.title = element_blank()) &
  plot_annotation(tag_levels = list(c("B", "C", "D")))

ggsave("part_Monien_MnO_disMn.pdf", MS_comb, width = 116, height = 65, units = "mm")
ggsave("part_Monien_MnO_disMn.jpg", MS_comb, width = 116, height = 65, units = "mm")
