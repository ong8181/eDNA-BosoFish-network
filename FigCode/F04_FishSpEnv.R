####
#### Boso Peninsula project
#### Temperature dependence of fish interactions
#### - Species-specific pattern
#### 2022.12.14 GAM revised, Ushio
####

#setwd("FigCode/")

# Load workspace
load("../09_StatisticalFishSpOut/09_StatisticalFishSpOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.11
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.17
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(mgcv); packageVersion("mgcv") # 1.8.40, 2022.5.17
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"

# Prepare color palette
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
get_palette2 <- colorRampPalette(brewer.pal(8, "Set1"))
palette_custom <- get_palette(11)[c(1:4,7:11,6:5)]
palette_custom3 <- get_palette2(28)


# <---------------------------------------------> #
#  Check statistical tests
# <---------------------------------------------> #
water_temp_gamm1
water_temp_gamm2
richness_gamm1 # For supplement
richness_gamm2 # For supplement
total_edna_gamm1 # For supplement
total_edna_gamm2 # For supplement
# Set threshol p-value for the visualization in the main figures
pth <- 0.0001
# Fish sp name labels
fsp <- tax_df$scientific_name[match(top_taxa, rownames(tax_df))]
names(fsp) <- top_taxa
# Sort fish name alphabetically
sdf$effect_var <-  factor(sdf$effect_var, levels = names(sort(fsp)))
sdf$cause_var <-  factor(sdf$cause_var, levels = names(sort(fsp)))


# <---------------------------------------------> #
#  Visualization: Water temperature effects
# <---------------------------------------------> #
# Water temperature effects
# (only significant species for the main figure)
fishsp_water1 <- water_temp_gamm1 %>% filter(gam_pval < pth | lme_pval < pth) %>% 
  select(fish_taxa) %>% pull; length(fishsp_water1)
fishsp_water2 <- water_temp_gamm2 %>% filter(gam_pval < pth | lme_pval < pth) %>% 
  select(fish_taxa) %>% pull; length(fishsp_water2)
# In-strength
g1_1 <- sdf %>% 
  filter(effect_var %in% fishsp_water1) %>% droplevels() %>% 
  ggplot(aes(x = water_temp, y = abs(IS), group = effect_var, color = site_code)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "gam", se = FALSE, color = "gray30") +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  facet_wrap(~effect_var, ncol = 5, scales = "free", labeller = labeller(effect_var = fsp[fishsp_water1])) +
  scale_color_manual(values = palette_custom, name = "Site ID") +
  panel_border() +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Interaction strength (in-strength)") +
  xlab(expression(paste("Water temperature (", degree, "C)"))) +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic", size = 10)) +
  NULL
# Out-strength
g1_2 <- sdf %>% 
  filter(cause_var %in% fishsp_water2) %>% droplevels() %>% 
  ggplot(aes(x = cause_temp, y = abs(IS), group = cause_var, color = site_code)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "gam", se = FALSE, color = "gray30") +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  facet_wrap(~cause_var, ncol = 5, scales = "free", labeller = labeller(cause_var = fsp[fishsp_water2])) +
  scale_color_manual(values = palette_custom, name = "Site ID") +
  panel_border() +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Interaction strength (out-strength)") +
  xlab(expression(paste("Water temperature (", degree, "C)"))) +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic", size = 10)) +
  NULL


# <---------------------------------------------> #
#  Visualization: Richness effects
# <---------------------------------------------> #
# (only significant species for the main figure)
fishsp_rich1 <- richness_gamm1 %>% filter(gam_pval < pth | lme_pval < pth) %>% 
  select(fish_taxa) %>% pull; length(fishsp_rich1)
fishsp_rich2 <- richness_gamm2 %>% filter(gam_pval < pth | lme_pval < pth) %>% 
  select(fish_taxa) %>% pull; length(fishsp_rich2)
# In-strength
g2_1 <- sdf %>% 
  filter(effect_var %in% fishsp_rich1) %>% droplevels() %>%  
  ggplot(aes(x = sp_richness, y = abs(IS), group = effect_var, color = site_code)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "gam", se = FALSE, color = "gray30") +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  facet_wrap(~effect_var, ncol = 5, scales = "free", labeller = labeller(effect_var = fsp[fishsp_rich1])) +
  scale_color_manual(values = palette_custom, name = "Site ID") +
  panel_border() +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Interaction strength (in-strength)") +
  xlab("Species richness") +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic", size = 10)) +
  NULL
# Out-strength
g2_2 <- sdf %>% 
  filter(cause_var %in% fishsp_rich2) %>% droplevels() %>%  
  ggplot(aes(x = cause_richness, y = abs(IS), group = cause_var, color = site_code)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "gam", se = FALSE, color = "gray30") +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  facet_wrap(~cause_var, ncol = 5, scales = "free", labeller = labeller(cause_var = fsp[fishsp_rich2])) +
  scale_color_manual(values = palette_custom, name = "Site ID") +
  panel_border() +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Interaction strength (out-strength)") +
  xlab("Species richness") +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic", size = 10)) +
  NULL


# <---------------------------------------------> #
#  Visualization: Total abundance effects
# <---------------------------------------------> #
# Water temperature effects
# (only significant species for the main figure)
fishsp_dna1 <- total_edna_gamm1 %>% filter(gam_pval < pth | lme_pval < pth) %>% 
  select(fish_taxa) %>% pull; length(fishsp_dna1)
fishsp_dna2 <- total_edna_gamm2 %>% filter(gam_pval < pth | lme_pval < pth) %>% 
  select(fish_taxa) %>% pull; length(fishsp_dna2)
# In-strength
g3_1 <- sdf %>% 
  filter(effect_var %in% fishsp_dna1) %>% droplevels() %>%  
  ggplot(aes(x = total_dna_per_ml, y = abs(IS), group = effect_var, color = site_code)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "gam", se = FALSE, color = "gray30") +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  facet_wrap(~effect_var, ncol = 5, scales = "free", labeller = labeller(effect_var = fsp[fishsp_dna1])) +
  scale_color_manual(values = palette_custom, name = "Site ID") +
  panel_border() +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Interaction strength (in-strength)") +
  xlab("Total eDNA conc. (copies/ml water)") +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic", size = 10)) +
  NULL
# Out-strength
g3_2 <- sdf %>% 
  filter(cause_var %in% fishsp_dna2) %>% droplevels() %>%  
  ggplot(aes(x = cause_totdna, y = abs(IS), group = cause_var, color = site_code)) +
  geom_point(alpha = 0.3) + stat_smooth(method = "gam", se = FALSE, color = "gray30") +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  facet_wrap(~cause_var, ncol = 5, scales = "free", labeller = labeller(cause_var = fsp[fishsp_dna2])) +
  scale_color_manual(values = palette_custom, name = "Site ID") +
  panel_border() +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  ylab("Interaction strength (out-strength)") +
  xlab("Total eDNA conc. (copies/ml water)") +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic", size = 10)) +
  NULL


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save figures
## Water temperature effects
ggsave(file = sprintf("%s/PDF_ISvsTemp_inIS.pdf", fig_folder1), plot = g1_1, width = 14, height = 11)
ggsave(file = sprintf("%s/PDF_ISvsTemp_outIS.pdf", fig_folder1), plot = g1_2, width = 10, height = 5)
## Species richness effects
ggsave(file = sprintf("%s/PDF_ISvsRich_inIS.pdf", fig_folder1), plot = g2_1, width = 14, height = 9)
ggsave(file = sprintf("%s/PDF_ISvsRich_outIS.pdf", fig_folder1), plot = g2_2, width = 4, height = 5)
## Total DNA concentration effects
ggsave(file = sprintf("%s/PDF_ISvsDNA_inIS.pdf", fig_folder1), plot = g3_1, width = 14, height = 11)
ggsave(file = sprintf("%s/PDF_ISvsDNA_outIS.pdf", fig_folder1), plot = g3_2, width = 10, height = 5)
## Output legend
legend01 <- get_legend(g1_1)
ggsave(file = sprintf("%s/Site_legend.pdf", fig_folder1),
       plot = legend01, width = 3.6, height = 1)

# Save objects
saveRDS(list(g1_1, g1_2), sprintf("%s/Fig_ISvsTemp.obj", fig_folder1))
saveRDS(list(g2_1, g2_2), sprintf("%s/Fig_ISvsRich.obj", fig_folder1))
saveRDS(list(g3_1, g3_2), sprintf("%s/Fig_ISvsDNA.obj", fig_folder1))

write.csv(cbind(opt_mdr_res[,1:11],
                tax_df[match(opt_mdr_res$effect_var, rownames(tax_df)),],
                water_temp_gamm1),
          sprintf("%s/CSV_MDRperformance.csv", fig_folder1), row.names = F)

