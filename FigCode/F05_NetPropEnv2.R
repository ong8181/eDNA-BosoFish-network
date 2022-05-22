####
#### Boso Peninsula project
#### No.6: Network properties and environmental variables
####

# Load workspace
load("../06_VisualizeNetwork2Out/06_VisualizeNetwork2Out.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.10.13
library(GGally); packageVersion("GGally") # 2.1.2, 2021.11.1
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.10.13
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
options(tibble.print_min = Inf)
options(tibble.width = Inf)
theme_set(theme_bw())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"
fig_subfolder1 <- sprintf("%s/Fig_SiteNetworks", fig_folder1)
dir.create(fig_subfolder1)

# Prepare color palette
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
get_palette2 <- colorRampPalette(brewer.pal(8, "Set1"))
palette_custom <- rev(get_palette(11))
palette_custom2 <- get_palette(28)
palette_custom3 <- get_palette2(28)


# <---------------------------------------------> #
# Prepare information for all sites
# <---------------------------------------------> #
# Add edge strength information
same_site_id <- which(sapply(str_split(uic_edna_st_strong$cause_var, "_"), "[", 2) == sapply(str_split(uic_edna_st_strong$effect_var, "_"), "[", 2))
uic_edna_same_st <- uic_edna_st_strong[same_site_id,]

# Assign information for uic_edna_same_st
## Site ID
uic_edna_same_st$site_id <- sapply(str_split(uic_edna_same_st$cause_var, "_"), "[", 2)
## eDNA abundance for each species and site
uic_edna_same_st$effect_sp_abundance <- colSums(asv_df_sitesp[,uic_edna_same_st$effect_var], na.rm = T)
uic_edna_same_st$cause_sp_abundance <- colSums(asv_df_sitesp[,uic_edna_same_st$cause_var], na.rm = T)

# Generate site-specific summary table
all(names(rowSums(asv_df_conv > 0)) == rownames(sample_df))
sample_df$total_div <- rowSums(asv_df_conv > 0)
site_env_df <- sample_df %>% group_by(site_code) %>%
  summarize(max_temp = max(water_temp, na.rm = T),
            med_temp = median(water_temp, na.rm = T),
            mean_temp = mean(water_temp, na.rm = T),
            min_temp = min(water_temp, na.rm = T),
            mean_total_edna = mean(total_dna_per_ml, na.rm = T),
            sd_total_edna = sd(total_dna_per_ml, na.rm = T),
            mean_div_per_sample = mean(total_div, na.rm = T),
            sd_div_per_sample = sd(total_div, na.rm = T),
            NULL)
# Add site richness information
site_summary_df <- readRDS("../07_EnvironmentNetworkOut/site_summary_df.obj")

# Add the number of link, and total eDNA information
tmp <- uic_edna_same_st %>% group_by(site_id) %>%
  summarize(n_link = n(),
            sum_cause_edna = sum(cause_sp_abundance),
            sum_effect_edna = sum(effect_sp_abundance))
site_env_df$n_link <- tmp$n_link
site_env_df$sum_cause_edna <- tmp$sum_cause_edna
site_env_df$sum_effect_edna <- tmp$sum_effect_edna; rm(tmp)
site_env_df$site_total_richness <- site_summary_df$richness
site_env_df$site_total_edna <- site_summary_df$total_edna

# ----------------------------------------------- #
# Visualize results
# ----------------------------------------------- #
# Extract subset data frame
site_env_df1 <- site_env_df %>% select(med_temp,
                                       mean_total_edna,
                                       mean_div_per_sample,
                                       site_total_richness,
                                       n_link)

# Column one: median water temperature
tx_size <- 10
g1_1 <- ggplot(site_env_df1, aes(x = med_temp, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab("Scaled density") + xlab(NULL) + theme(text = element_text(size = tx_size))
g1_2 <- ggplot(site_env_df1, aes(x = med_temp, y = mean_total_edna)) +
  geom_point() + xlab(NULL) +
  ylab("Mean total eDNA concentrations\n(copies/ml water/sample)") +
  theme(text = element_text(size = tx_size)) +
  NULL
g1_3 <- ggplot(site_env_df1, aes(x = med_temp, y = mean_div_per_sample)) +
  geom_point() + xlab(NULL) +
  stat_smooth(method = "lm", color = "black", size = 0.5) +
  ylab("The number of\ndetected species per sample") +
  theme(text = element_text(size = tx_size)) +
  ylim(30, 55) +
  NULL
g1_4 <- ggplot(site_env_df1, aes(x = med_temp, y = site_total_richness)) +
  geom_point() + xlab(NULL) +
  stat_smooth(method = "lm", color = "black", linetype = 2, size = 0.5) +
  ylab("The number of\ndetected species per site") +
  theme(text = element_text(size = tx_size)) +
  ylim(220, 320) +
  NULL
g1_5 <- ggplot(site_env_df1, aes(x = med_temp, y = n_link)) +
  geom_point() +
  xlab(expression(paste("Median water temperature (", degree, "C)"))) +
  ylab("The number of\ninterspecific interactions") +
  theme(text = element_text(size = tx_size)) +
  ylim(600, 850) +
  NULL

# Column two: mean total eDNA concentrations
g2_2 <- ggplot(site_env_df1, aes(x = mean_total_edna, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab(NULL) + xlab(NULL) + theme(text = element_text(size = tx_size))
  NULL
g2_3 <- ggplot(site_env_df1, aes(x = mean_total_edna, y = mean_div_per_sample)) +
  geom_point() + xlab(NULL) +
  ylab(NULL) +
  theme(text = element_text(size = tx_size)) +
  ylim(30, 55) +
  NULL
g2_4 <- ggplot(site_env_df1, aes(x = mean_total_edna, y = site_total_richness)) +
  geom_point() + xlab(NULL) +
  ylab(NULL) +
  theme(text = element_text(size = tx_size)) +
  ylim(220, 320) +
  NULL
g2_5 <- ggplot(site_env_df1, aes(x = mean_total_edna, y = n_link)) +
  geom_point() +
  xlab("Mean total eDNA concentrations\n(copies/ml water/sample)") +
  ylab(NULL) +
  theme(text = element_text(size = tx_size)) +
  ylim(600, 850) +
  NULL

# Column three: the number of detected species per sample
g3_3 <- ggplot(site_env_df1, aes(x = mean_div_per_sample, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab(NULL) + xlab(NULL) + theme(text = element_text(size = tx_size)) +
  NULL
g3_4 <- ggplot(site_env_df1, aes(x = mean_div_per_sample, y = site_total_richness)) +
  geom_point() + xlab(NULL) +
  ylab(NULL) +
  theme(text = element_text(size = tx_size)) +
  ylim(220, 320) +
  NULL
g3_5 <- ggplot(site_env_df1, aes(x = mean_div_per_sample, y = n_link)) +
  geom_point() +
  xlab("The number of\ndetected species per sample") +
  ylab(NULL) +
  theme(text = element_text(size = tx_size)) +
  ylim(600, 850) +
  NULL

# Column four: the number of detected species per sample
g4_4 <- ggplot(site_env_df1, aes(x = site_total_richness, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab(NULL) + xlab(NULL) + theme(text = element_text(size = tx_size)) +
  NULL
g4_5 <- ggplot(site_env_df1, aes(x = site_total_richness, y = n_link)) +
  geom_point() +
  xlab("The number of\ndetected species per site") +
  ylab(NULL) +
  stat_smooth(method = "lm", color = "black", linetype = 1, size = 0.5) +
  theme(text = element_text(size = tx_size)) +
  ylim(600, 850) +
  NULL

# Column five: the number of detected species per sample
g5_5 <- ggplot(site_env_df1, aes(x = n_link, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab(NULL) + xlab("The number of\ninterspecific interactions") + theme(text = element_text(size = tx_size)) +
  NULL

# Combine all panels
if(F){
  g1_1 <- g1_1 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
  g1_2 <- g1_2 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
  g1_3 <- g1_3 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
  g1_4 <- g1_4 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
  g1_5 <- g1_5 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
  
  g2_2 <- g2_2 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g2_3 <- g2_3 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g2_4 <- g2_4 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g2_5 <- g2_5 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  g3_3 <- g3_3 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g3_4 <- g3_4 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g3_5 <- g3_5 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  g4_4 <- g4_4 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g4_5 <- g4_5 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  g5_5 <- g5_5 + theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

# Combine all figures
g_all <- plot_grid(g1_1, g1_2, g1_3, g1_4, g1_5,
                   NULL, g2_2, g2_3, g2_4, g2_5,
                   NULL, NULL, g3_3, g3_4, g3_5,
                   NULL, NULL, NULL, g4_4, g4_5,
                   NULL, NULL, NULL, NULL, g5_5,
                   ncol = 5, byrow = F,
                   #rel_widths = c(1,1,1,1,2),
                   align = "hv",
                   axis = "lrbt")

# Save figure and object
ggsave(sprintf("%s/PDF_CorPairs.pdf", fig_folder1),
       plot = g_all, width = 12, height = 12)
saveRDS(g_all, sprintf("%s/Fig_CorPairs.obj", fig_folder1))

