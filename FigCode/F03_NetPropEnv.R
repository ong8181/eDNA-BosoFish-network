####
#### Boso Peninsula project
#### No. 4 Figure: Environmental variables and network patterns
#### 2021.10.13 Ushio
#### 2021.10.29 Ushio (R4.1.0)
#### 2021.11.10 Ushio (R4.1.2)
#### 2022.11.10 Ushio (R4.2.1)
#### 2022.11.14 Ushio: Do SMA
#### 2022.12.14 GAM revised, Ushio
#### R 4.2.1
####

# Load workspace
load("../06_CompileSmapCoefOut/06_CompileSmapCoefOut.RData")

# Load saved objects
gamm_sums <- readRDS("../08_StatisticalGeneralOut/gamm_summaries.obj")
pred_all <- readRDS("../08_StatisticalGeneralOut/pred_all.obj")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.10
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(mgcv); packageVersion("mgcv") # 1.8.40, 2022.11.10
library(GGally); packageVersion("GGally") # 2.1.2, 2021.11.1
library(glue); packageVersion("glue") # 1.6.2, 2022.12.14
library(ggtext); packageVersion("ggtext") # 0.1.2, 2022.12.14
options(tibble.print_min = 20)
options(tibble.width = Inf)
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
# Overall pattern (Main figures)
# <---------------------------------------------> #
# How to generate predictions
# # Get predicted effects of explaining variables
# # In-strength (main figure)
# pred01 <- predict(gamm01$gam, type = "response", se = T, newdata = data.frame(water_temp = seq(min(sdf$water_temp, na.rm = T), max(sdf$water_temp, na.rm = T), length.out = 200)))
# pred02 <- predict(gamm02$gam, type = "response", se = T, newdata = data.frame(sp_richness = seq(min(sdf$sp_richness, na.rm = T), max(sdf$sp_richness, na.rm = T), length.out = 200)))
# pred03 <- predict(gamm03$gam, type = "response", se = T, newdata = data.frame(total_dna_per_ml = seq(min(sdf$total_dna_per_ml, na.rm = T), max(sdf$total_dna_per_ml, na.rm = T), length.out = 200)))
# # Out-strength (main figure)
# pred04 <- predict(gamm04$gam, type = "response", se = T, newdata = data.frame(cause_temp = seq(min(sdf$cause_temp, na.rm = T), max(sdf$cause_temp, na.rm = T), length.out = 200)))
# pred05 <- predict(gamm05$gam, type = "response", se = T, newdata = data.frame(cause_richness = seq(min(sdf$cause_richness, na.rm = T), max(sdf$cause_richness, na.rm = T), length.out = 200)))
# pred06 <- predict(gamm06$gam, type = "response", se = T, newdata = data.frame(cause_totdna = seq(min(sdf$cause_totdna, na.rm = T), max(sdf$cause_totdna, na.rm = T), length.out = 200)))
# # In-strength (SI figure)
# pred07 <- predict(gamm07$gam, type = "response", se = T, newdata = data.frame(salinity = seq(min(sdf$salinity, na.rm = T), max(sdf$salinity, na.rm = T), length.out = 200)))
# pred08 <- predict(gamm08$gam, type = "response", se = T, newdata = data.frame(tide = seq(min(sdf$tide, na.rm = T), max(sdf$tide, na.rm = T), length.out = 200)))
# pred09 <- predict(gamm09$gam, type = "response", se = T, newdata = data.frame(wave_m = seq(min(sdf$wave_m, na.rm = T), max(sdf$wave_m, na.rm = T), length.out = 200)))
# # Out-strength (SI figure)
# pred10 <- predict(gamm10$gam, type = "response", se = T, newdata = data.frame(cause_salinity = seq(min(sdf$cause_salinity, na.rm = T), max(sdf$cause_salinity, na.rm = T), length.out = 200)))
# pred11 <- predict(gamm11$gam, type = "response", se = T, newdata = data.frame(cause_tide = seq(min(sdf$cause_tide, na.rm = T), max(sdf$cause_tide, na.rm = T), length.out = 200)))
# pred12 <- predict(gamm12$gam, type = "response", se = T, newdata = data.frame(cause_wave = seq(min(sdf$cause_wave, na.rm = T), max(sdf$cause_wave, na.rm = T), length.out = 200)))

# Main figures
new_var01 <- seq(min(sdf$water_temp, na.rm = T), max(sdf$water_temp, na.rm = T), length.out = 200)
new_var02 <- seq(min(sdf$sp_richness, na.rm = T), max(sdf$sp_richness, na.rm = T), length.out = 200)
new_var03 <- seq(min(sdf$total_dna_per_ml, na.rm = T), max(sdf$total_dna_per_ml, na.rm = T), length.out = 200)
new_var04 <- seq(min(sdf$cause_temp, na.rm = T), max(sdf$cause_temp, na.rm = T), length.out = 200)
new_var05 <- seq(min(sdf$cause_richness, na.rm = T), max(sdf$cause_richness, na.rm = T), length.out = 200)
new_var06 <- seq(min(sdf$cause_totdna, na.rm = T), max(sdf$cause_totdna, na.rm = T), length.out = 200)
# SI figures
new_var07 <- seq(min(sdf$salinity, na.rm = T), max(sdf$salinity, na.rm = T), length.out = 200)
new_var08 <- seq(min(sdf$tide, na.rm = T), max(sdf$tide, na.rm = T), length.out = 200)
new_var09 <- seq(min(sdf$wave_m, na.rm = T), max(sdf$wave_m, na.rm = T), length.out = 200)
new_var10 <- seq(min(sdf$cause_salinity, na.rm = T), max(sdf$cause_salinity, na.rm = T), length.out = 200)
new_var11 <- seq(min(sdf$cause_tide, na.rm = T), max(sdf$cause_tide, na.rm = T), length.out = 200)
new_var12 <- seq(min(sdf$cause_wave, na.rm = T), max(sdf$cause_wave, na.rm = T), length.out = 200)

# Generate plot for the effects of environmental variables
## water_temp; LME p = 0.2419(+), GAM p < 2e-16
p1 <- data.frame(x = new_var01, y = as.numeric(pred_all[[1]]$fit), z = as.numeric(pred_all[[1]]$se.fit)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y - 1.96*z, ymax = y + 1.96*z), alpha=0.2) +
  geom_line() + ylab("Overall effects") + 
  xlab(expression(paste("Water temperature (", degree, "C)"))) +
  ylim(0.007, 0.0200)
## sp_richness; LME p = 0.4549(+), GAM p = 0.0631
p2 <- (p1 + xlab("Species richness")) %+%
  data.frame(x = new_var02, y = as.numeric(pred_all[[2]]$fit), z = as.numeric(pred_all[[2]]$se.fit))
## total_dna; LME p = 0.153(-), GAM p < 2e-16
p3 <- (p1 + xlab("Total DNA conc. (copies/ml water)")) %+%
  data.frame(x = new_var03, y = as.numeric(pred_all[[3]]$fit), z = as.numeric(pred_all[[3]]$se.fit))
## water_temp; LME p = 0.4793(+), GAM p < 2e-16
p4 <- (p1) %+%
  data.frame(x = new_var04, y = as.numeric(pred_all[[4]]$fit), z = as.numeric(pred_all[[4]]$se.fit))
## sp_richness; LME p = 0(-), GAM p = 8.24e-07
p5 <- (p1 + xlab("Species richness")) %+%
  data.frame(x = new_var05, y = as.numeric(pred_all[[5]]$fit), z = as.numeric(pred_all[[5]]$se.fit))
## total_dna; LME p = 0.0042(-), GAM p < 2e-16
p6 <- (p1 + xlab("Total DNA conc. (copies/ml water)")) %+%
  data.frame(x = new_var06, y = as.numeric(pred_all[[6]]$fit), z = as.numeric(pred_all[[6]]$se.fit))
### Add statistical clarity information
p1 <- p1 + geom_textbox(x = min(new_var01), y = 0.0195, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* < 2.0 &times; 10<sup>-16</sup>")
p2 <- p2 + geom_textbox(x = min(new_var02), y = 0.0195, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* > 0.05")
p3 <- p3 + geom_textbox(x = min(new_var03), y = 0.0195, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* < 2.0 &times; 10<sup>-16</sup>")
p4 <- p4 + geom_textbox(x = min(new_var04), y = 0.0095, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* < 2.0 &times; 10<sup>-16</sup>")
p5 <- p5 + geom_textbox(x = min(new_var05), y = 0.0095, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* = 8.1 &times; 10<sup>-7</sup><br>
                                 GAM: *P* = 8.2 &times; 10<sup>-7</sup>")
p6 <- p6 + geom_textbox(x = min(new_var06), y = 0.0095, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* = 0.0042<br>
                                 GAM: *P* < 2.0 &times; 10<sup>-16</sup>")
## For supplementary figures
## salinity; LME p = 0.6566(+), GAM p = 2.18e-05
q1 <- data.frame(x = new_var07, y = as.numeric(pred_all[[7]]$fit), z = as.numeric(pred_all[[7]]$se.fit)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y - 1.96*z, ymax = y + 1.96*z), alpha=0.2) +
  geom_line() + ylab("Overall effects") + 
  xlab(expression(paste("Salinity (\u2030)"))) +
  ylim(0.010, 0.030)
## tide; LME p = 0(-), GAM p = 8.36e-07
q2 <- (q1 + xlab("Tide (cm)")) %+%
  data.frame(x = new_var08, y = as.numeric(pred_all[[8]]$fit), z = as.numeric(pred_all[[8]]$se.fit))
## wave_m; LME p = 0.1713(+), GAM p = 0.00234
q3 <- (q1 + xlab("Wave (m)")) %+%
  data.frame(x = new_var09, y = as.numeric(pred_all[[9]]$fit), z = as.numeric(pred_all[[9]]$se.fit))
## salinity; LME p = 0.3703(1), GAM p = 0.000619
q4 <- (q1) %+%
  data.frame(x = new_var10, y = as.numeric(pred_all[[10]]$fit), z = as.numeric(pred_all[[10]]$se.fit))
## tide; LME p = 0.5173(+), GAM p = 0.517
q5 <- (q1 + xlab("Tide (cm)")) %+%
  data.frame(x = new_var11, y = as.numeric(pred_all[[11]]$fit), z = as.numeric(pred_all[[11]]$se.fit))
## wave_m; LME p = 0.0049(+), GAM p = 6.53e-06
q6 <- (q1 + xlab("Wave (m)")) %+%
  data.frame(x = new_var12, y = as.numeric(pred_all[[12]]$fit), z = as.numeric(pred_all[[12]]$se.fit))
### Add statistical clarity information
q1 <- q1 + geom_textbox(x = min(new_var07), y = 0.0295, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* = 2.2 &times; 10<sup>-5</sup>")
q2 <- q2 + geom_textbox(x = min(new_var08), y = 0.0295, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* = 8.2 &times; 10<sup>-7</sup>
                             <br>GAM: *P* = 8.4 &times; 10<sup>-7</sup>")
q3 <- q3 + geom_textbox(x = min(new_var09), y = 0.0295, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* = 0.00234")
q4 <- q4 + geom_textbox(x = min(new_var10), y = 0.0295, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* = 6.2 &times; 10<sup>-4</sup>")
q5 <- q5 + geom_textbox(x = min(new_var11), y = 0.0295, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* > 0.05<br>
                                 GAM: *P* > 0.05")
q6 <- q6 + geom_textbox(x = min(new_var12), y = 0.0295, hjust = 0, vjust = 1, box.size = NA,
                        label = "LME: *P* = 0.0049<br>
                                 GAM: *P* = 6.5 &times; 10<sup>-6</sup>")

# Save as R object
saveRDS(list(p1, p2, p3, p4, p5, p6), sprintf("%s/Fig_OverallEffect1.obj", fig_folder1))
saveRDS(list(q1, q2, q3, q4, q5, q6), sprintf("%s/Fig_OverallEffect2.obj", fig_folder1))
# ggsave(filename = sprintf("%s/Fig_Salinity.pdf", fig_folder1), plot = q1,
#        width = 6, height = 6, device = cairo_pdf)

# Output statistic tables
gamm_stats <- data.frame(model = rep(NA,12), linear_Tvalue = rep(NA,12), linear_Pvalue = rep(NA,12),
                         gam_edf = rep(NA,12), gam_Fvalue = rep(NA,12), gam_Pvalue = rep(NA,12))
gamm_stats[,1] <- c("In-IS, water_temp", "In-IS, richness", "In-IS, total DNA",
                    "Out-IS, water_temp", "Out-IS, richness", "Out-IS, total DNA",
                    "In-IS, salinity", "In-IS, tide", "In-IS, wave",
                    "Out-IS, salinity", "Out-IS, tide", "Out-IS, wave")
for(i in 1:12) {
  gamm_stats[i,2:3] <- gamm_sums[[i]][[1]]$tTable[2,c("t-value","p-value")]
  gamm_stats[i,4:6] <- gamm_sums[[i]][[2]]$s.table[1,c("edf","F", "p-value")]
}
write.csv(gamm_stats, sprintf("%s/Table_GAMMstats.csv", fig_folder1))


# <---------------------------------------------> #
# Overall pattern (For supplementary figures)
# <---------------------------------------------> #
# Water temperature, species richness, total eDNA concentrations
## In-strength
s1 <- sdf %>% ggplot(aes(x = water_temp, y = abs(IS), color = site_code)) +
  geom_point(alpha = 0.05) + scale_color_manual(values = palette_custom, name = "Site ID") +
  stat_smooth(method = "gam", se = FALSE) + coord_cartesian(ylim = c(1e-06, 1e01)) +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  xlab(expression(paste("Water temperature (", degree, "C)"))) +
  ylab("Interaction strength") +
  theme_cowplot() +
  #guides(fill = guide_legend(override.aes = list(fill = NA))) +
  NULL
s2 <- s1 + aes(x = sp_richness) + xlab("Species richness") 
s3 <- s1 + aes(x = total_dna_per_ml) +  
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  stat_smooth(method = "gam", se = FALSE) + coord_cartesian(ylim = c(1e-06, 1e01)) +
  xlab("Total eDNA conc. (copies/ml water)")
## Out-strength
s4 <- s1 + aes(x = cause_temp) + xlab(expression(paste("Water temperature (", degree, "C)")))
s5 <- s1 + aes(x = cause_richness) + xlab("Species richness") + theme_cowplot()
s6 <- s1 + aes(x = cause_totdna) +  
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  stat_smooth(method = "gam", se = FALSE) + coord_cartesian(ylim = c(1e-06, 1e01)) +
  xlab("Total eDNA conc. (copies/ml water)")
s_legend <- get_legend(s1)

## Prepare sub-titles
title_inst <- ggdraw() + draw_label("In-strength", fontface = 'bold', size = 15, x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
title_otst <- ggdraw() + draw_label("Out-strength", fontface = 'bold', size = 15, x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

## Assemble panels
legend_s <- get_legend(s1)
s_sub_1 <- plot_grid(s1 + theme(legend.position = "none"),
                     s2 + theme(legend.position = "none"),
                     s3 + theme(legend.position = "none"),
                     align = "hv", nrow = 1, labels = c("a", "b", "c"))
s_sub_2 <- plot_grid(s4 + theme(legend.position = "none"),
                     s5 + theme(legend.position = "none"),
                     s6 + theme(legend.position = "none"),
                     align = "hv", nrow = 1, labels = c("d", "e", "f"))
s_all_1 <- plot_grid(title_inst, s_sub_1, title_otst, s_sub_2,
                     nrow = 4, rel_heights = c(0.2, 1, 0.2, 1))
s_all <- plot_grid(s_all_1, NULL, s_legend, nrow = 1, rel_widths =  c(1,0.02,0.1))
ggsave(sprintf("%s/JPG_Overall_pattern_SI1.jpg", fig_folder1), s_all,
       width = 14, height = 9, dpi = 200)




## Overall pattern (For SI)
## Salinity, tide, and wave
## In-strength
t1 <- s1 + aes(x = salinity) + xlab(expression(paste("Salinity (\u2030)")))
t2 <- s1 + aes(x = tide) + xlab("Tide (cm)")
t3 <- s1 + aes(x = wave_m) + xlab("Wave (m)")
## Out-strength
t4 <- s1 + aes(x = cause_salinity) + xlab(expression(paste("Salinity (\u2030)")))
t5 <- s1 + aes(x = cause_tide) + xlab("Tide (cm)")
t6 <- s1 + aes(x = cause_wave) + xlab("Wave (m)")
# Assemble figures
t_sub_1 <- plot_grid(t1 + theme(legend.position = "none"),
                     t2 + theme(legend.position = "none"),
                     t3 + theme(legend.position = "none"),
                     align = "hv", nrow = 1, labels = c("a", "b", "c"))
t_sub_2 <- plot_grid(t4 + theme(legend.position = "none"),
                     t5 + theme(legend.position = "none"),
                     t6 + theme(legend.position = "none"),
                     align = "hv", nrow = 1, labels = c("d", "e", "f"))
t_all_1 <- plot_grid(title_inst, t_sub_1, title_otst, t_sub_2,
                     nrow = 4, rel_heights = c(0.2, 1, 0.2, 1))
t_all <- plot_grid(t_all_1, NULL, s_legend, nrow = 1, rel_widths =  c(1,0.02,0.1))
ggsave(sprintf("%s/JPG_Overall_pattern_SI2.jpg", fig_folder1), t_all,
       width = 14, height = 9, dpi = 200)



# <---------------------------------------------> #
# Correlations among properties
# <---------------------------------------------> #
# Correlations between network properties
sdf_sample0 <- sdf %>% group_by(sample_id) %>% summarize_all(mean)
sdf_sample <- sdf %>% select(sample_id, time_code, site_code, site_name, unique_code, date, time,
                             weather, water_color, water_q_tool, sampling_person) %>%
  group_by(sample_id) %>% summarize_all(unique)
sdf_sample <- cbind(sdf_sample, sdf_sample0 %>%
                      select(lat_n, long_e, water_temp, salinity,
                             tide, linear_dist, total_reads,
                             total_dna_estimated, total_dna_per_ml,
                             effect_var_edna, sp_richness,
                             IS, cause_var_edna, n_int))
rm(sdf_sample0)
sdf_comp <- sdf_sample %>% select(IS, n_int, total_dna_per_ml, water_temp, sp_richness)

# Check correlations
## SMA
## 1st column
sma_lm01 <- sdf_comp %>% smatr::sma("log10(total_dna_per_ml) ~ water_temp", data = .)
slope01 <- sma_lm01$coef[[1]][2,1]; intct01 <- sma_lm01$coef[[1]][1,1] # p <2e-16
sma_lm02 <- sdf_comp %>% smatr::sma("sp_richness ~ water_temp", data = .)
slope02 <- sma_lm02$coef[[1]][2,1]; intct02 <- sma_lm02$coef[[1]][1,1] # p <2e-16
sma_lm03 <- sdf_comp %>% smatr::sma("log10(abs(IS)) ~ water_temp", data = .)
slope03 <- sma_lm03$coef[[1]][2,1]; intct03 <- sma_lm03$coef[[1]][1,1] # p = 0.184
sma_lm04 <- sdf_comp %>% smatr::sma("n_int ~ water_temp", data = .)
slope04 <- sma_lm04$coef[[1]][2,1]; intct04 <- sma_lm04$coef[[1]][1,1] # p <2e-16
## 2nd column
sma_lm05 <- sdf_comp %>% smatr::sma("sp_richness ~ log10(total_dna_per_ml)", data = .)
slope05 <- sma_lm05$coef[[1]][2,1]; intct05 <- sma_lm05$coef[[1]][1,1] # p = 0.0159
sma_lm06 <- sdf_comp %>% smatr::sma("log10(abs(IS)) ~ log10(total_dna_per_ml)", data = .)
slope06 <- sma_lm06$coef[[1]][2,1]; intct06 <- sma_lm06$coef[[1]][1,1] # p = 0.0201
sma_lm07 <- sdf_comp %>% smatr::sma("n_int ~ log10(total_dna_per_ml)", data = .)
slope07 <- sma_lm07$coef[[1]][2,1]; intct07 <- sma_lm07$coef[[1]][1,1] # p = 6.28e-06
## 3rd column
sma_lm08 <- sdf_comp %>% smatr::sma("log10(abs(IS)) ~ sp_richness", data = .)
slope08 <- sma_lm08$coef[[1]][2,1]; intct08 <- sma_lm08$coef[[1]][1,1] # p = 6.23e-05
sma_lm09 <- sdf_comp %>% smatr::sma("n_int ~ sp_richness", data = .)
slope09 <- sma_lm09$coef[[1]][2,1]; intct09 <- sma_lm09$coef[[1]][1,1] # p <2e-16
## 4th column
sma_lm10 <- sdf_comp %>% smatr::sma("n_int ~  log10(abs(IS))", data = .)
slope10 <- sma_lm10$coef[[1]][2,1]; intct10 <- sma_lm10$coef[[1]][1,1] # p = 0.0023

# Column one: Water temperature
tx_size <- 7.5
g1_1 <- ggplot(sdf_comp, aes(x = water_temp, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab("Scaled density") + xlab(NULL) + theme(text = element_text(size = tx_size))
g1_2 <- ggplot(sdf_comp, aes(x = water_temp, y = total_dna_per_ml)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = intct01, slope = slope01, linetype = 1, color = "red3") + xlab(NULL) +
  ylab("Total eDNA concentrations\n(copies/ml water/sample)") +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  theme(text = element_text(size = tx_size)) +
  NULL
g1_3 <- ggplot(sdf_comp, aes(x = water_temp, y = sp_richness)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = intct02, slope = slope02, linetype = 1, color = "red3") + xlab(NULL) +
  ylab("Species richness") +
  theme(text = element_text(size = tx_size)) +
  NULL
g1_4 <- ggplot(sdf_comp, aes(x = water_temp, y = abs(IS))) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = intct03, slope = slope03, linetype = 1, color = "gray80") + xlab(NULL) +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  ylab("Mean interaction strength") +
  theme(text = element_text(size = tx_size)) +
  NULL
g1_5 <- ggplot(sdf_comp, aes(x = water_temp, y = n_int)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = intct04, slope = slope04, linetype = 1, color = "red3") +
  xlab(expression(paste("Water temperature (", degree, "C)"))) +
  ylab("The number of\ninterspecific interactions") +
  theme(text = element_text(size = tx_size)) +
  NULL

# Column two: Total eDNA concentrations
g2_2 <- ggplot(sdf_comp, aes(x = total_dna_per_ml, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  ylab(NULL) + xlab(NULL) + theme(text = element_text(size = tx_size))
NULL
g2_3 <- ggplot(sdf_comp, aes(x = total_dna_per_ml, y = sp_richness)) +
  geom_point(alpha = 0.3) +
  xlab(NULL) + ylab(NULL) +
  geom_abline(intercept = intct05, slope = slope05, linetype = 1, color = "red3") +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  theme(text = element_text(size = tx_size)) +
  NULL
g2_4 <- ggplot(sdf_comp, aes(x = total_dna_per_ml, y = abs(IS))) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  xlab(NULL) + ylab(NULL) +
  geom_abline(intercept = intct06, slope = slope06, linetype = 1, color = "red3") +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  theme(text = element_text(size = tx_size)) +
  NULL
g2_5 <- ggplot(sdf_comp, aes(x = total_dna_per_ml, y = n_int)) +
  geom_point(alpha = 0.3) +
  xlab("Total eDNA concentrations\n(copies/ml water/sample)") +
  ylab(NULL) +
  geom_abline(intercept = intct07, slope = slope07, linetype = 1, color = "red3") +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  theme(text = element_text(size = tx_size)) +
  NULL

# Column three: the number of detected species per sample
g3_3 <- ggplot(sdf_comp, aes(x = sp_richness, ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  ylab(NULL) + xlab(NULL) + theme(text = element_text(size = tx_size)) +
  NULL
g3_4 <- ggplot(sdf_comp, aes(x = sp_richness, y = abs(IS))) +
  geom_point(alpha = 0.3) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
  geom_abline(intercept = intct08, slope = slope08, linetype = 1, color = "red3") +
  theme(text = element_text(size = tx_size)) +
  NULL
g3_5 <- ggplot(sdf_comp, aes(x = sp_richness, y = n_int)) +
  geom_point(alpha = 0.3) +
  xlab("Species richness") +
  ylab(NULL) +
  geom_abline(intercept = intct09, slope = slope09, linetype = 1, color = "red3") +
  theme(text = element_text(size = tx_size)) +
  NULL

# Column four: the number of detected species per sample
g4_4 <- ggplot(sdf_comp, aes(x = abs(IS), ..scaled..)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  ylab(NULL) + xlab(NULL) + theme(text = element_text(size = tx_size)) +
  NULL
g4_5 <- ggplot(sdf_comp, aes(x = abs(IS), y = n_int)) +
  geom_point(alpha = 0.3) +
  xlab("Mean interaction strength") +
  ylab(NULL) +
  scale_x_continuous(trans = "log10", labels = macam::label_10_to_power) +
  geom_abline(intercept = intct10, slope = slope10, linetype = 1, color = "red3") +
  theme(text = element_text(size = tx_size)) +
  NULL

# Column five: the number of detected species per sample
g5_5 <- ggplot(sdf_comp, aes(x = n_int, ..scaled..)) +
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

