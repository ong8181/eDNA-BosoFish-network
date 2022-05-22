####
#### Boso Peninsula project
#### No. 5 Figure: Environmental variables and network patterns
#### - Species-specific pattern
####

# Load workspace
load("../08_EnvironmentNetwork_spOut/08_EnvironmentNetwork_spOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.17
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(mgcv); packageVersion("mgcv") # 1.8.40, 2022.5.17
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"

# Prepare color palette
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
get_palette2 <- colorRampPalette(brewer.pal(8, "Set1"))
palette_custom <- rev(get_palette(11))
palette_custom3 <- get_palette2(28)


# <---------------------------------------------> #
#  Add information
# <---------------------------------------------> #
tax_tbl$scientific_name
tax_tbl$common_jp_name

sci_name_id1 <- match(uic_same_dup_long2$cause_var2, tax_tbl$tax_id)
sci_name_id2 <- match(uic_same_dup_long3$cause_var2, tax_tbl$tax_id)
uic_same_dup_long2$fish_sci <- factor(as.character(tax_tbl$scientific_name[sci_name_id1]))
uic_same_dup_long3$fish_sci <- factor(as.character(tax_tbl$scientific_name[sci_name_id2]))


# <---------------------------------------------> #
#  Visualization: Overview
# <---------------------------------------------> #
## Effects of temperature on fish species interactions
t1 <- uic_same_dup_long2 %>%
  ggplot(aes(x = value, y = te, color = env_vars), color = "gray60") +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray80", method = "loess", size = 1) +
  facet_wrap(.~ env_vars, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette2(12), name = "Abiotic/biotic variables") +
  scale_y_log10() + panel_border() +
  ylab("Total median TE") +
  NULL


# <---------------------------------------------> #
#  Visualization: GLM
# <---------------------------------------------> #
## Compile GLM results
for(i in 1:50) {
  pred_all[[i]]$tax_id <- names(pred_all)[i]
  pred_all[[i]]$fish_sci <- as.character(tax_tbl$scientific_name[which(tax_tbl$tax_id == names(pred_all)[i])])
}

# Extract fish species with significant positive and negative slopes
df_all$tax_id <- rownames(df_all)
sci_name_id <- match(rownames(df_all) , tax_tbl$tax_id)
df_all$fish_sci <- as.character(tax_tbl$scientific_name[sci_name_id])
df_all$signif <- "N.S."
df_all$signif[df_all$p.value < 0.05] <- "P < 0.05"
df_all$signif[df_all$p.value >= 0.05 & df_all$p.value <= 0.10 ] <- "P < 0.10"
signif_tax <- df_all[df_all$signif == "P < 0.05", "tax_id"]
signif_sci <- as.character(df_all[df_all$signif == "P < 0.05", "fish_sci"])
signif_id <- match(signif_tax, names(pred_all))

pred_all_lm <- pred_all[[signif_id[1]]]
pred_all_lm$slope <- if(df_all[signif_id[1],"Estimate"] >= 0) "positive" else "negative"
for (i in signif_id[2:length(signif_id)]) {
  pred_all_lm_i <- pred_all[[i]]
  pred_all_lm_i$slope <- if(df_all[i,"Estimate"] >= 0) "positive" else "negative"
  pred_all_lm <- rbind(pred_all_lm, pred_all_lm_i)
}

uic_same_dup_long4 <- uic_same_dup_long3 %>% filter(fish_sci %in% signif_sci)

t2 <- uic_same_dup_long4 %>%
  ggplot(aes(x = value, y = te), color = "gray60") +
  geom_line(data = pred_all_lm, aes(x = value, y = pred, color = slope), size = 1.5) +
  scale_color_manual(values = c("royalblue", "red3")) +
  geom_point(alpha = 0.1) +
  facet_wrap(.~ fish_sci, ncol = 5, scales = "free") +
  scale_y_log10() + panel_border() +
  xlab(expression(paste("Median water temperature (", degree, "C)"))) +
  ylab("Interaction strength (TE)") +
  theme(legend.position = "none", strip.text = element_text(face = "italic", size = 10)) +
  NULL


# <---------------------------------------------> #
#  Visualization: GAM
# <---------------------------------------------> #
## Compile GAM results
for(i in 1:50) {
  gam_pred_all[[i]]$tax_id <- names(gam_pred_all)[i]
  gam_pred_all[[i]]$fish_sci <- as.character(tax_tbl$scientific_name[which(tax_tbl$tax_id == names(gam_pred_all)[i])])
}

#edf   Ref.df           F      p.value
#Taxa0479 7.443692 8.255424  3.52169559 0.0005255407

# Extract fish species with significant positive and negative slopes
gam_df_all$tax_id <- rownames(gam_df_all)
gam_sci_name_id <- match(rownames(gam_df_all) , tax_tbl$tax_id)
gam_df_all$fish_sci <- as.character(tax_tbl$scientific_name[gam_sci_name_id])
gam_df_all$signif <- "N.S."
gam_df_all$signif[gam_df_all$p.value < 0.05] <- "P < 0.05"
gam_df_all$signif[gam_df_all$p.value >= 0.05 & gam_df_all$p.value <= 0.10 ] <- "P < 0.10"
gam_signif_tax <- gam_df_all[gam_df_all$signif == "P < 0.05", "tax_id"]
gam_signif_sci <- as.character(gam_df_all[gam_df_all$signif == "P < 0.05", "fish_sci"])
gam_signif_id <- match(gam_signif_tax, names(gam_pred_all))
#length(gam_signif_id)
pred_all_gamcomb <- gam_pred_all[[gam_signif_id[1]]]
## Compare the initial and last values
pred_all_gamcomb$slope <- if((rev(gam_pred_all[[1]]$pred)[1] - gam_pred_all[[1]]$pred[1]) >= 0) "positive" else "negative"
for (i in gam_signif_id[2:length(gam_signif_id)]) {
  pred_all_gamcomb_i <- gam_pred_all[[i]]
  pred_all_gamcomb_i$slope <- if((rev(gam_pred_all[[i]]$pred)[1] - gam_pred_all[[i]]$pred[1]) >= 0) "positive" else "negative"
  pred_all_gamcomb <- rbind(pred_all_gamcomb, pred_all_gamcomb_i)
}
## Add the difference between the initial and last values for histogram
gam_df_all$Estimate_diff <- NA
gam_df_all$slope <- NA
for (i in 1:nrow(gam_df_all)) {
  gam_df_all$Estimate_diff[i] <- (rev(gam_pred_all[[i]]$pred)[1] - gam_pred_all[[i]]$pred[1])
  gam_df_all$slope <- if (gam_df_all$Estimate_diff[i] >= 0) "positive" else "negative"
}

# Select GAM significant fish species
uic_same_dup_long5 <- uic_same_dup_long3 %>% filter(fish_sci %in% gam_signif_sci)

t3 <- uic_same_dup_long5 %>%
  ggplot(aes(x = value, y = te), color = "gray60") +
  geom_line(data = pred_all_gamcomb, aes(x = value, y = pred, color = slope), size = 1.5) +
  scale_color_manual(values = c("royalblue", "red3")) +
  geom_point(alpha = 0.1) +
  facet_wrap(.~ fish_sci, ncol = 8, scales = "free") +
  scale_y_log10() + panel_border() +
  xlab(expression(paste("Median water temperature (", degree, "C)"))) +
  ylab("Interaction strength (TE)") +
  theme(legend.position = "none", strip.text = element_text(face = "italic", size = 8)) +
  NULL

## Quick GAM version
# if(F){
#   t3_0 <- uic_same_dup_long5 %>%
#     ggplot(aes(x = value, y = te), color = "gray60") +
#     stat_smooth(color = "gray10", method = "gam", se = F, size = 1) +
#     scale_color_manual(values = c("royalblue", "red3")) +
#     geom_point(alpha = 0.1) +
#     facet_wrap(.~ fish_sci, ncol = 8, scales = "free") +
#     scale_y_log10() + panel_border() +
#     xlab(expression(paste("Median water temperature (", degree, "C)"))) +
#     ylab("Interaction strength (TE)") +
#     theme(legend.position = "none", strip.text = element_text(face = "italic", size = 10)) +
#     NULL
# }


## Histogram
### GLM
h1 <- ggplot(df_all, aes(x = Estimate, fill = signif)) +
  scale_fill_manual(values = c("gray30", "red3", "lightsalmon1"),
                    name = NULL) +
  geom_histogram(breaks = seq(-0.6, 0.6, by = 0.075), alpha = 0.8,
                 color = "white", size = 0.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(breaks = seq(1,13, by = 1)) +
  xlab("Temperature dependence of interactions\n(regression slope)") +
  ylab("Count") + theme(legend.position = c(0.8,0.8)) +
  NULL

### GAM
h2 <- ggplot(gam_df_all, aes(x = Estimate_diff, fill = signif)) +
  scale_fill_manual(values = c("gray30", "red3", "lightsalmon1"),
                    name = NULL) +
  geom_histogram(breaks = seq(-0.225, 0.225, by = 0.025), alpha = 0.8,
                 color = "white", size = 0.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(breaks = seq(1,13, by = 1)) +
  xlab("Temperature dependence of interactions") +
  ylab("Count") + theme(legend.position = c(0.8,0.8)) +
  NULL
h2


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save figures
## Effect TE
ggsave(file = sprintf("%s/PDF_TotalMedTE_env.pdf", fig_folder1),
       plot = t1, width = 12, height = 10)
ggsave(file = sprintf("%s/PDF_TotalMedTE_temp.pdf", fig_folder1),
       plot = t2, width = 12, height = 10)
ggsave(file = sprintf("%s/PDF_TotalMedTE_temp_GAM.pdf", fig_folder1),
       plot = t3, width = 18, height = 14)

ggsave(file = sprintf("%s/PDF_Histogram.pdf", fig_folder1),
       plot = h1, width = 6, height = 4)
ggsave(file = sprintf("%s/PDF_Histogram_GAM.pdf", fig_folder1),
       plot = h2, width = 6, height = 4)

# Save objects
saveRDS(list(t1, t2, h1), sprintf("%s/Fig_TotalMedTE_env.obj", fig_folder1))
saveRDS(list(t3, h2), sprintf("%s/Fig_TotalMedTE_env_GAM.obj", fig_folder1))
saveRDS(uic_same_dup_long3, sprintf("%s/FigData_TotalMedTE_env.obj", fig_folder1))

