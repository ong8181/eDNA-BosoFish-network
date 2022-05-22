####
#### Boso Peninsula project
#### No. 7 Environmental variables and network patterns
####

# Load workspace
load("06_VisualizeNetwork2Out/06_VisualizeNetwork2Out.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.8.25
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(mgcv); packageVersion("mgcv") # 1.8.36, 2021.8.26
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)


# <---------------------------------------------> #
#  Check UIC results
# <---------------------------------------------> #
#uic_edna_env_strong
#uic_env_edna_strong
#uic_edna_sp_strong
#uic_edna_st_strong


# <---------------------------------------------> #
#  Brief check
# <---------------------------------------------> #
## Extract site-specific UICs
effect_site <- str_sub(uic_edna_st_strong$effect_var, start = 10, end = -1)
cause_site <- str_sub(uic_edna_st_strong$cause_var, start = 10, end = -1)
same_site <- effect_site == cause_site
uic_same_st <- uic_edna_st_strong[same_site,]
uic_same_st$site <- str_sub(uic_same_st$effect_var, start = 10, end = -1)
uic_same_st <- uic_same_st %>%
  mutate(effect_var2 = str_sub(effect_var, start = 1, end = 8)) %>%
  mutate(cause_var2 = str_sub(cause_var, start = 1, end = 8))

## Calculate diversity
asv_site_df <- aggregate(asv_df_conv, by = list(sample_df$site_code), sum) %>% .[,-1]
dim(asv_site_df)
site_diversity <- rowSums(asv_site_df > 0)
asv_quantile <- quantile(as.numeric(unlist(c(asv_site_df))), na.rm = T, probs = seq(0,1,0.1))
site_diversity_dom <- rowSums(asv_site_df > asv_quantile["90%"])
site_abundance <- rowSums(asv_site_df, na.rm = T) # = total DNA copy / ml for each site


## Site-grouped summary
site_summary_df <- sample_df %>%
  group_by(site_code) %>%
  summarize(temp_mean = mean(water_temp, na.rm = T),
            temp_med = median(water_temp, na.rm = T),
            temp_min = min(water_temp, na.rm = T),
            temp_max = max(water_temp, na.rm = T),
            lat = mean(lat_n, na.rm = T),
            lon = mean(long_e, na.rm = T),
            linear_dist = mean(linear_dist, na.rm = T),
            salinity_mean = mean(salinity, na.rm = T),
            wave_m_mean = mean(wave_m, na.rm = T),
            tide_mean = mean(tide, na.rm = T)) %>%
  mutate(richness = site_diversity, 
         richness_dom = site_diversity_dom,
         total_edna = site_abundance) # = total DNA copy / ml for each site
saveRDS(site_summary_df, sprintf("%s/site_summary_df.obj", output_folder))

## Collect information
stat_te_effect <- uic_same_st %>% group_by(site, effect_var2) %>%
  summarize(te_med = median(te), te_min = min(te), te_max = max(te), te_mean = mean(te), n_int = n())
stat_te_cause <- uic_same_st %>% group_by(site, cause_var2) %>%
  summarize(te_med = median(te), te_min = min(te), te_max = max(te), te_mean = mean(te), n_int = n())
stat_te_effect$site_sp <- paste0(stat_te_effect$site, "_", stat_te_effect$effect_var2)
stat_te_cause$site_sp <- paste0(stat_te_cause$site, "_", stat_te_cause$cause_var2)

## Assign site parameters to te_df
### Effect var
site_position_id1 <- match(stat_te_effect$site, site_summary_df$site_code)
effect_te_df <- data.frame(site = stat_te_effect$site,
                           site_sp = paste0(stat_te_effect$site, "_", stat_te_effect$effect_var2),
                           lat = site_summary_df[site_position_id1, "lat"],
                           lon = site_summary_df[site_position_id1, "lon"],
                           temp_mean = site_summary_df[site_position_id1, "temp_mean"],
                           temp_med = site_summary_df[site_position_id1, "temp_med"],
                           temp_min = site_summary_df[site_position_id1, "temp_min"],
                           temp_max = site_summary_df[site_position_id1, "temp_max"],
                           linear_dist = site_summary_df[site_position_id1, "linear_dist"],
                           salinity_mean = site_summary_df[site_position_id1, "salinity_mean"],
                           wave_m_mean = site_summary_df[site_position_id1, "wave_m_mean"],
                           tide_mean = site_summary_df[site_position_id1, "tide_mean"],
                           richness = site_summary_df[site_position_id1, "richness"],
                           richness_dom = site_summary_df[site_position_id1, "richness_dom"],
                           total_edna = site_summary_df[site_position_id1, "total_edna"],
                           sp = stat_te_effect$effect_var2,
                           n_int = stat_te_effect$n_int,
                           te_med = stat_te_effect$te_med,
                           te_mean = stat_te_effect$te_mean,
                           te_min = stat_te_effect$te_min,
                           te_max = stat_te_effect$te_max)
effect_te_long <- effect_te_df %>%
  select(site_sp, lat, linear_dist,
         temp_med, temp_mean, temp_min, temp_max,
         salinity_mean, wave_m_mean, tide_mean,
         richness, richness_dom, total_edna,
         n_int,
         te_med, te_mean, te_min, te_max) %>%
  pivot_longer(cols = -c(site_sp, te_med, te_mean, te_min, te_max))


### Cause var
site_position_id2 <- match(stat_te_cause$site, site_summary_df$site_code)
cause_te_df <- data.frame(site = stat_te_cause$site,
                          site_sp = paste0(stat_te_cause$site, "_", stat_te_cause$cause_var2),
                          lat = site_summary_df[site_position_id2, "lat"],
                          lon = site_summary_df[site_position_id2, "lon"],
                          temp_mean = site_summary_df[site_position_id2, "temp_mean"],
                          temp_med = site_summary_df[site_position_id2, "temp_med"],
                          temp_min = site_summary_df[site_position_id2, "temp_min"],
                          temp_max = site_summary_df[site_position_id2, "temp_max"],
                          linear_dist = site_summary_df[site_position_id2, "linear_dist"],
                          salinity_mean = site_summary_df[site_position_id2, "salinity_mean"],
                          wave_m_mean = site_summary_df[site_position_id2, "wave_m_mean"],
                          tide_mean = site_summary_df[site_position_id2, "tide_mean"],
                          richness = site_summary_df[site_position_id2, "richness"],
                          richness_dom = site_summary_df[site_position_id2, "richness_dom"],
                          total_edna = site_summary_df[site_position_id2, "total_edna"],
                          sp = stat_te_cause$cause_var2,
                          n_int = stat_te_cause$n_int,
                          te_med = stat_te_cause$te_med,
                          te_mean = stat_te_cause$te_mean,
                          te_min = stat_te_cause$te_min,
                          te_max = stat_te_cause$te_max)
cause_te_long <- cause_te_df %>%
  select(site_sp, lat, linear_dist,
         temp_med, temp_mean, temp_min, temp_max,
         salinity_mean, wave_m_mean, tide_mean,
         richness, richness_dom, total_edna,
         n_int,
         te_med, te_mean, te_min, te_max) %>%
  pivot_longer(cols = -c(site_sp, te_med, te_mean, te_min, te_max))


### TE sums up
unique_site_sp <- sort(unique(c(cause_te_df$site_sp, effect_te_df$site_sp)))
effect_site_sp_id <- match(effect_te_df$site_sp, unique_site_sp)
cause_site_sp_id <- match(cause_te_df$site_sp, unique_site_sp)
uic_same_st_cause_effect_rev <- uic_same_st
uic_same_st_cause_effect_rev$effect_var2 <- uic_same_st$cause_var2
uic_same_st_cause_effect_rev$cause_var2 <- uic_same_st$effect_var2
uic_same_duplicate <- rbind(uic_same_st, uic_same_st_cause_effect_rev)

te_dupli_sum <- uic_same_duplicate %>%
  group_by(site, effect_var2) %>%
  summarize(te_med = median(te), te_mean = mean(te), te_min = min(te), te_max = max(te), n_int = n())
te_dupli_sum$site_sp <-  paste0(te_dupli_sum$site, "_", te_dupli_sum$effect_var2)

site_position_id3 <- match(te_dupli_sum$site, site_summary_df$site_code)
site_position_id4 <- match(te_dupli_sum$site_sp, stat_te_effect$site_sp)
site_position_id5 <- match(te_dupli_sum$site_sp, stat_te_cause$site_sp)
sum_te_df <- data.frame(site = te_dupli_sum$site,
                        site_sp = paste0(te_dupli_sum$site, "_", te_dupli_sum$effect_var2),
                        lat = site_summary_df[site_position_id3, "lat"],
                        lon = site_summary_df[site_position_id3, "lon"],
                        temp_mean = site_summary_df[site_position_id3, "temp_mean"],
                        temp_med = site_summary_df[site_position_id3, "temp_med"],
                        temp_min = site_summary_df[site_position_id3, "temp_min"],
                        temp_max = site_summary_df[site_position_id3, "temp_max"],
                        linear_dist = site_summary_df[site_position_id3, "linear_dist"],
                        salinity_mean = site_summary_df[site_position_id3, "salinity_mean"],
                        wave_m_mean = site_summary_df[site_position_id3, "wave_m_mean"],
                        tide_mean = site_summary_df[site_position_id3, "tide_mean"],
                        richness = site_summary_df[site_position_id3, "richness"],
                        richness_dom = site_summary_df[site_position_id3, "richness_dom"],
                        total_edna = site_summary_df[site_position_id3, "total_edna"],
                        sp = te_dupli_sum$effect_var2,
                        n_int = te_dupli_sum$n_int,
                        n_int_effect = stat_te_effect$n_int[site_position_id4],
                        n_int_cause = stat_te_cause$n_int[site_position_id5],
                        te_med = te_dupli_sum$te_med,
                        te_mean = te_dupli_sum$te_mean,
                        te_min = te_dupli_sum$te_min,
                        te_max = te_dupli_sum$te_max)
sum_te_long <- sum_te_df %>%
  select(site_sp, lat, linear_dist,
         temp_med, temp_mean, temp_min, temp_max,
         salinity_mean, wave_m_mean, tide_mean,
         richness, richness_dom, total_edna,
         n_int, n_int_effect, n_int_cause,
         te_med, te_mean, te_min, te_max) %>%
  pivot_longer(cols = -c(site_sp, te_med, te_mean, te_min, te_max))




# <---------------------------------------------> #
#  Visualization
# <---------------------------------------------> #
## Effect TE patterns
e1 <- ggplot(effect_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_min, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Effect TE") +
  panel_border() + ggtitle("min TE") +
  NULL
e2 <- ggplot(effect_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_med, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Effect TE") +
  panel_border() + ggtitle("Median TE") +
  NULL
e3 <- ggplot(effect_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_mean, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Effect TE") +
  panel_border() + ggtitle("Mean TE") +
  NULL
e4 <- ggplot(effect_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_max, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Effect TE") +
  panel_border() + ggtitle("Max TE") +
  NULL

## Causal TE patterns
c1 <- ggplot(cause_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_min, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Causal TE") +
  panel_border() + ggtitle("min TE") +
  NULL
c2 <- ggplot(cause_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_med, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Causal TE") +
  panel_border() + ggtitle("Median TE") +
  NULL
c3 <- ggplot(cause_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_mean, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Causal TE") +
  panel_border() + ggtitle("Mean TE") +
  NULL
c4 <- ggplot(cause_te_long %>% filter(name != "n_int"),
             aes(x = value, y = te_max, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Causal TE") +
  panel_border() + ggtitle("Max TE") +
  NULL


## total TE patterns
t1 <- ggplot(sum_te_long %>% filter(name != "n_int" & name != "n_int_effect" & name != "n_int_cause"),
             aes(x = value, y = te_min, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Total TE") +
  panel_border() + ggtitle("min TE") +
  NULL
t2 <- ggplot(sum_te_long  %>% filter(name != "n_int" & name != "n_int_effect" & name != "n_int_cause"),
             aes(x = value, y = te_med, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Total TE") +
  panel_border() + ggtitle("Median TE") +
  NULL
t3 <- ggplot(sum_te_long  %>% filter(name != "n_int" & name != "n_int_effect" & name != "n_int_cause"),
             aes(x = value, y = te_mean, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Total TE") +
  panel_border() + ggtitle("Mean TE") +
  NULL
t4 <- ggplot(sum_te_long %>% filter(name != "n_int" & name != "n_int_effect" & name != "n_int_cause"),
             aes(x = value, y = te_max, color = name)) +
  geom_point(alpha = 0.5) +
  stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  facet_wrap(.~ name, ncol = 3, scales = "free") +
  scale_color_manual(values = get_palette(12)) +
  scale_y_log10() +
  xlab("Environmental variables") + ylab("Total TE") +
  panel_border() + ggtitle("Max TE") +
  NULL


# <---------------------------------------------> #
#                    Save results                 #
# <---------------------------------------------> #
# Save figures
## Effect TE
ggsave(sprintf("%s/effect_minTE.pdf", output_folder), plot = e1, width = 12, height = 10)
ggsave(sprintf("%s/effect_medTE.pdf", output_folder), plot = e2, width = 12, height = 10)
ggsave(sprintf("%s/effect_meanTE.pdf", output_folder), plot = e3, width = 12, height = 10)
ggsave(sprintf("%s/effect_maxTE.pdf", output_folder), plot = e4, width = 12, height = 10)
## Cause TE
ggsave(sprintf("%s/cause_minTE.pdf", output_folder), plot = c1, width = 12, height = 10)
ggsave(sprintf("%s/cause_medTE.pdf", output_folder), plot = c2, width = 12, height = 10)
ggsave(sprintf("%s/cause_meanTE.pdf", output_folder), plot = c3, width = 12, height = 10)
ggsave(sprintf("%s/cause_maxTE.pdf", output_folder), plot = c4, width = 12, height = 10)
## Total TE
ggsave(sprintf("%s/total_minTE.pdf", output_folder), plot = t1, width = 12, height = 10)
ggsave(sprintf("%s/total_medTE.pdf", output_folder), plot = t2, width = 12, height = 10)
ggsave(sprintf("%s/total_meanTE.pdf", output_folder), plot = t3, width = 12, height = 10)
ggsave(sprintf("%s/total_maxTE.pdf", output_folder), plot = t4, width = 12, height = 10)

# Save results
# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))
