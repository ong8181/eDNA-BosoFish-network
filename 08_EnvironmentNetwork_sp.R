####
#### Boso Peninsula project
#### No. 8 Environmental variables and network patterns
#### - Species-specific pattern
####

# Load workspace
load("07_EnvironmentNetworkOut/07_EnvironmentNetworkOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.5.17
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(mgcv); packageVersion("mgcv") # 1.8.40, 2022.5.17
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

# Set random seeds (for reproduction)
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)


# <---------------------------------------------> #
#  Compile data
# <---------------------------------------------> #
## Starting objects
uic_same_st
## Add site information
site_id_match <- match(uic_same_st$site, site_summary_df$site_code)
uic_same_st2 <- site_summary_df[site_id_match,] %>% 
  bind_cols(
    uic_same_st %>% select(effect_var, cause_var, tp, rmse, te,
                           cause_pair, effect_var2, cause_var2)
    )
## Switch effect and causal variable
uic_same_st3 <- uic_same_st2
uic_same_st3$effect_var <- uic_same_st2$cause_var
uic_same_st3$cause_var <- uic_same_st2$effect_var
uic_same_st3$effect_var2 <- uic_same_st2$cause_var2
uic_same_st3$cause_var2 <- uic_same_st2$effect_var2
uic_same_st2$switched <- "N"
uic_same_st3$switched <- "Y"
## Bind rows
uic_same_st_dupli <- rbind(uic_same_st2, uic_same_st3)
## pivot_longer
uic_same_dup_long <- uic_same_st_dupli %>% 
  select(site_code, cause_pair,
         cause_var, effect_var, cause_var2, effect_var2,
         switched, te,
         temp_mean, temp_med, temp_min, temp_max,
         lat, linear_dist, salinity_mean,
         wave_m_mean, tide_mean, richness,
         total_edna) %>%
  pivot_longer(cols = -c(site_code, cause_pair,
                         cause_var, effect_var, cause_var2, effect_var2,
                         switched, te),
               names_to = "env_vars",
               values_to = "value")


# ------------------------------------------- #
# Select dominant species
# ------------------------------------------- #
taxa_orders_conc <- names(sort(colSums(asv_df_conv), decreasing = T))
taxa_orders_freq <- names(sort(colSums(asv_df_conv > 0), decreasing = T))


# --------------------------------------------- #
#  Visualize pattern
# --------------------------------------------- #
uic_same_dup_long2 <- uic_same_dup_long %>% filter(cause_var2 %in% taxa_orders_freq[1:12])

## Example
t1 <- uic_same_dup_long2 %>% filter(env_vars == "temp_med") %>%
  ggplot(aes(x = value, y = te, color = cause_var2)) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap(. ~ cause_var2, ncol = 3) +
  scale_y_log10() +
  panel_border()


# --------------------------------------------- #
# General additive model (GAM)
# --------------------------------------------- #
uic_same_dup_long3 <- uic_same_dup_long %>% filter(env_vars == "temp_med")

# Fish species-specific GAMM [TE v.s. temperature]
# Some fish species do not have sufficient data point ==> use glm
gam_df_all <- data.frame("edf" = rep(NaN, 50),
                          "Ref.df" = rep(NaN, 50),
                          "F" = rep(NaN, 50),
                          "p-value" = rep(NaN, 50))
gam_pred_all <- list(NULL)
all(unique(uic_same_dup_long3$cause_var2) %in% taxa_orders_freq[1:50])
for(sp_i in 1:50){
  if (sp_i == 23 | sp_i == 40 | sp_i == 47 | sp_i == 49) {
    gam_i <- gam(te ~ s(value, k = 9), family = Gamma(link = "log"),
                   data = uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]))
  } else {
    gam_i <- gam(te ~ s(value), family = Gamma(link = "log"),
                   data = uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]))
  }
  gam_df_all[sp_i,] <- as.data.frame(summary(gam_i)$s.table)
  rownames(gam_df_all)[sp_i] <- taxa_orders_freq[sp_i]
  x_min_i <- min(uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]) %>% .$value)
  x_max_i <- max(uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]) %>% .$value)
  new_d <- data.frame(value = seq(x_min_i, x_max_i, length.out = 100))
  pred_i <- exp(predict.gam(gam_i, newdata = new_d))
  gam_pred_all[[sp_i]] <- data.frame(value = new_d$value, pred = pred_i)
  names(gam_pred_all)[sp_i] <- taxa_orders_freq[sp_i]
}

# Some fish species do not have sufficient data point ==> use glm
df_all <- data.frame("Estimate" = rep(NaN, 50),
                     "Std. Error" = rep(NaN, 50),
                     "t value" = rep(NaN, 50),
                     "p-value" = rep(NaN, 50))
pred_all <- list(NULL)
for(sp_i in 1:50){
  glm_i <- glm(te ~ value, family = Gamma(link = "log"),
               data = uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]))
  df_all[sp_i,] <- matrix(summary(glm_i)$coefficients[2,], ncol = 4)
  rownames(df_all)[sp_i] <- taxa_orders_freq[sp_i]
  x_min_i <- min(uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]) %>% .$value)
  x_max_i <- max(uic_same_dup_long3 %>% filter(cause_var2 == taxa_orders_freq[sp_i]) %>% .$value)
  new_d <- data.frame(value = seq(x_min_i, x_max_i, length.out = 100))
  pred_i <- exp(predict(glm_i, newdata = new_d))
  pred_all[[sp_i]] <- data.frame(value = new_d$value, pred = pred_i)
  names(pred_all)[sp_i] <- taxa_orders_freq[sp_i]
}


# --------------------------------------------- #
# Save results
# --------------------------------------------- #
# Save figures
## Total TE
quartz(file = sprintf("%s/total_medTE_temp.pdf", output_folder),
       type = "pdf", family = "HiraKakuPro-W3", width = 12, height = 10)
t1; dev.off()

# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))

