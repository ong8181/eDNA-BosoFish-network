####
#### Boso Peninsula project
#### No. 9 Statistical test
#### 2022.12.12 Ushio (R4.2.1)
####

# Load workspace
load("06_CompileSmapCoefOut/06_CompileSmapCoefOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
library(mgcv); packageVersion("mgcv") # 1.8.40
theme_set(theme_cowplot())

# Custom library
library(macam); packageVersion("macam") # 0.0.10

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
(output_folder <- outdir_create())


# <---------------------------------------------> #
# Check effects of temperature on S-map coefficients 
# <---------------------------------------------> #
# ------------------------------------------------------- #
# Water temperature effects (in-strength + out-strength)
# ------------------------------------------------------- #
## Generate a summary object
water_temp_gamm1 <- water_temp_gamm2 <- data.frame(fish_taxa = top_taxa,
                                                   lme_value = NA, lme_DF = NA, lme_tval = NA, lme_pval = NA,
                                                   gam_edf = NA, gam_F = NA, gam_pval = NA)
## Check the dependence of interaction strength for top 50 fish species
for (i in 1:length(top_taxa)) {
  # Perform GAMM (in-strength)
  gamm_res <- sdf %>% filter(effect_var == top_taxa[i]) %>% 
    gamm(abs(IS) ~ s(water_temp), data = ., family = Gamma(link = "log"),
         random=list(site_code=~1), niterPQL = 100)
  # Assign values to data.frame
  water_temp_gamm1[i,c("gam_edf", "gam_F", "gam_pval")] <- summary(gamm_res$gam)$s.table[c(1,3,4)]
  water_temp_gamm1[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <- summary(gamm_res$lme)$tTable[2,c(1,3,4,5)]
  
  # Perform GAMM (out-strength)
  gamm_res <- sdf %>% filter(cause_var == top_taxa[i]) %>% 
    gamm(abs(IS) ~ s(cause_temp), data = ., family = Gamma(link = "log"),
         random=list(site_code=~1), niterPQL = 100)
  # Assign values to data.frame
  water_temp_gamm2[i,c("gam_edf", "gam_F", "gam_pval")] <- summary(gamm_res$gam)$s.table[c(1,3,4)]
  water_temp_gamm2[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <- summary(gamm_res$lme)$tTable[2,c(1,3,4,5)]
}
# Assign significance
water_temp_gamm1$lme_sig <- (water_temp_gamm1$lme_pval < 0.05)
water_temp_gamm1$gam_sig <- (water_temp_gamm1$gam_pval < 0.05)
sum(water_temp_gamm1$lme_sig); sum(water_temp_gamm1$gam_sig)
water_temp_gamm2$lme_sig <- (water_temp_gamm2$lme_pval < 0.05)
water_temp_gamm2$gam_sig <- (water_temp_gamm2$gam_pval < 0.05)
sum(water_temp_gamm2$lme_sig); sum(water_temp_gamm2$gam_sig)


# ------------------------------------------------------- #
# Total eDNA effects (in-strength + out-strength)
# ------------------------------------------------------- #
# Generate a summary object
total_edna_gamm1 <- total_edna_gamm2 <- data.frame(fish_taxa = top_taxa,
                                                   lme_value = NA, lme_DF = NA, lme_tval = NA, lme_pval = NA,
                                                   gam_edf = NA, gam_F = NA, gam_pval = NA)
# Check the dependence of interaction strength for top 50 fish species
for (i in 1:length(top_taxa)) {
  # Perform GAMM (in-strength)
  gamm_res <- sdf %>% filter(effect_var == top_taxa[i]) %>% 
    gamm(abs(IS) ~ s(total_dna_per_ml), data = ., family = Gamma(link = "log"), random=list(site_code=~1),
         niterPQL = 100)
  # Assign values to data.frame
  total_edna_gamm1[i,c("gam_edf", "gam_F", "gam_pval")] <- summary(gamm_res$gam)$s.table[c(1,3,4)]
  total_edna_gamm1[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <- summary(gamm_res$lme)$tTable[2,c(1,3,4,5)]
  
  # Perform GAMM (out-strength)
  gamm_res <- try(sdf %>% filter(cause_var == top_taxa[i]) %>% 
                    gamm(abs(IS) ~ s(cause_totdna), data = ., family = Gamma(link = "log"), random=list(site_code=~1),
                         niterPQL = 100))
  # Assign values to data.frame
  if(class(gamm_res)[[1]] != "try-error") {
    total_edna_gamm2[i,c("gam_edf", "gam_F", "gam_pval")] <- summary(gamm_res$gam)$s.table[c(1,3,4)]
    total_edna_gamm2[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <- summary(gamm_res$lme)$tTable[2,c(1,3,4,5)]
  } else {
    total_edna_gamm2[i,c("gam_edf", "gam_F", "gam_pval")] <- c(NA, NA, NA)
    total_edna_gamm2[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <-  c(NA, NA, NA, NA)
  }
}
# Assign significance
total_edna_gamm1$lme_sig <- (total_edna_gamm1$lme_pval < 0.05)
total_edna_gamm1$gam_sig <- (total_edna_gamm1$gam_pval < 0.05)
sum(total_edna_gamm1$lme_sig); sum(total_edna_gamm1$gam_sig)
total_edna_gamm2$lme_sig <- (total_edna_gamm2$lme_pval < 0.05)
total_edna_gamm2$gam_sig <- (total_edna_gamm2$gam_pval < 0.05)
sum(total_edna_gamm2$lme_sig, na.rm = T); sum(total_edna_gamm2$gam_sig, na.rm = T)


# ------------------------------------------------------- #
# Species richness effects (in-strength + out-strength)
# ------------------------------------------------------- #
richness_gamm1 <- richness_gamm2 <- data.frame(fish_taxa = top_taxa,
                                               lme_value = NA, lme_DF = NA, lme_tval = NA, lme_pval = NA,
                                               gam_edf = NA, gam_F = NA, gam_pval = NA)
# Check the dependence of interaction strength for top 50 fish species
for (i in 1:length(top_taxa)) {
  # Perform GAMM (in-strength)
  gamm_res <- sdf %>% filter(effect_var == top_taxa[i]) %>% 
    gamm(abs(IS) ~ s(sp_richness), data = ., family = Gamma(link = "log"), random=list(site_code=~1),
         niterPQL = 100)
  # Assign values to data.frame
  richness_gamm1[i,c("gam_edf", "gam_F", "gam_pval")] <- summary(gamm_res$gam)$s.table[c(1,3,4)]
  richness_gamm1[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <- summary(gamm_res$lme)$tTable[2,c(1,3,4,5)]
  
  # Perform GAMM (out-strength)
  gamm_res <- sdf %>% filter(cause_var == top_taxa[i]) %>% 
    gamm(abs(IS) ~ s(cause_richness), data = ., family = Gamma(link = "log"), random=list(site_code=~1),
         niterPQL = 100)
  # Assign values to data.frame
  richness_gamm2[i,c("gam_edf", "gam_F", "gam_pval")] <- summary(gamm_res$gam)$s.table[c(1,3,4)]
  richness_gamm2[i,c("lme_value", "lme_DF", "lme_tval", "lme_pval")] <- summary(gamm_res$lme)$tTable[2,c(1,3,4,5)]
}
# Assign significance
richness_gamm1$lme_sig <- (richness_gamm1$lme_pval < 0.05)
richness_gamm1$gam_sig <- (richness_gamm1$gam_pval < 0.05)
sum(richness_gamm1$lme_sig); sum(richness_gamm1$gam_sig)
# Assign significance
richness_gamm2$lme_sig <- (richness_gamm2$lme_pval < 0.05)
richness_gamm2$gam_sig <- (richness_gamm2$gam_pval < 0.05)
sum(richness_gamm2$lme_sig); sum(richness_gamm2$gam_sig)


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))
# Save session info
macam::save_session_info(session_info_dir = "00_SessionInfo")
