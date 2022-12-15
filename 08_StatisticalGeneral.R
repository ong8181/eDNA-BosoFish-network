####
#### Boso Peninsula project
#### No. 8 Statistics for general pattern
#### 2022.12.12 Ushio (R4.2.1)
####

# Load workspace
load("06_CompileSmapCoefOut/06_CompileSmapCoefOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
library(mgcv); packageVersion("mgcv") # 1.8.40
library(gamm4); packageVersion("gamm4") # 0.2.6
theme_set(theme_cowplot())

# Custom library
library(macam); packageVersion("macam") # 0.0.10

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
(output_folder <- outdir_create())


# <---------------------------------------------> #
# Testing overall pattern
# <---------------------------------------------> #
class(sdf$IS); class(sdf$water_temp); class(sdf$sp_richness); class(sdf$total_dna_per_ml) 
class(sdf$salinity); class(sdf$tide); class(sdf$wave_m)
sdf$effect_var <- factor(sdf$effect_var)
sdf$cause_var <- factor(sdf$cause_var)

# Water temperature, species richness, total eDNA concentrations
gamm01 <- sdf %>% gamm(abs(IS) ~ s(water_temp), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, effect_var = ~1), niterPQL = 20)
gamm02 <- sdf %>% gamm(abs(IS) ~ s(sp_richness), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, effect_var = ~1), niterPQL = 20)
gamm03 <- sdf %>% gamm(abs(IS) ~ s(total_dna_per_ml), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, effect_var = ~1), niterPQL = 20)
gamm04 <- sdf %>% gamm(abs(IS) ~ s(cause_temp), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, cause_var = ~1), niterPQL = 20)
gamm05 <- sdf %>% gamm(abs(IS) ~ s(cause_richness), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, cause_var = ~1), niterPQL = 20)
gamm06 <- sdf %>% gamm(abs(IS) ~ s(cause_totdna), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, cause_var = ~1), niterPQL = 20)
# Salinity, tide and wave
gamm07 <- sdf %>% gamm(abs(IS) ~ s(salinity), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, effect_var = ~1), niterPQL = 20)
gamm08 <- sdf %>% gamm(abs(IS) ~ s(tide), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, effect_var = ~1), niterPQL = 20)
gamm09 <- sdf %>% gamm(abs(IS) ~ s(wave_m), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, effect_var = ~1), niterPQL = 20)
gamm10 <- sdf %>% gamm(abs(IS) ~ s(cause_salinity), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, cause_var = ~1), niterPQL = 20)
gamm11 <- sdf %>% gamm(abs(IS) ~ s(cause_tide), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, cause_var = ~1), niterPQL = 20)
gamm12 <- sdf %>% gamm(abs(IS) ~ s(cause_wave), data = ., family = Gamma(link = "log"),
                       random=list(site_code = ~1, cause_var = ~1), niterPQL = 20)

## Check results
### Effect IS
summary(gamm01$lme); summary(gamm01$gam) # water_temp; LME p = 0.2419(+), GAM p < 2e-16
summary(gamm02$lme); summary(gamm02$gam) # sp_richness; LME p = 0.4549(+), GAM p = 0.0631
summary(gamm03$lme); summary(gamm03$gam) # total_dna; LME p = 0.153(-), GAM p < 2e-16
### Causal IS
summary(gamm04$lme); summary(gamm04$gam) # water_temp; LME p = 0.4793(+), GAM p < 2e-16
summary(gamm05$lme); summary(gamm05$gam) # sp_richness; LME p = 0(-), GAM p = 8.24e-07
summary(gamm06$lme); summary(gamm06$gam) # total_dna; LME p = 0.0042(-), GAM p < 2e-16
### Effect IS
summary(gamm07$lme); summary(gamm07$gam) # salinity; LME p = 0.6566(+), GAM p = 2.18e-05
summary(gamm08$lme); summary(gamm08$gam) # tide; LME p = 0(-), GAM p = 8.36e-07
summary(gamm09$lme); summary(gamm09$gam) # wave_m; LME p = 0.1713(+), GAM p = 0.00234
### Caussal IS
summary(gamm10$lme); summary(gamm10$gam) # salinity; LME p = 0.3703(1), GAM p = 0.000619
summary(gamm11$lme); summary(gamm11$gam) # tide; LME p = 0.5173(+), GAM p = 0.517
summary(gamm12$lme); summary(gamm12$gam) # wave_m; LME p = 0.0049(+), GAM p = 6.53e-06

# Save summaries
gamm01_sum <- list(summary(gamm01$lme), summary(gamm01$gam))
gamm02_sum <- list(summary(gamm02$lme), summary(gamm02$gam))
gamm03_sum <- list(summary(gamm03$lme), summary(gamm03$gam))
gamm04_sum <- list(summary(gamm04$lme), summary(gamm04$gam))
gamm05_sum <- list(summary(gamm05$lme), summary(gamm05$gam))
gamm06_sum <- list(summary(gamm06$lme), summary(gamm06$gam))

gamm07_sum <- list(summary(gamm07$lme), summary(gamm07$gam))
gamm08_sum <- list(summary(gamm08$lme), summary(gamm08$gam))
gamm09_sum <- list(summary(gamm09$lme), summary(gamm09$gam))
gamm10_sum <- list(summary(gamm10$lme), summary(gamm10$gam))
gamm11_sum <- list(summary(gamm11$lme), summary(gamm11$gam))
gamm12_sum <- list(summary(gamm12$lme), summary(gamm12$gam))



# Get predicted effects of explaining variables
# In-strength (main figure)
pred01 <- predict(gamm01$gam, type = "response", se = T, newdata = data.frame(water_temp = seq(min(sdf$water_temp, na.rm = T), max(sdf$water_temp, na.rm = T), length.out = 200)))
pred02 <- predict(gamm02$gam, type = "response", se = T, newdata = data.frame(sp_richness = seq(min(sdf$sp_richness, na.rm = T), max(sdf$sp_richness, na.rm = T), length.out = 200)))
pred03 <- predict(gamm03$gam, type = "response", se = T, newdata = data.frame(total_dna_per_ml = seq(min(sdf$total_dna_per_ml, na.rm = T), max(sdf$total_dna_per_ml, na.rm = T), length.out = 200)))
# Out-strength (main figure)
pred04 <- predict(gamm04$gam, type = "response", se = T, newdata = data.frame(cause_temp = seq(min(sdf$cause_temp, na.rm = T), max(sdf$cause_temp, na.rm = T), length.out = 200)))
pred05 <- predict(gamm05$gam, type = "response", se = T, newdata = data.frame(cause_richness = seq(min(sdf$cause_richness, na.rm = T), max(sdf$cause_richness, na.rm = T), length.out = 200)))
pred06 <- predict(gamm06$gam, type = "response", se = T, newdata = data.frame(cause_totdna = seq(min(sdf$cause_totdna, na.rm = T), max(sdf$cause_totdna, na.rm = T), length.out = 200)))
# In-strength (SI figure)
pred07 <- predict(gamm07$gam, type = "response", se = T, newdata = data.frame(salinity = seq(min(sdf$salinity, na.rm = T), max(sdf$salinity, na.rm = T), length.out = 200)))
pred08 <- predict(gamm08$gam, type = "response", se = T, newdata = data.frame(tide = seq(min(sdf$tide, na.rm = T), max(sdf$tide, na.rm = T), length.out = 200)))
pred09 <- predict(gamm09$gam, type = "response", se = T, newdata = data.frame(wave_m = seq(min(sdf$wave_m, na.rm = T), max(sdf$wave_m, na.rm = T), length.out = 200)))
# Out-strength (SI figure)
pred10 <- predict(gamm10$gam, type = "response", se = T, newdata = data.frame(cause_salinity = seq(min(sdf$cause_salinity, na.rm = T), max(sdf$cause_salinity, na.rm = T), length.out = 200)))
pred11 <- predict(gamm11$gam, type = "response", se = T, newdata = data.frame(cause_tide = seq(min(sdf$cause_tide, na.rm = T), max(sdf$cause_tide, na.rm = T), length.out = 200)))
pred12 <- predict(gamm12$gam, type = "response", se = T, newdata = data.frame(cause_wave = seq(min(sdf$cause_wave, na.rm = T), max(sdf$cause_wave, na.rm = T), length.out = 200)))


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Delete large files
rm(sdf_all); rm(intmat_all); rm(intmat_only_sp); rm(intmat_only_sp_abs)
rm(gamm01); rm(gamm02); rm(gamm03); rm(gamm04); rm(gamm05); rm(gamm06)
rm(gamm07); rm(gamm08); rm(gamm09); rm(gamm10); rm(gamm11); rm(gamm12)

# Save workspace and objects
saveRDS(list(gamm01_sum, gamm02_sum, gamm03_sum, gamm04_sum, gamm05_sum, gamm06_sum,
             gamm07_sum, gamm08_sum, gamm09_sum, gamm10_sum, gamm11_sum, gamm12_sum),
        sprintf("%s/gamm_summaries.obj", output_folder))
saveRDS(list(pred01, pred02, pred03, pred04, pred05, pred06,
             pred07, pred08, pred09, pred10, pred11, pred12),
        sprintf("%s/pred_all.obj", output_folder))
save(list = ls(all.names = TRUE),
      file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
macam::save_session_info(session_info_dir = "00_SessionInfo")
