####
#### Boso Peninsula project
#### No. 7 Summarize and Visualize UIC results
#### 2022.11.9 Ushio (R4.2.1)
####

# Load workspace
load("05_MDRSmap_SpOut/05_MDRSmap_SpOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.8
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
library(ggraph); packageVersion("ggraph") # 2.1.0, 2022.11.8
library(igraph); packageVersion("igraph") # 1.3.5, 2022.11.8
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

# Custom package
library(macam); packageVersion("macam") # 0.0.10, 2022.11.9


# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
(output_folder <- macam::outdir_create())


# <---------------------------------------------> #
# Check MDR S-map results
# <---------------------------------------------> #
asv_df_occur <- asv_df_conv_sort %>% select(top_taxa) %>% 
  mutate(site_code = sample_df_sort$site_code) %>%
  pivot_longer(cols = -site_code, names_to = "species") %>%
  group_by(species, site_code) %>% 
  summarize(mean = mean(value), n = sum(value > 0)) %>% 
  data.frame


# <---------------------------------------------> #
# Prepare library indices (TRUE/FALSE) and time indices
# <---------------------------------------------> #
# Add NA to separate study sites
valid_lib_idx <- uic_lib_NAadd
valid_idx <- c()
for (i in 1:nrow(valid_lib_idx)) valid_idx <- c(valid_idx, valid_lib_idx[i,1]:valid_lib_idx[i,2])
valid_idx_df <- cbind(data.frame(valid_idx = valid_idx), sample_df_sort)


# <---------------------------------------------> #
# Preparation for reconstructing interaction matrix
# <---------------------------------------------> #
## Structure of the interaction matrix
##                 <--- Causal species --->
##        |        x  x  x  x  x ........ x 
##        |        x  x  x  x  x ........ x 
## Effect species  x  x  x  x  x ........ x 
##        |        x  x  x  x  x ........ x 
##        |        x  x  x  x  x ........ x 
##
# Prepare an empty interaction matrix
env_vars <- c("water_temp","salinity","wave_m","tide")
intmat_all <- list()
intmat_i <- matrix(0, ncol = length(top_taxa) + length(env_vars), nrow = length(top_taxa))
rownames(intmat_i) <- top_taxa
colnames(intmat_i) <- c(top_taxa, env_vars)
# Append empty matrices
for (i in 1:length(valid_idx)) {
  intmat_all <- c(intmat_all, list(intmat_i))
}
# Assign name
names(intmat_all) <- paste0(sample_df_sort$site_code, "_", sample_df_sort$date)


# <---------------------------------------------> #
# Summarize S-map coefficients to data.frame
# <---------------------------------------------> #
# Assign interaction strength to each element
for (i in 1:nrow(opt_mdr_res)) {
  ## Extract an effect variable
  tax_i <- opt_mdr_res$effect_var[i]
  ## Extract embedding information
  block_vars_i <- opt_mdr_res$block_mvd[i][[1]] %>%
    colnames %>% str_split(pattern = "_tp") %>% sapply(`[`, 1)
  # Remove time-delayed self-interactions
  tax_valid_idx <- c(1, which(block_vars_i != tax_i))
  # Extract tp information
  tp_vars_i <- opt_mdr_res$block_mvd[i][[1]] %>%
    colnames %>% str_split(pattern = "_tp") %>% sapply(`[`, 2) %>% 
    as.numeric
  ## Extract S-map coefficients column names
  smapc_colnames <- paste0(sprintf("c_%s", 1:length(block_vars_i)))
  ## Extract S-map coefficients
  smapc_df <- opt_mdr_res$smap_coefficients[i][[1]] %>%
    select(all_of(smapc_colnames)) %>% .[valid_idx_df$valid_idx,]

  ## Make delayed block to check whether a focal species pair is
  block_delay_edna <- block_vars_i[tax_valid_idx][!(block_vars_i[tax_valid_idx] %in% env_vars)] %>% 
    asv_df_NAadd[.]
  block_delay_env <- block_vars_i[tax_valid_idx][(block_vars_i[tax_valid_idx] %in% env_vars)] %>% 
    sample_df_NAadd[.]
  block_delay <- cbind(block_delay_edna, block_delay_env)
  for (block_i in 1:ncol(block_delay)) {
    tp_i <- tp_vars_i[tax_valid_idx][block_i]
    block_delay[block_i] <- dplyr::lag(block_delay[block_i], n = abs(tp_i))
  }
  ## Make delayed environmental variables (temp, richness, dna, wave, salinity, tide)
  block_delay_temp <- block_delay_wave <- block_delay_salinity <-
    block_delay_tide <- block_delay_richness <- block_delay_totdna <- block_delay
  for (block_i in 1:ncol(block_delay)) {
    tp_i <- tp_vars_i[tax_valid_idx][block_i]
    block_delay_temp[block_i] <- dplyr::lag(sample_df_NAadd$water_temp, n = abs(tp_i))
    block_delay_wave[block_i] <- dplyr::lag(sample_df_NAadd$wave_m, n = abs(tp_i))
    block_delay_salinity[block_i] <- dplyr::lag(sample_df_NAadd$salinity, n = abs(tp_i))
    block_delay_tide[block_i] <- dplyr::lag(sample_df_NAadd$tide, n = abs(tp_i))
    block_delay_richness[block_i] <- dplyr::lag(apply(asv_df_NAadd, 1, function(x) sum(x>0)), n = abs(tp_i))
    block_delay_totdna[block_i] <- dplyr::lag(sample_df_NAadd$total_dna_per_ml, n = abs(tp_i))
  }
  # Trim block_delay_xxxx
  block_delay <- block_delay[valid_idx,]
  block_delay_temp <- block_delay_temp[valid_idx,]
  block_delay_wave <- block_delay_wave[valid_idx,]
  block_delay_salinity <- block_delay_salinity[valid_idx,]
  block_delay_tide <- block_delay_tide[valid_idx,]
  block_delay_richness <- block_delay_richness[valid_idx,]
  block_delay_totdna <- block_delay_totdna[valid_idx,]
    
  ## Create one tidy data.frame that includes all information
  smapc_df_tmp <- smapc_df[,tax_valid_idx]
  colnames(smapc_df_tmp) <- sprintf("effect_from_%s", block_vars_i[tax_valid_idx])
  #colnames(smapc_df_tmp) <- block_vars_i[tax_valid_idx]
  smapc_df_tmp <- cbind(smapc_df_tmp, sample_df_sort)
  smapc_df_tmp$sample_id <- rownames(sample_df_sort)
  smapc_df_tmp$effect_var <- tax_i
  smapc_df_tmp$effect_var_edna <- block_delay[,tax_i]
  smapc_df_tmp$sp_richness <- apply(asv_df_conv_sort, 1, function(x) sum(x > 0))
  smapc_df_tmp_long <- pivot_longer(smapc_df_tmp,
                                    cols = -c(c(colnames(sample_df_sort), "sample_id", "effect_var", "effect_var_edna",
                                                "sp_richness")),
                                    names_to = "cause_var", values_to = "IS")
  # Assign causal var's edna conc.
  smapc_df_tmp_long$cause_var_edna <- NA
  smapc_df_tmp_long$cause_temp <- NA
  smapc_df_tmp_long$cause_wave <- NA
  smapc_df_tmp_long$cause_salinity <- NA
  smapc_df_tmp_long$cause_tide <- NA
  smapc_df_tmp_long$cause_richness <- NA
  smapc_df_tmp_long$cause_totdna <- NA
  for (j in 1:nrow(smapc_df_tmp_long)) {
    cause_var_i <- smapc_df_tmp_long$cause_var[j] %>% 
      str_split(pattern = "effect_from_") %>% .[[1]] %>% .[2]
    sample_id_i <- smapc_df_tmp_long$sample_id[j]
    # Assign values to the data.frame
    smapc_df_tmp_long$cause_var_edna[j] <- block_delay[sample_id_i, cause_var_i]
    smapc_df_tmp_long$cause_temp[j] <- block_delay_temp[sample_id_i, cause_var_i]
    smapc_df_tmp_long$cause_wave[j] <- block_delay_wave[sample_id_i, cause_var_i]
    smapc_df_tmp_long$cause_salinity[j] <- block_delay_salinity[sample_id_i, cause_var_i]
    smapc_df_tmp_long$cause_tide[j] <- block_delay_tide[sample_id_i, cause_var_i]
    smapc_df_tmp_long$cause_richness[j] <- block_delay_richness[sample_id_i, cause_var_i]
    smapc_df_tmp_long$cause_totdna[j] <- block_delay_totdna[sample_id_i, cause_var_i]
  }
  valid_smapc_cond1 <- smapc_df_tmp_long$cause_var_edna > 0 & !is.na(smapc_df_tmp_long$cause_var_edna)
  valid_smapc_cond2 <- smapc_df_tmp_long$effect_var_edna > 0 & !is.na(smapc_df_tmp_long$effect_var_edna)
  valid_smapc_cond3 <- !is.na(smapc_df_tmp_long$IS)
  
  if (i == 1) {
    smapc_df_long <- smapc_df_tmp_long[valid_smapc_cond1 & valid_smapc_cond2 & valid_smapc_cond3,] %>%
      as.data.frame
  } else {
    smapc_df_tmp_long <- smapc_df_tmp_long[valid_smapc_cond1 & valid_smapc_cond2 & valid_smapc_cond3,] %>% 
      as.data.frame
    row.names(smapc_df_tmp_long) <- as.character((nrow(smapc_df_long)+1):(nrow(smapc_df_long)+nrow(smapc_df_tmp_long)))
    smapc_df_long <- rbind(smapc_df_long, smapc_df_tmp_long)
  }
}
# Delete temporal objects
rm(smapc_df_tmp); rm(tax_i); rm(cause_var_i); rm(sample_id_i)
rm(smapc_df_tmp_long); rm(smapc_colnames); rm(tp_i); rm(tp_vars_i)
rm(valid_smapc_cond1); rm(valid_smapc_cond2); rm(valid_smapc_cond3)
rm(i); rm(j)
rm(block_delay); rm(block_delay_edna); rm(block_delay_env)
rm(block_delay_richness); rm(block_delay_salinity)
rm(block_delay_temp); rm(block_delay_tide); rm(block_delay_totdna)

# Rename smapc_df_long to sdf
sdf_all <- smapc_df_long; rm(smapc_df_long)
# Rename cause_var
sdf_all$cause_var <- sapply(sdf_all$cause_var, function(x) {str_split(x, pattern = "effect_from_") %>% .[[1]] %>% .[2]})
# Remove self-regulation effects and environment effects
sdf <- sdf_all[sdf_all$effect_var != sdf_all$cause_var,]
sdf <- sdf[!(sdf$cause_var %in% env_df_var),]
# Check dim
dim(sdf_all); dim(sdf)


# <---------------------------------------------> #
# Assign values to interaction matrix
# <---------------------------------------------> #
site_date_code <- paste0(sdf$site_code, "_", sdf$date)
for (i in 1:nrow(sdf)) {
  site_date_i <- site_date_code[i]
  intmat_all[[site_date_i]][sdf_all$effect_var[i], sdf_all$cause_var[i]] <- sdf_all$IS[i]
}

# Assign NA to diagonal elements and trim environmental effects
intmat_only_sp <- list()
for (i in 1:length(intmat_all)) {
  intmat_only_sp <- c(intmat_only_sp, list(intmat_all[[i]][top_taxa, top_taxa]))
  diag(intmat_only_sp[[i]]) <- NA
}
# Assign NA to interaction strength = 0
for (i in 1:length(intmat_only_sp)) intmat_only_sp[[i]][intmat_only_sp[[i]] == 0] <- NA
intmat_only_sp_abs <- map(intmat_only_sp, abs)

# Add the total number of interaction detected
n_int_df <- sdf %>% group_by(sample_id) %>% summarize(n_int = n()) %>% data.frame
sdf_all$n_int <- sdf$n_int <- 0
for (i in 1:nrow(sdf_all)) sdf_all$n_int[i] <- n_int_df[n_int_df$sample_id == sdf_all$sample_id[i],"n_int"]
for (i in 1:nrow(sdf)) sdf$n_int[i] <- n_int_df[n_int_df$sample_id == sdf$sample_id[i],"n_int"]


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save results
# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))

