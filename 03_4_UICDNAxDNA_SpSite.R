####
#### Boso Peninsula project
#### No.3-4: Detecting causal pairs using UIC (DNA xmap DNA, species x site)
####

# Load workspace
load("01_CompileDataOut/01_CompileDataOut.RData")

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.10.3
library(rUIC); packageVersion("rUIC") # 0.1.5, 2020.10.30
library(pforeach); packageVersion("pforeach") # 1.3, 2021.2.10
library(cowplot); packageVersion("cowplot") # 1.1.0, 2021.1.7
theme_set(theme_cowplot())

# Set random seeds and create output folder
ran_seed <- 8181
set.seed(ran_seed)
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
output_folder03 <- "03_UICEnvDNAOut"
dir.create(output_folder03)


# <---------------------------------------------> #
#                 Preprocessing                   #
# <---------------------------------------------> #
## Extract dominant top fish spp based on detection frequency
taxa_orders_conc <- names(sort(colSums(asv_df_conv), decreasing = T))
top_taxa <- names(sort(colSums(asv_df_conv > 0), decreasing = T)) %>% .[1:50]
tax_df[top_taxa,]; asv_df_conv[,top_taxa]

## Add sampling time index for time series analysis
sample_df$time_index <- as.numeric(str_sub(sample_df$time_code, start = 3, end = 4))

# Sort data by site
sample_df_sort <- sample_df %>% dplyr::arrange(site_code)
asv_df_conv_sort <- asv_df_conv %>%
  dplyr::arrange(sample_df$site_code) %>%
  select(top_taxa); dim(asv_df_conv_sort)
all(rownames(sample_df_sort) == rownames(asv_df_conv_sort))

# Sort ASV data by site x species
asv_df_sitesp <- tibble(time_index = unique(sample_df_sort$time_index))
for(fish_sp_i in top_taxa){
  # Collect eDNA data
  asv_df_sitesp_tmp <- asv_df_conv_sort[,fish_sp_i] %>%
    as_tibble() %>%
    bind_cols(site_code = sample_df_sort$site_code) %>%
    bind_cols(time_index = sample_df_sort$time_index) %>%
    pivot_wider(id_cols = time_index, names_from = site_code) %>%
    .[,-1]
  ## Rename colnames
  colnames(asv_df_sitesp_tmp) <- colnames(asv_df_sitesp_tmp) %>%
    sprintf("%s_%s", fish_sp_i, .)
  asv_df_sitesp <- asv_df_sitesp %>% bind_cols(asv_df_sitesp_tmp)
}
# Remove time_index
asv_df_sitesp <- asv_df_sitesp %>% select(!time_index)

# Sort sample_data by site x species
env_var <- c("water_temp", "salinity", "wave_m", "tide")
sample_df_sitesp <- tibble(time_index = unique(sample_df_sort$time_index))
for(env_i in env_var){
  # Collect eDNA data
  sample_df_sitesp_tmp <- sample_df_sort[,env_i] %>%
    as_tibble() %>%
    bind_cols(site_code = sample_df_sort$site_code) %>%
    bind_cols(time_index = sample_df_sort$time_index) %>%
    pivot_wider(id_cols = time_index, names_from = site_code) %>%
    .[,-1]
  ## Rename colnames
  colnames(sample_df_sitesp_tmp) <- colnames(sample_df_sitesp_tmp) %>%
    sprintf("%s_%s", env_i, .)
  sample_df_sitesp <- sample_df_sitesp %>% bind_cols(sample_df_sitesp_tmp)
}
# Remove time_index
sample_df_sitesp <- sample_df_sitesp %>% select(!time_index)


# <---------------------------------------------> #
#          Load causally-related env data        #
# <---------------------------------------------> #
# Load eDNA uic Climate data
uic_edna_env <- readRDS(sprintf("%s/uic_edna_env.obj", output_folder03))
uic_edna_env_signif <- uic_edna_env %>% filter(pval <= 0.05) %>% filter(tp <= 0)


# <---------------------------------------------> #
#               Set basic parameters              #
# <---------------------------------------------> #
# Basic parameters
TP_RANGE <- seq(2, -6, by = -1)
E_RANGE <- 0:10


# <---------------------------------------------> #
#                     Main loop                   #
# <---------------------------------------------> #
# Prepare output object
uic_edna_edna_site0 <- data.frame()
total_process <- ncol(asv_df_sitesp) #* length(TP_RANGE)
process_i <- 1

for (edna_i in 1:ncol(asv_df_sitesp)){
  # Record start time
  start_time <- proc.time()[3]
  
  # Set potential causal variable
  x_std <- data.frame(x = as.numeric(scale(asv_df_sitesp[,edna_i])))
  
  # Prepare conditional variables
  taxa_i <- str_sub(colnames(asv_df_sitesp)[edna_i], end = -6)
  site_i <- str_sub(colnames(asv_df_sitesp)[edna_i], start = 10)
  env_cause_var <- unique((uic_edna_env_signif %>% filter(effect_var == taxa_i))$cause_var)
  env_tp_df <- NULL
  if(length(env_cause_var) > 0){
    env_tp_df <- data.frame(NULL)
    for(env_cause_i in env_cause_var){
      env_tp <- uic_edna_env_signif %>% filter(effect_var == taxa_i) %>%
        filter(cause_var == env_cause_i) %>% .[which.max(.[,"te"]), "tp"] %>% as.integer()
      env_tp_df <- rbind(env_tp_df, data.frame(env_var = env_cause_i, tp = env_tp))
      rm(env_tp)
    }
  }
  # Generate conditional variables data.frame
  z_std <- NULL
  if(!is.null(env_tp_df)){
    z_std <- tibble(nan_col = rep(NaN, nrow(sample_df_sitesp)))
    for(z_i in 1:nrow(env_tp_df)){
      env_cause_i <- paste0(env_tp_df[z_i, "env_var"], "_", site_i)
      env_tp_i <- env_tp_df[z_i, "tp"]
      if(env_tp_i > -1){
        z_std_tmp <- data.frame(z_tmp = dplyr::lead(as.numeric(scale(sample_df_sitesp[,env_cause_i])), n = env_tp_i + 1))
      }else{
        z_std_tmp <- data.frame(z_tmp = dplyr::lag(as.numeric(scale(sample_df_sitesp[,env_cause_i])), n = - env_tp_i - 1))
      }
      z_std <- z_std %>% bind_cols(z_std_tmp)
      colnames(z_std)[ncol(z_std)] <- sprintf("z%01d", z_i)
    }
    z_std <- z_std[,-1]
  }

  # UIC main loop
  # Parallel computing for eDNA taxa
  uic_edna_edna_site_tmp <- pforeach(edna_j = 1:ncol(asv_df_sitesp), .c=rbind)({
    # Set potential causal variable
    y_std <- data.frame(y = as.numeric(scale(asv_df_sitesp[,edna_j])))
    
    # Determine an optimal E
    if(!is.null(env_tp_df)){
      simp_res <- rUIC::simplex(cbind(x_std, y_std, z_std), lib_var = "x",
                                cond_var = c("y", colnames(z_std)),
                                lib = c(1,nrow(cbind(x_std, y_std, z_std))),
                                E = E_RANGE, tp = 1, tau = 1, 
                                Enull = "adaptive", n_boot = 2000, seed = 12345)
      Eopt <- with(simp_res, max(c(0, E[pval < 0.05])))
      # Calculate UIC
      uic_res_tmp <- rUIC::uic(cbind(x_std, y_std, z_std), lib_var = "x",
                               tar_var = "y", cond_var = colnames(z_std),
                               E = Eopt + 1, tp = TP_RANGE,
                               lib = c(1,nrow(cbind(x_std, y_std, z_std))),
                               n_boot = 2000, seed = 1234) %>%
        cbind(data.frame(effect_var = rep(colnames(asv_df_sitesp)[edna_i], length(TP_RANGE)),
                         cause_var =  rep(colnames(asv_df_sitesp)[edna_j], length(TP_RANGE))), .,
              data.frame(cond1 = NA, cond2 = NA, cond3 = NA, cond4 = NA))
      # Add conditional variables
      for(env_cause_i in 1:nrow(env_tp_df)){
        uic_res_tmp[,c("cond1","cond2","cond3","cond4")[env_cause_i]] <- 
          sprintf("%s; tp = %d", env_tp_df[env_cause_i, "env_var"], env_tp_df[env_cause_i, "tp"])
      }
    
    } else {
      simp_res <- rUIC::simplex(cbind(x_std, y_std), lib_var = "x",
                                cond_var = c("y"),
                                lib = c(1, nrow(cbind(x_std, y_std))),
                                E = E_RANGE, tp = 1, tau = 1, 
                                Enull = "adaptive", n_boot = 2000, seed = 12345)
      Eopt <- with(simp_res, max(c(0, E[pval < 0.05])))
      # Calculate UIC
      uic_res_tmp <- rUIC::uic(cbind(x_std, y_std), lib_var = "x",
                               tar_var = "y",
                               E = Eopt + 1, tp = TP_RANGE,
                               lib = c(1, nrow(cbind(x_std, y_std))),
                               n_boot = 2000, seed = 1234) %>%
        cbind(data.frame(effect_var = rep(colnames(asv_df_sitesp)[edna_i], length(TP_RANGE)),
                         cause_var =  rep(colnames(asv_df_sitesp)[edna_j], length(TP_RANGE))), .,
              data.frame(cond1 = NA, cond2 = NA, cond3 = NA, cond4 = NA))
    }
    
    uic_res_tmp
  })
  
  # Combine results
  uic_edna_edna_site0 <- rbind(uic_edna_edna_site0, uic_edna_edna_site_tmp)
  
  # Delete temporal objects
  #rm(start_time); rm(elapsed_time)
  rm(uic_edna_edna_site_tmp)
  rm(x_std) #; rm(y_std)
  rm(z_std)
  #rm(Eopt); rm(simp_res)
  
  # Output process message
  elapsed_time <- round(proc.time()[3] - start_time, 2)
  message(paste("Process", process_i, "/", total_process, "finished:", elapsed_time, "sec elapsed"))
  process_i <- process_i + 1
}

# Delete temporal objects
rm(start_time); rm(elapsed_time)
rm(process_i); rm(total_process)
#rm(uic_edna_edna_tmp)


# <---------------------------------------------> #
#                   Save results                  #
# <---------------------------------------------> #
# Save and output results
saveRDS(uic_edna_edna_site0, sprintf("%s/uic_edna_edna_site.obj", output_folder03))
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder03, output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_%s.txt",  output_folder, substr(Sys.time(), 1, 10)))

