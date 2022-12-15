####
#### Boso Peninsula project
#### No.3-1: Detecting causal pairs using UIC (DNA xmap Env)
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
# Check data structure
head(asv_df_conv)
head(sample_df)
all(rownames(asv_df_conv) == rownames(sample_df))
dim(asv_df); dim(tax_df); dim(sample_df)

## Extract dominant top fish spp based on detection frequency
taxa_orders_conc <- names(sort(colSums(asv_df_conv), decreasing = T))
top_taxa <- names(sort(colSums(asv_df_conv > 0), decreasing = T)) %>% .[1:50]
tax_df[top_taxa,]; asv_df_conv[,top_taxa]

## Add sampling time index for time series analysis
sample_df$time_index <- as.numeric(str_sub(sample_df$time_code, start = 3, end = 4))

# Sort data.frame
sample_df_sort <- sample_df %>% dplyr::arrange(site_code)
asv_df_conv_sort <- asv_df_conv %>% dplyr::arrange(sample_df$site_code)
all(rownames(sample_df_sort) == rownames(asv_df_conv_sort))
# Identify library index for asv_df_conv_sort
start_time_id <- which(sample_df_sort$time_code == "BP01")
end_time_id <- c(start_time_id[2:length(start_time_id)] - 1, nrow(sample_df))
uic_lib <- matrix(c(start_time_id, end_time_id), ncol = 2, byrow = F)

# Basic paramters
TP_RANGE <- seq(2, -6, by = -1)
E_RANGE <- 0:10
env_var <- c("water_temp", "salinity", "wave_m", "tide")

# <---------------------------------------------> #
#                     Main loop                   #
# <---------------------------------------------> #
# Prepare output object
uic_edna_env0 <- data.frame()
total_process <- length(env_var) #* length(TP_RANGE)
process_i <- 1

for (env_i in env_var){
  # Record start time
  start_time <- proc.time()[3]
  
  # Set potential causal variable
  y_std <- data.frame(y = as.numeric(scale(sample_df_sort[env_i])))
  
  # Parallel computing for eDNA taxa
  uic_edna_env_tmp <- pforeach(col_j = 1:length(top_taxa), .c=rbind)({
    ## Prepare output object
    uic_edna_env00 <- data.frame()
    fish_ts_i <- asv_df_conv_sort %>% select(top_taxa[col_j])
    ## Prepare causal variable (eDNA data)
    x_std <- data.frame(x = as.numeric(scale(fish_ts_i)))
    ## Determine an optimal E
    simp_res <- rUIC::simplex(cbind(x_std, y_std), lib_var = "x", cond_var = "y",
                              lib = uic_lib, E = E_RANGE, tp = 1, tau = 1, 
                              Enull = "adaptive", n_boot = 2000, seed = 1234)
    Eopt <- with(simp_res, max(c(0, E[pval < 0.05])))
    ## Calculate UIC
    for (tp_i in TP_RANGE){
      uic_res_tmp <- rUIC::uic(cbind(x_std, y_std), lib_var = "x", tar_var = "y",
                               E = Eopt + 1, tp = tp_i, lib = uic_lib,
                               n_boot = 2000, seed = 1234) %>%
        cbind(data.frame(effect_var = top_taxa[col_j],
                         cause_var = env_i), .)
      uic_edna_env00 <- rbind(uic_edna_env00, uic_res_tmp)
    }
    
    uic_edna_env00
  })
  
  # Combine results
  uic_edna_env0 <- rbind(uic_edna_env0, uic_edna_env_tmp)
  
  # Output process message
  elapsed_time <- round(proc.time()[3] - start_time, 2)
  message(paste("Process", process_i, "/", total_process, "finished:", elapsed_time, "sec elapsed"))
  process_i <- process_i + 1
}

# Delete temporal objects
rm(start_time); rm(elapsed_time)
rm(y_std)#; rm(x_std); rm(uic_res_tmp)
rm(process_i); rm(total_process)
rm(uic_edna_env_tmp)


# <---------------------------------------------> #
#                   Save results                  #
# <---------------------------------------------> #
# Save and output results
saveRDS(top_taxa, sprintf("%s/top_taxa.obj", output_folder03))
saveRDS(uic_edna_env0, sprintf("%s/uic_edna_env.obj", output_folder03))
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder03, output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_%s.txt",  output_folder, substr(Sys.time(), 1, 10)))

