####
#### Boso Peninsula project
#### No. 5 MDR S-map 
#### 2022.11.04 Ushio (R4.2.1)
#### 2022.11.08 
####

# Load workspace
load("04_CompileUICresOut/04_CompileUICresOut_Minimal.RData") 

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.10.23
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
theme_set(theme_cowplot())

# Load custom library
library(macam); packageVersion("macam") # 0.0.10, 2022.11.08

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
(output_folder <- macam::outdir_create())


# ----------------------------------------------------- #
# Load UIC summary objects
# ----------------------------------------------------- #
edna_env_summary <- read_csv("04_CompileUICresOut/uic_edna_env_strong.csv")
env_edna_summary <- read_csv("04_CompileUICresOut/uic_env_edna_strong.csv")
edna_sp_summary <- read_csv("04_CompileUICresOut/uic_edna_sp_strong.csv")


# ----------------------------------------------------- #
# Do MDR S-map: Species-level analysis (sites aggregated)
# ----------------------------------------------------- #
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
# Add NA to separate study sites
NA_nrow <- 10
NA_block_asv <- data.frame(matrix(rep(NA, ncol(asv_df_conv_sort)*NA_nrow), nrow = NA_nrow))
NA_block_sample <- data.frame(matrix(rep(NA, ncol(sample_df_sort)*NA_nrow), nrow = NA_nrow))
colnames(NA_block_asv) <- colnames(asv_df_conv_sort)
colnames(NA_block_sample) <- colnames(sample_df_sort)
uic_lib_NAadd <- uic_lib # To update idx
for (i in 1:nrow(uic_lib)) {
  if (i == 1) {
    asv_df_tmp <- asv_df_conv_sort[uic_lib[i,1]:uic_lib[i,2],]
    sample_df_tmp <- sample_df_sort[uic_lib[i,1]:uic_lib[i,2],]
    asv_df_NAadd <- rbind(asv_df_tmp, NA_block_asv)
    sample_df_NAadd <- rbind(sample_df_tmp, NA_block_sample)
  } else {
    uic_lib_NAadd[i,1] <- nrow(asv_df_NAadd) + 1
    asv_df_tmp <- asv_df_conv_sort[uic_lib[i,1]:uic_lib[i,2],]
    sample_df_tmp <- sample_df_sort[uic_lib[i,1]:uic_lib[i,2],]
    asv_df_NAadd <- rbind(asv_df_NAadd, asv_df_tmp, NA_block_asv)
    sample_df_NAadd <- rbind(sample_df_NAadd, sample_df_tmp, NA_block_sample)
    uic_lib_NAadd[i,2] <- nrow(asv_df_NAadd) - NA_nrow
  }
}


# <---------------------------------------------> #
# Main loop
# <---------------------------------------------> #
# Normalize data
asv_df_NAadd_std <- as.data.frame(apply(asv_df_NAadd, 2, function(x) as.numeric(scale(x))))

# Tentative parameters
theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
alpha_range <- 0 #c(0,1)
regularized_i <- TRUE
step_i <- 1
total_step <- length(top_taxa)*length(theta_range)*length(lambda_range)*length(alpha_range)
step_for_taxa <- length(theta_range)*length(lambda_range)*length(alpha_range)

for (tax_i in top_taxa) {
  # Estimate optimal embedding dimension for MVD
  simp_x <- rUIC::simplex(asv_df_NAadd_std, lib_var = tax_i, E = 2:10, tp = 1)
  (Ex <- simp_x[which.min(simp_x$rmse),"E"])
  
  for (theta_i in theta_range) {
    for (lambda_i in lambda_range) {
      for (alpha_i in alpha_range) {
        # Record starting time
        start_time <- proc.time()[3]
        
        # Extract causal environmental variables
        env_df_tmp <- as.data.frame(edna_sp_summary[edna_sp_summary$effect_var == tax_i,14:17][1,])
        env_df_var <- sapply(str_split(env_df_tmp, pattern = ";"), `[`, 1)
        env_df_var <- env_df_var[!is.na(env_df_var)]
        env_df_tp <- sapply(str_split(env_df_tmp, pattern = "; tp = "), `[`, 2)
        env_df_tp <- as.numeric(env_df_tp[!is.na(env_df_tp)])
        if (length(env_df_var) > 0) {
          for (env_i in 1:length(env_df_var)) {
            if (env_i == 1) {
              env_block_mvd <- as.data.frame(dplyr::lag(sample_df_NAadd[,env_df_var[env_i]], n = abs(env_df_tp[env_i])))
            } else {
              env_block_tmp <- as.data.frame(dplyr::lag(sample_df_NAadd[,env_df_var[env_i]], n = abs(env_df_tp[env_i])))
              env_block_mvd <- cbind(env_block_mvd, env_block_tmp)
            }
          }
          colnames(env_block_mvd) <- sprintf("%s_tp%s", env_df_var, env_df_tp)
          # Normalize environemntal variables
          env_block_mvd <- as.data.frame(apply(env_block_mvd, 2, function(x) as.numeric(scale(x)))) 
        }
        
        # Make MVD block (only strongest causality is included)
        block_mvd <- macam::make_block_mvd(asv_df_NAadd_std,
                                           as.data.frame(edna_sp_summary[edna_sp_summary$effect_var == tax_i,]),
                                           effect_var = tax_i,
                                           E_effect_var = Ex,
                                           cause_var_colname = "cause_var",
                                           include_var = "strongest_only",
                                           p_threshold = 0.05,
                                           sort_tp = T, silent = TRUE)
        if (length(env_df_var) > 0) {
          # Merge environmental variables and ASV block
          block_mvd <- cbind(block_mvd, env_block_mvd)
        }
        
        # Compute multiview distance
        dist_mvd_list <- try(macam::compute_mvd(block_mvd, effect_var = tax_i, Ex, tp = 1,
                                                make_block_method = "naive",
                                                n_ssr = 1000, k = floor(sqrt(1000)),
                                                random_seed = 1234, distance_only = FALSE,
                                                silent = TRUE))
        if (class(dist_mvd_list) == "try-error") {
          dist_mvd <- dist_mvd_list$parms <- dist_mvd_list$top_embeddings <- NA
          mdr_res <- list(stats = data.frame(N=NA, rho=NA, mae=NA, rmse=NA, pval=NA))
        } else {
          dist_mvd <- dist_mvd_list$multiview_dist
          
          # Perform MDR S-map
          mdr_res <- try(macam::s_map_mdr(block_mvd, dist_w = dist_mvd, tp = 1,
                                          theta = theta_i, lambda = lambda_i,
                                          regularized = regularized_i,
                                          alpha = alpha_i, glmnet_parallel = TRUE,
                                          save_smap_coefficients = TRUE, random_seed = 1234),
                         silent = TRUE)
          if (class(mdr_res) == "try-error") {
            mdr_res <- list(stats = data.frame(N=NA, rho=NA, mae=NA, rmse=NA, pval=NA))
          } else {
            mdr_res$stats$pval <- cor.test(mdr_res$model_output$obs, mdr_res$model_output$pred, use = "complete.obs")$p.value
          }
        }

        # Collect results
        all_mdr_res_i <- data.frame(effect_var = tax_i, E = Ex)
        all_mdr_res_i$N <- mdr_res$stats$N
        all_mdr_res_i$rho <- mdr_res$stats$rho
        all_mdr_res_i$mae <- mdr_res$stats$mae
        all_mdr_res_i$rmse <- mdr_res$stats$rmse
        all_mdr_res_i$pval <- mdr_res$stats$pval
        all_mdr_res_i$theta <- theta_i
        all_mdr_res_i$lambda <- lambda_i
        all_mdr_res_i$alpha <- alpha_i
        all_mdr_res_i$regularized <- regularized_i
        all_mdr_res_i$block_mvd <- I(list(block_mvd))
        all_mdr_res_i$mvd_parms <- I(list(dist_mvd_list$parms))
        all_mdr_res_i$top_embeddings <- I(list(dist_mvd_list$top_embeddings))
        
        # Collect results
        if (step_i == 1) {
          all_mdr_res <- all_mdr_res_i
        } else {
          all_mdr_res <- rbind(all_mdr_res, all_mdr_res_i)
        }
        
        # Record time elapsed
        elapsed_time <- proc.time()[3] - start_time
        # Report progress
        message(sprintf("Step %s/%s, %s, theta = %s, lambda = %s, alpha = %s finished: %.1f sec elapsed",
                        step_i, total_step, tax_i, theta_i, lambda_i, alpha_i, elapsed_time))
        # Counting steps
        step_i <- step_i + 1
      }
    }
  }
  # Save temporal file
  saveRDS(all_mdr_res[all_mdr_res$effect_var == tax_i,],
            sprintf("%s/top%02d_%s.obj", output_folder, which(colnames(asv_df_NAadd_std) == tax_i), tax_i))
  write.csv(all_mdr_res[all_mdr_res$effect_var == tax_i,1:11],
            sprintf("%s/top%02d_%s.csv", output_folder, which(colnames(asv_df_NAadd_std) == tax_i), tax_i),
            row.names = FALSE)
}

# Delete temporal paramters
rm(simp_x); rm(theta_i); rm(alpha_i); rm(regularized_i); rm(lambda_i)
rm(step_i); rm(all_mdr_res_i); rm(mdr_res); rm(Ex); rm(env_i)
rm(start_time); rm(elapsed_time); rm(tax_i); rm(i)


# <---------------------------------------------> #
# Save all MDR S-map results
# <---------------------------------------------> #
saveRDS(all_mdr_res, sprintf("%s/all_mdr_res_sp.obj", output_folder))
all_mdr_res_min <- all_mdr_res[,1:11]
saveRDS(all_mdr_res_min, sprintf("%s/all_mdr_res_sp_min.obj", output_folder))
rm(all_mdr_res) # Delete this because this is too large


# <---------------------------------------------> #
# Choose optimal conditions for each taxon
# <---------------------------------------------> #
# Extract strongest causality (Environmental variables ==> eDNA)
all_mdr_best <- tibble()
for(effect_group in unique(all_mdr_res_min$effect_var)){
  mdr_tmp <- all_mdr_res_min %>% filter(effect_var == effect_group) %>% .[which.min(.$rmse),]
  all_mdr_best <- rbind(all_mdr_best, mdr_tmp); rm(mdr_tmp)
}
# Add taxa information
all_mdr_best_sp <- cbind(all_mdr_best, tax_df[all_mdr_best$effect_var,])
#hist(all_mdr_best_sp$rho)


# <---------------------------------------------> #
# Perform MDR S-map with optimal condition
# <---------------------------------------------> #
# Tentative parameters
for (i in 1:nrow(all_mdr_best_sp)) {
  # Record starting time
  start_time <- proc.time()[3]
  
  # Collect optimal parameters
  tax_i <- all_mdr_best_sp$effect_var[i]
  Ex <- all_mdr_best_sp$E[i]
  theta_i <- all_mdr_best_sp$theta[i]
  alpha_i <- all_mdr_best_sp$alpha[i]
  lambda_i <- all_mdr_best_sp$lambda[i]
  regularized_i <- all_mdr_best_sp$regularized[i]
  
  # Extract causal environmental variables
  env_df_tmp <- as.data.frame(edna_sp_summary[edna_sp_summary$effect_var == tax_i,14:17][1,])
  env_df_var <- sapply(str_split(env_df_tmp, pattern = ";"), `[`, 1)
  env_df_var <- env_df_var[!is.na(env_df_var)]
  env_df_tp <- sapply(str_split(env_df_tmp, pattern = "; tp = "), `[`, 2)
  env_df_tp <- as.numeric(env_df_tp[!is.na(env_df_tp)])
  if (length(env_df_var) > 0) {
    for (env_i in 1:length(env_df_var)) {
      if (env_i == 1) {
        env_block_mvd <- as.data.frame(dplyr::lag(sample_df_NAadd[,env_df_var[env_i]], n = abs(env_df_tp[env_i])))
      } else {
        env_block_tmp <- as.data.frame(dplyr::lag(sample_df_NAadd[,env_df_var[env_i]], n = abs(env_df_tp[env_i])))
        env_block_mvd <- cbind(env_block_mvd, env_block_tmp)
      }
    }
    colnames(env_block_mvd) <- sprintf("%s_tp%s", env_df_var, env_df_tp)
  }
  
  # Make MVD block (only strongest causality is included)
  block_mvd <- macam::make_block_mvd(asv_df_NAadd_std,
                                     as.data.frame(edna_sp_summary[edna_sp_summary$effect_var == tax_i,]),
                                     effect_var = tax_i,
                                     E_effect_var = Ex,
                                     cause_var_colname = "cause_var",
                                     include_var = "strongest_only",
                                     p_threshold = 0.05,
                                     sort_tp = T, silent = TRUE)
  if (length(env_df_var) > 0) {
    # Merge environmental variables and ASV block
    block_mvd <- cbind(block_mvd, env_block_mvd)
  }
  
  # Compute multiview distance
  dist_mvd_list <- macam::compute_mvd(block_mvd, effect_var = tax_i, Ex, tp = 1,
                                      make_block_method = "naive",
                                      n_ssr = 1000, k = floor(sqrt(1000)),
                                      random_seed = 1234, distance_only = FALSE,
                                      silent = TRUE)
  dist_mvd <- dist_mvd_list$multiview_dist
  
  # Perform MDR S-map
  mdr_res <- try(macam::s_map_mdr(block_mvd, dist_w = dist_mvd, tp = 1,
                                  theta = theta_i, lambda = lambda_i,
                                  regularized = regularized_i,
                                  alpha = alpha_i, glmnet_parallel = TRUE,
                                  save_smap_coefficients = TRUE, random_seed = 1234),
                 silent = TRUE)
  if (class(mdr_res) == "try-error") {
    mdr_res <- list(stats = data.frame(N=NA, rho=NA, mae=NA, rmse=NA))
  } else {
    mdr_res$stats$pval <- cor.test(mdr_res$model_output$obs, mdr_res$model_output$pred, use = "complete.obs")$p.value
  }
  
  # Collect results
  opt_mdr_res_i <- data.frame(effect_var = tax_i, E = Ex)
  opt_mdr_res_i$N <- mdr_res$stats$N
  opt_mdr_res_i$rho <- mdr_res$stats$rho
  opt_mdr_res_i$mae <- mdr_res$stats$mae
  opt_mdr_res_i$rmse <- mdr_res$stats$rmse
  opt_mdr_res_i$pval <- mdr_res$stats$pval
  opt_mdr_res_i$theta <- theta_i
  opt_mdr_res_i$lambda <- lambda_i
  opt_mdr_res_i$alpha <- alpha_i
  opt_mdr_res_i$regularized <- regularized_i
  opt_mdr_res_i$block_mvd <- I(list(block_mvd))
  opt_mdr_res_i$embeddings <- I(list(dist_mvd_list$embeddings))
  #opt_mdr_res_i$multiview_dist <- I(list(dist_mvd_list$multiview_dist))
  opt_mdr_res_i$parms <- I(list(dist_mvd_list$parms))
  opt_mdr_res_i$model_output <- I(list(mdr_res$model_output))
  opt_mdr_res_i$smap_coefficients <- I(list(mdr_res$smap_coefficients))
  
  # Collect results
  if (i == 1) { opt_mdr_res <- opt_mdr_res_i } else { opt_mdr_res <- rbind(opt_mdr_res, opt_mdr_res_i)}
  
  # Record time elapsed
  elapsed_time <- proc.time()[3] - start_time
  # Report progress
  message(sprintf("Step %s/%s, %s, theta = %s, lambda = %s, alpha = %s finished: %.1f sec elapsed",
                  i, nrow(all_mdr_best_sp), tax_i, theta_i, lambda_i, alpha_i, elapsed_time))
}
# Delete temporal objects
rm(i); rm(opt_mdr_res_i); rm(mdr_res)
rm(tax_i); rm(theta_i); rm(regularized_i); rm(alpha_i)
rm(lambda_i); rm(Ex); rm(dist_mvd); rm(dist_mvd_list)
rm(block_mvd); rm(asv_df_tmp); rm(sample_df_tmp)
rm(env_i); rm(env_block_tmp); rm(env_block_mvd)
rm(env_df_tmp); rm(NA_block_asv); rm(NA_block_sample)
rm(NA_nrow); rm(start_time); rm(elapsed_time)
rm(total_step); rm(ran_seed)


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save the MDR S-map results
saveRDS(opt_mdr_res, sprintf("%s/opt_mdr_res.obj", output_folder))

# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))

