####
#### Boso Peninsula project
#### No. 4 Compile UIC results
####

# Load workspace
load("01_CompileDataOut/01_CompileDataOut.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.1.7
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
theme_set(theme_cowplot())

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)
output_folder_load <- "03_UICEnvDNAOut"


# <---------------------------------------------> #
#                Load UIC results                 #
# <---------------------------------------------> #
top_taxa <- readRDS(sprintf("%s/top_taxa.obj", output_folder_load))
uic_edna_env <- readRDS(sprintf("%s/uic_edna_env.obj", output_folder_load))
uic_env_edna <- readRDS(sprintf("%s/uic_env_edna.obj", output_folder_load))
uic_edna_edna_sp <- readRDS(sprintf("%s/uic_edna_edna_sp.obj", output_folder_load))
uic_edna_edna_st <- readRDS(sprintf("%s/uic_edna_edna_site.obj", output_folder_load))


# <---------------------------------------------> #
#       Extract significant UIC results           #
# <---------------------------------------------> #
# Extract causal pairs
## uic_signif01: from uic_edna_env
## uic_signif02: from uic_env_edna
## uic_signif03: from uic_edna_edna_sp
## uic_signif04: from uic_edna_edna_st
pval_threshold = 0.05
## Influence of environmental variables to eDNA
uic_signif01 <- uic_edna_env %>%
  filter(pval <= pval_threshold) %>%
  filter(tp <= 0) %>% mutate(cause_pair = paste0(effect_var, "_CausedBy_", cause_var))
dim(uic_signif01); hist(uic_signif01$te)
## Influence of eDNA to eDNA environemtal variables
uic_signif02 <- uic_env_edna %>%
  filter(pval <= pval_threshold) %>%
  filter(tp <= 0) %>% mutate(cause_pair = paste0(effect_var, "_CausedBy_", cause_var))
dim(uic_signif02); hist(uic_signif02$te)
## Interactions between eDNA species
uic_signif03 <- uic_edna_edna_sp %>%
  filter(pval <= pval_threshold) %>%
  filter(tp <= 0) %>% mutate(cause_pair = paste0(effect_var, "_CausedBy_", cause_var))
dim(uic_signif03); hist(uic_signif03$te)
hist(uic_edna_edna_sp$te, breaks = seq(-1,3,by = 0.1))
## Interactions between site-by-site eDNA species
uic_signif04 <- uic_edna_edna_st %>%
  filter(pval <= pval_threshold) %>%
  filter(tp <= 0) %>% mutate(cause_pair = paste0(effect_var, "_CausedBy_", cause_var))
dim(uic_signif04); hist(uic_signif04$te)
hist(uic_edna_edna_st$te)


# <---------------------------------------------> #
#       Extract strongest UIC results             #
# <---------------------------------------------> #
# Extract strongest causality (Environmental variables ==> eDNA)
uic_strong01 <- tibble()
for(cause_group in unique(uic_signif01$cause_pair)){
  uic_tmp <- uic_signif01 %>% filter(cause_pair == cause_group) %>% .[which.max(.$te),]
  uic_strong01 <- rbind(uic_strong01, uic_tmp)
  rm(uic_tmp)
}
edna_env_summary <- cbind(uic_strong01, tax_df[uic_strong01$effect_var,])

# Extract strongest causality (eDNA ==> Environmental variables)
uic_strong02 <- tibble()
for(cause_group in unique(uic_signif02$cause_pair)){
  uic_tmp <- uic_signif02 %>% filter(cause_pair == cause_group) %>% .[which.max(.$te),]
  uic_strong02 <- rbind(uic_strong02, uic_tmp)
  rm(uic_tmp)
}
env_edna_summary <- cbind(uic_strong02, tax_df[uic_strong02$cause_var,])

# Extract strongest causality (eDNA <==> eDNA; site-averaged)
uic_strong03 <- tibble()
for(cause_group in unique(uic_signif03$cause_pair)){
  uic_tmp <- uic_signif03 %>% filter(cause_pair == cause_group) %>% .[which.max(.$te),]
  uic_strong03 <- rbind(uic_strong03, uic_tmp)
  rm(uic_tmp)
}
edna_sp_summary <- cbind(uic_strong03[,1],
                         tax_df[uic_strong03$effect_var,"common_jp_name"],
                         uic_strong03[,2],
                         tax_df[uic_strong03$cause_var,"common_jp_name"],
                         uic_strong03[,3:ncol(uic_strong03)])
colnames(edna_sp_summary)[c(1,3)] <- c("effect_var", "cause_var")
colnames(edna_sp_summary)[c(2,4)] <- c("effect_var_jp", "cause_var_jp")
## Remove rows with cause_var == effect_var
edna_sp_summary <- edna_sp_summary %>% filter(cause_var != effect_var)

# Extract strongest causality (eDNA <==> eDNA; site-by-site)
uic_strong04 <- tibble()
for(cause_group in unique(uic_signif04$cause_pair)){
  uic_tmp <- uic_signif04 %>% filter(cause_pair == cause_group) %>% .[which.max(.$te),]
  uic_strong04 <- rbind(uic_strong04, uic_tmp)
  rm(uic_tmp)
}
edna_st_summary <- cbind(uic_strong04[,1],
                         tax_df[str_sub(uic_strong04$effect_var, end = -6),"common_jp_name"],
                         uic_strong04[,2],
                         tax_df[str_sub(uic_strong04$cause_var, end = -6),"common_jp_name"],
                         uic_strong04[,3:ncol(uic_strong03)])
colnames(edna_st_summary)[c(1,3)] <- c("effect_var", "cause_var")
colnames(edna_st_summary)[c(2,4)] <- c("effect_var_jp", "cause_var_jp")
## Remove rows with cause_var == effect_var
edna_st_summary <- edna_st_summary %>% filter(cause_var != effect_var)


# <---------------------------------------------> #
#                  Save UIC results               #
# <---------------------------------------------> #
# Save results
write_excel_csv(edna_env_summary, sprintf("%s/uic_edna_env_strong.csv", output_folder))
write_excel_csv(env_edna_summary, sprintf("%s/uic_env_edna_strong.csv", output_folder))
write_excel_csv(edna_sp_summary, sprintf("%s/uic_edna_sp_strong.csv", output_folder))
write_excel_csv(edna_st_summary, sprintf("%s/uic_edna_st_strong.csv", output_folder))

# Delete temporal objects
rm(edna_env_summary); rm(env_edna_summary); rm(edna_sp_summary); rm(edna_st_summary)
rm(uic_strong01); rm(uic_strong02); rm(uic_strong03); rm(uic_strong04)
rm(uic_signif01); rm(uic_signif02); rm(uic_signif03); rm(uic_signif04)
rm(output_folder_load); rm(pval_threshold)

# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))

