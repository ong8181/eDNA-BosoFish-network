####
#### Boso Peninsula project
#### No. 6 Visualize UIC results site-by-site
####

# Load workspace
load("03_UICEnvDNAOut/03_4_UICDNAxDNA_SpSiteOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.11.9
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
library(ggraph); packageVersion("ggraph") # 2.0.5, 2021.3.1
library(igraph); packageVersion("igraph") # 1.2.8, 2021.11.9
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

## Set color here
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.10.13
get_palette2 <- colorRampPalette(brewer.pal(8, "Set1"))
palette_custom3 <- get_palette2(28)

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)
output_folder_load <- "04_CompileUICresOut"

# <---------------------------------------------> #
#                Load UIC results                 #
# <---------------------------------------------> #
uic_edna_env_strong <- read_csv(sprintf("%s/uic_edna_env_strong.csv", output_folder_load))
uic_env_edna_strong <- read_csv(sprintf("%s/uic_env_edna_strong.csv", output_folder_load))
uic_edna_st_strong <- read_csv(sprintf("%s/uic_edna_st_strong.csv", output_folder_load),
                               col_types = cols(cond1 = col_character(),
                                                cond2 = col_character(),
                                                cond3 = col_character(),
                                                cond4 = col_character()))

# <---------------------------------------------> #
#            Prepare node information
# <---------------------------------------------> #
# Select top taxa that was used for UIC
tax_df[top_taxa,]; asv_df_sitesp
# Prepare tax data (+ change character to factor)
tax_tbl <- as_tibble(tax_df[top_taxa,])
tax_tbl <- tax_tbl %>% mutate_if(is.character, as.factor)
tax_tbl$tax_id <- rownames(tax_df[top_taxa,])


# <---------------------------------------------> #
#            Compile site-specific data
# <---------------------------------------------> #
#suspend_time <- 10
#for(site_i in unique(sample_df$site_code)){
#site_list <- "St01"
#for(site_i in site_list){
  # Choose site
  site_i <- "St11"
  
  # Select site-specific data
  asv_df_site_i <- asv_df_sitesp[,grep(site_i, colnames(asv_df_sitesp))]
  # Rename colnames
  colnames(asv_df_site_i) <- sapply(str_split(colnames(asv_df_site_i), "_"), "[", 1)
  
  # Site-specific UIC res
  effect_var_site_i <- which(sapply(str_split(uic_edna_st_strong$effect_var, "_"), "[", 2) == site_i)
  cause_var_site_i <- which(sapply(str_split(uic_edna_st_strong$cause_var, "_"), "[", 2) == site_i)
  uic_edna_st_i <- uic_edna_st_strong[intersect(effect_var_site_i, cause_var_site_i),]
  
  # Remove absent species
  absent_sp_id <- names(which(!(colSums(asv_df_site_i) > 0)))
  absent_sp_effect_i <- which(!is.na(match(str_sub(uic_edna_st_i$effect_var, end = -6), absent_sp_id)))
  absent_sp_cause_i <- which(!is.na(match(str_sub(uic_edna_st_i$cause_var, end = -6), absent_sp_id)))
  absent_sp_id <- unique(absent_sp_effect_i, absent_sp_cause_i)
  if(length(absent_sp_id) > 0) uic_edna_st_i <- uic_edna_st_i[-absent_sp_id,]
  
  
  # <---------------------------------------------> #
  #           Visualize information flow            #
  #                  between species                #
  # <---------------------------------------------> #
  # Prepare taxa information
  tax_tbl$order_family <- paste0(tax_tbl$order, "_", tax_tbl$family)
  network_taxa_name <- sort(as.character(unique(tax_tbl$order_family)))
  taxa_edge <- data.frame(from = "origin", to = network_taxa_name)
  for(taxa_name_i in network_taxa_name){
    to_i <- tax_tbl[tax_tbl$order_family == taxa_name_i,]$tax_id
    taxa_edge <- rbind(taxa_edge, data.frame(from = taxa_name_i, to = to_i))
  }
  # Create a vertices data frame. One line per object of our hierarchy
  taxa_vertices <- data.frame(name = c("origin", network_taxa_name),
                              abundance = NA,
                              log_abundance = NA,
                              group1 = c(NA, rep("origin", length(network_taxa_name))),
                              group2 = c(NA, rep("origin", length(network_taxa_name))),
                              group3 = c(NA, rep("origin", length(network_taxa_name))),
                              species = NA) %>%
    rbind(data.frame(name = tax_tbl$tax_id,
                     abundance = ceiling(colSums(asv_df_site_i[,top_taxa], na.rm = T)),
                     log_abundance = log(ceiling(colSums(asv_df_site_i[,top_taxa] + 1, na.rm = T))),
                     group1 = tax_tbl$family,
                     group2 = tax_tbl$order,
                     group3 = tax_tbl$order_family,
                     species = as.character(tax_tbl$common_jp_name)))
  taxa_vertices$group1 <- as.factor(taxa_vertices$group1)
  taxa_vertices$group2 <- as.factor(taxa_vertices$group2)
  taxa_vertices$group3 <- as.factor(taxa_vertices$group3)
  taxa_vertices$node_shape <- as.integer(21)
  taxa_vertices$node_shape[which(taxa_vertices$abundance == 0)] <- as.integer(4)
  
  # Add information of the label we are going to add: angle, horizontal adjustment and potential flip
  taxa_vertices$id <- NA
  my_leaves <- which(is.na(match(taxa_vertices$name, taxa_edge$from)))
  n_leaves <- length(my_leaves)
  taxa_vertices$id[my_leaves] <- seq(1:n_leaves)
  taxa_vertices$angle <- 90 - 360 * (taxa_vertices$id) / n_leaves
  taxa_vertices$hjust <- 1

  # Create network matrix
  ## Prepare empty matrix
  te_mat <- matrix(rep(0, length(unique(tax_tbl$tax_id))^2),
                   ncol = length(unique(tax_tbl$tax_id)))
  colnames(te_mat) <- rownames(te_mat) <- tax_tbl$tax_id
  ## Assigning TE values to the empty matrix
  ### Rows = effect_var, Cols = cause_var
  for(i in 1:nrow(uic_edna_st_i)){
    effect_var_i <- str_sub(as.character(uic_edna_st_i[i,"effect_var"]), end = -6)
    cause_var_i <- str_sub(as.character(uic_edna_st_i[i,"cause_var"]), end = -6)
    te_mat[effect_var_i,cause_var_i] <- as.numeric(uic_edna_st_i[i,"te"])
  }
  ### Apply threshold for network visualization
  te_threshold <- quantile(c(te_mat[te_mat != 0]), probs = seq(0,1,0.1))
  te_mat[te_mat < te_threshold['80%']] <- 0
  ### Sort te_mat
  te_mat <- te_mat[rownames(taxa_vertices)[30:79], rownames(taxa_vertices)[30:79]]
  
  # Creat edge list
  edge_list <- as.data.frame(get.edgelist(graph.adjacency(te_mat, weighted = TRUE, diag = FALSE, mode = "directed")))
  colnames(edge_list) <- c("from", "to")
  edge_list$value <- c(te_mat[te_mat != 0])
  
  # Create a plot
  edge_graph <- igraph::graph_from_data_frame(taxa_edge, vertices = taxa_vertices)
  edge_from <- match(edge_list$from, taxa_vertices$name)
  edge_to <- match(edge_list$to, taxa_vertices$name)
  edge_strength <- as.numeric(unlist(infotheo::discretize(abs(edge_list$value), nbins = 10)))

  # Generate network plot
  gn1 <- ggraph(edge_graph, layout = 'dendrogram', circular = TRUE)
  
  # Compile node labels information
  node_coord <- data.frame(x = gn1$data$x, y = gn1$data$y,
                           order_family = gn1$data$group3,
                           species = gn1$data$species,
                           family = gn1$data$group1)
  node_coord <- node_coord[!is.na(gn1$data$species),]
  node_coord$label_angle <- node_coord$hjust_val <- NaN
  for(id_i in 1:nrow(node_coord)){
    xy_values <- as.numeric(unlist(node_coord[id_i, c("x","y")]))
    # Calculate angle and hjust values
    node_coord$label_angle[id_i] <- 90 - atan2(xy_values[1], xy_values[2])/(2*pi)*360
    node_coord$hjust_val[id_i] <- 0
    # Rotate if x coordinate < 0
    if(xy_values[1] < 0){
      node_coord$label_angle[id_i] <- node_coord$label_angle[id_i] + 180
      node_coord$hjust_val[id_i] <- 1
    }
  }
  
  # Customize netowrk plot
  gn2 <- gn1 +
    # Set edge parameters
    geom_conn_bundle(data = get_con(from = edge_from, to = edge_to, value = edge_strength),
                     aes(x = x * 0.95, y = y * 0.95, color = value), edge_alpha = 0.8,
                     arrow = grid::arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")) +
    # Set node parameters
    geom_node_point(aes(filter = leaf, x = x*1, y = y*1,
                        size = abundance, fill = group3, shape = node_shape),
                    color = "gray10") +
    # Blank node parameters to enlarge plotting area (Otherwise part of node labels might disappear)
    geom_node_point(aes(filter = leaf, x = x*1.5, y = y*1.5), size = NA, shape = NA) +
    # Add label information
    annotate("text",
             x = node_coord$x * 1.1,
             y = node_coord$y * 1.1,
             label = as.character(node_coord$species),
             family = "HiraKakuProN-W3",
             size = 2, angle = node_coord$label_angle, hjust = node_coord$hjust_val) +
    # Set color parameters
    #scale_fill_igv() +
    scale_fill_manual(values = palette_custom3) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_identity() +
    scale_edge_color_gradient2(low = "azure", mid = "royalblue", high = "darkorchid4", midpoint = 6) + 
    # Set titles
    #ggtitle(sprintf("房総半島 %s (%s)・魚類群集ネットワーク (弱い因果は除去)",
    ggtitle(sprintf("Boso Peninsula %s (%s), fish-fish interaction network (Only strong interactions are shown)",
                    site_i, unique(sample_df[sample_df$site_code == site_i, "site_name"]))) +
    labs(fill = "Order_Family", size = "Abudance", edge_color = "Interaction strength") +
    # Set theme parameters
    theme_graph(title_family = "HiraKakuProN-W3", title_margin = 5, title_size = 10) +
    theme(text = element_text(family = "HiraKakuProN-W3"),
          legend.position="right",
          legend.key.height = grid::unit(0.6, "cm"),
          legend.key.width = grid::unit(0.6, "cm"),
          legend.text = element_text(size = 7),
          plot.title = element_text(vjust = 8)) +
    NULL
  #gn2
  #Sys.sleep(suspend_time)
  
  quartz(file = sprintf("%s/FishNetwork_%s.jpg", output_folder, site_i),
         type = "jpg", family = "HiraKakuPro-W3", width = 11.5, height = 7.5, dpi = 300)
  gn2; dev.off()
  #Sys.sleep(suspend_time)
#}


# <---------------------------------------------> #
#                    Save results                 #
# <---------------------------------------------> #
# Save results
# Delete large files
rm(uic_edna_edna_site0)
# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))

