####
#### Boso Peninsula project
#### No.2: Visualize fish interaction network (all sites combined)
#### 2021.10.13 Ushio
#### 2021.11.10 Ushio
#### 2022.11.10 Ushio
####

# Load workspace
load("../07_VisualizeNetworkOut/07_VisualizeNetworkOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.11.10
library(cols4all); packageVersion("cols4all") # 0.4, 2022.11.10, c4a_gui()
library(ggsci); packageVersion("ggsci") # 2.9, 2021.1.7
library(ggraph); packageVersion("ggraph") # 2.1.0, 2022.11.10
library(igraph); packageVersion("igraph") # 1.3.5, 2022.11.10
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"

# Prepare color palette
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
get_palette2 <- colorRampPalette(brewer.pal(8, "Set1"))
palette_custom <- rev(get_palette(11))
palette_custom2 <- get_palette(28)
palette_custom3 <- get_palette2(28)
palette_c4a <- c4a("palette36", 28)
#c4a_gui()
#"palette36"
#"polychrome36"

#------------------------------------------#
# Visualize each fish dynamics
#------------------------------------------#
# Original data
#taxa_vertices
# Select color palettes
#c4a_gui()

# Convert jp name to en name
tax_tbl$scientific_name
tax_tbl$common_jp_name
jp_name_id <- match(node_coord$species, tax_tbl$common_jp_name)
node_coord$species_sci <- tax_tbl$scientific_name[jp_name_id]

# Leverage all site information
gn3 <- gn1 +
  # Set edge parameters
  geom_conn_bundle(data = get_con(from = edge_from, to = edge_to, value = edge_strength),
                   aes(x = x * 0.95, y = y * 0.95, color = value), edge_alpha = 0.8,
                   arrow = grid::arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")) +
  # Set node parameters
  geom_node_point(aes(filter = leaf, x = x*1, y = y*1, size = abundance, fill = group3, shape = 21),
                  color = "gray10") +
  # Blank node parameters to enlarge plotting area (Otherwise part of node labels might disappear)
  geom_node_point(aes(filter = leaf, x = x*1.5, y = y*1.5), size = NA, shape = NA) +
  # Add label information
  annotate("text",
           x = node_coord$x * 1.1,
           y = node_coord$y * 1.1,
           label = as.character(node_coord$species_sci),
           family = "Arial",
           fontface = "italic",
           size = 2, angle = node_coord$label_angle, hjust = node_coord$hjust_val) +
  # Set color parameters
  scale_size_binned(limits = c(1,300), breaks = c(10,50,100,200)) +
  scale_fill_manual(values = palette_c4a) +
  guides(fill = guide_legend(override.aes = list(size = 3, shape = 21,
                                                 fill = palette_c4a))) +
  scale_shape_identity() +
  scale_edge_color_gradient2(low = "azure", mid = "royalblue", high = "darkorchid4",
                             midpoint = 6, limits = c(0,10)) + 
  # Set titles
  ggtitle("Fish interaction network (All sites, only strong interactions are shown)") +
  labs(fill = "Family", size = "Abudance", edge_color = "Transfer Entropy (scaled)") +
  # Set theme parameters
  theme_graph(title_family = "Arial", title_margin = 5, title_size = 10) +
  theme(text = element_text(family = "Arial"),
        legend.position="right",
        legend.key.height = grid::unit(0.6, "cm"),
        legend.key.width = grid::unit(0.6, "cm"),
        legend.text = element_text(size = 7),
        plot.title = element_text(vjust = 8)) +
  #theme(plot.background = element_rect("white")) +
  NULL
#gn3
ggsave(file = sprintf("%s/JPG_FishNetworkAll.jpg", fig_folder1),
       plot = gn3, width = 11.5, height = 8, dpi = 450)
ggsave(file = sprintf("%s/JPG_FishNetworkAll2.jpg", fig_folder1),
       plot = gn3 + theme(legend.position = "none") + ggtitle(NULL),
       width = 10, height = 10, dpi = 600)


# Visualize statistics
p1 <- ggplot(tax_tbl, aes(x = habitat, y = intcat_te_sum)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 <- ggplot(tax_tbl, aes(x = water_area, y = intcat_te_sum)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2)
p3 <- ggplot(tax_tbl, aes(x = edna_sum, y = intcat_te_sum)) + geom_point()
p4 <- ggplot(tax_tbl, aes(x = edna_sum, y = int_n)) + geom_point()
p5 <- ggplot(tax_tbl, aes(x = int_n, y = intcat_te_sum)) + geom_point()
p_all <- plot_grid(p1, p2, p3, p4, p5, ncol = 3, align = "hv", axis = "lrtb")
p_list <- list(p1, p2, p3, p4, p5)

# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save figures
quartz(file = sprintf("%s/PDF_FishNetworkAll.pdf", fig_folder1),
       type = "pdf", family = "Arial", width = 11.5, height = 8)
gn3; dev.off()
ggsave(file = sprintf("%s/PDF_FishNetworkProperty.pdf", fig_folder1),
       plot = p_all, width = 15, height = 12.5)

# Save R objects
saveRDS(gn3, sprintf("%s/Fig_NetworkAll.obj", fig_folder1))
saveRDS(p_list, sprintf("%s/Fig_NetworkProp.obj", fig_folder1))
#readRDS(sprintf("%s/Fig_NetworkAll.obj", fig_folder1))

# Save supplementary materials
## Output TE matrix for check
write.csv(te_mat, sprintf("%s/CHECK_TE_matrix.csv", fig_folder1))
write.csv(tax_tbl, sprintf("%s/CHECK_TaxTable.csv", fig_folder1))
te_mat2 <- te_mat
te_mat2_names <- as.character(unlist(tax_tbl[match(rownames(te_mat2), tax_tbl$tax_id),"scientific_name"]))
colnames(te_mat2) <- rownames(te_mat2) <- te_mat2_names
write.csv(te_mat2, sprintf("%s/CHECK_TaxTable2.csv", fig_folder1))

