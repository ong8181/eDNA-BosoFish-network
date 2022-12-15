####
#### Boso Peninsula project
#### Format figures
#### 2021.10.13 Ushio
#### 2021.10.29 Ushio
#### 2021.11.10 Ushio
#### 2022.04.01 Ushio
#### 2022.05.17 Ushio (R4.1.2)
#### 2022.11.11 Ushio (R4.2.1)
#### 2022.12.14 Ushio (R4.2.1)
####

#setwd("FigCode/")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.11
library(lubridate); packageVersion("lubridate") # 1.8.0, 2021.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggimage); packageVersion("ggimage") # 0.3.1, 2022.11.11
library(magick); packageVersion("magick") # 2.7.3, 2021.12.8
library(igraph); packageVersion("igraph") # 1.3.5, 2022.11.11
library(ggraph); packageVersion("ggraph") # 2.1.0, 2022.11.11
library(ggpubr); packageVersion("ggpubr") # 0.4.0, 2021.10.17
#library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
#library(patchwork); packageVersion("patchwork") # 1.1.2, 2022.11.14
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"

# Prepare color palette
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.11.11
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
palette_custom <- get_palette(11)[c(1:4,7:11,6:5)]


# <---------------------------------------------> #
#  Load figure objects
# <---------------------------------------------> #
# Site map
Fig_SiteMap <- image_read("0_Illustrations/SiteMap.png")
Fig_SiteMap <- ggdraw() + draw_image(Fig_SiteMap)
# General pattern
Fig_DNApattern01 <- readRDS(paste0(fig_folder1, "/Fig_DNApattern01.obj"))
Fig_DNApattern02 <- readRDS(paste0(fig_folder1, "/Fig_DNApattern02.obj"))
Fig_DNApatternCombined <- readRDS(paste0(fig_folder1, "/Fig_DNApatternCombined.obj"))
Fig_CorPairs <- readRDS(paste0(fig_folder1, "/Fig_CorPairs.obj"))

# UIC-based network structure
#Fig_NetworkAll_raw <- readRDS(paste0(fig_folder1, "/Fig_NetworkAll.obj")) # Vector version
Fig_NetworkAll <- image_read("0_Illustrations/FishNetworkAll_pls_Fish2.jpg") # Image version
Fig_NetworkAll <- ggdraw() + draw_image(Fig_NetworkAll)
Fig_NetworkProp <- readRDS(paste0(fig_folder1, "/Fig_NetworkProp.obj"))

# MDR S-map: Overall patterns, the effects of environmental variables
Fig_IS_EnvAll1 <- readRDS(paste0(fig_folder1, "/Fig_OverallEffect1.obj"))
Fig_IS_EnvAll2 <- readRDS(paste0(fig_folder1, "/Fig_OverallEffect2.obj"))

# MDR S-map: Overall patterns, all data points
Fig_ISoverall1 <- image_read(paste0(fig_folder1, "/JPG_Overall_pattern_SI1.jpg"))
Fig_ISoverall1 <- ggdraw() + draw_image(Fig_ISoverall1)
Fig_ISoverall2 <- image_read(paste0(fig_folder1, "/JPG_Overall_pattern_SI2.jpg"))
Fig_ISoverall2 <- ggdraw() + draw_image(Fig_ISoverall2)

# MDR S-map: Interaction strength for each fish species
Fig_ISvsDNA <- readRDS(paste0(fig_folder1, "/Fig_ISvsDNA.obj"))
Fig_ISvsRich <- readRDS(paste0(fig_folder1, "/Fig_ISvsRich.obj"))
Fig_ISvsTemp <- readRDS(paste0(fig_folder1, "/Fig_ISvsTemp.obj"))

# Load individual fish images
## Fish images will be in "ggfish" object.
#source("F05_LoadFishImages.R")
#saveRDS(ggfish, sprintf("%s/ggfish.obj", fig_folder1))
# Put fish images on the vector version networks
#Fig_NetworkAll +
#  patchwork::inset_element(ggfish$Acanthopagrus_schelegeli, 0.5, 0.5, 0.65, 0.65)


# <---------------------------------------------> #
#  Compile figures
# <---------------------------------------------> #
# Figure: Overall pattern
Fig_SiteDynamics <- plot_grid(Fig_SiteMap,
                              Fig_DNApatternCombined,
                              ncol = 1, rel_heights = c(1, 2.5),
                              labels = c("a", NA))

# Figure: Network figure
# Load network figures
# net_legend_all <- get_legend(Fig_NetworkAll_raw + theme(legend.text = element_text(size = 10)))
# net_legend1 <- ggpubr::as_ggplot(net_legend_all$grobs[[1]])
# net_legend2 <- ggpubr::as_ggplot(net_legend_all$grobs[[2]])
# net_legend3 <- ggpubr::as_ggplot(net_legend_all$grobs[[3]])
# net_legend <- plot_grid(NULL, plot_grid(net_legend1, net_legend2, NULL, ncol = 3),
#                         net_legend3, NULL, ncol = 1, rel_heights = c(1, 1.2, 3, 1))
Fig_NetworkAll


# Figure: IS overall petterns
## Prepare sub-titles
title_inst <- ggdraw() + draw_label("In-strength", fontface = 'bold', size = 15, x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
title_otst <- ggdraw() + draw_label("Out-strength", fontface = 'bold', size = 15, x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
## Assemble panels
Fig_ISenv_Main1 <- plot_grid(Fig_IS_EnvAll1[[1]], Fig_IS_EnvAll1[[2]], Fig_IS_EnvAll1[[3]],
                            align = "hv", nrow = 1, labels = c("a", "b", "c"))
Fig_ISenv_Main2 <- plot_grid(Fig_IS_EnvAll1[[4]], Fig_IS_EnvAll1[[5]], Fig_IS_EnvAll1[[6]],
                             align = "hv", nrow = 1, labels = c("d", "e", "f"))
Fig_ISenv_Main <- plot_grid(title_inst, Fig_ISenv_Main1,
                            title_otst, Fig_ISenv_Main2,
                            nrow = 4, rel_heights = c(0.2, 1, 0.2, 1))
Fig_ISenv_SI1 <- plot_grid(Fig_IS_EnvAll2[[1]], Fig_IS_EnvAll2[[2]], Fig_IS_EnvAll2[[3]],
                           align = "hv", nrow = 1, labels = c("a", "b", "c"))
Fig_ISenv_SI2 <- plot_grid(Fig_IS_EnvAll2[[4]], Fig_IS_EnvAll2[[5]], Fig_IS_EnvAll2[[6]],
                           align = "hv", nrow = 1, labels = c("d", "e", "f"))
Fig_ISenv_SI <- plot_grid(title_inst, Fig_ISenv_SI1,
                          title_otst, Fig_ISenv_SI2,
                          nrow = 4, rel_heights = c(0.2, 1, 0.2, 1))

# Figure
## Interaction strength + temperature
legend01 <- get_legend(Fig_ISvsTemp[[1]])
Fig_ISvsTemp_merge <- plot_grid(Fig_ISvsTemp[[1]] + ggtitle("In-strength") + ylab("Interaction strength") + theme(legend.position = "none"),
                                Fig_ISvsTemp[[2]] + ggtitle("Out-strength") + ylab("Interaction strength") + theme(legend.position = "none"),
                                legend01,
                                ncol = 1, rel_heights = c(3.4,1,0.3), labels = c("a", "b"))
#ggsave(sprintf("%s/PDF_ISvsTemp_merge.pdf", fig_folder1), Fig_ISvsTemp_merge,
#       width = 12, height = 12)


# <---------------------------------------------> #
#  Supplementary figures
# <---------------------------------------------> #
# Japanese black seabream eDNA copies
Fig_qPCR <- Fig_DNApattern01[[1]] +
  scale_color_manual(name = "Site ID", values = palette_custom) +
  theme(panel.grid.major.x = element_line(linetype = 2, linewidth = 0.4, color = "gray70")) +
  #scale_y_continuous(trans = "log10", labels = label_func, limits = c(0.2, 50)) +
  NULL

# Temperature, eDNA copies, and species richness
site_legend <- get_legend(Fig_DNApattern02[[1]])
Fig_Scatter <- plot_grid(Fig_DNApattern02[[1]] + theme(legend.position = "none"),
                         Fig_DNApattern02[[2]] + theme(legend.position = "none"),
                         site_legend, rel_widths = c(1,1,0.2), axis = "tblr",
                         labels = c("b","c", NULL), align = "hv", ncol = 3)
Fig_Pattern2 <- plot_grid(Fig_DNApattern02[[3]] +
                            theme(panel.grid.major.x = element_line(linetype = 2, color = "gray70"),
                                  panel.grid.minor.x = element_line(linetype = 2, color = "gray70")),
                          Fig_Scatter,
                          labels = c("a", NULL), ncol = 1)

# Correlation plot
#Fig_CorPairs

# Other parameters
## IS v.s. total DNA concentrations
Fig_ISvsDNA_merge <- plot_grid(Fig_ISvsDNA[[1]] + ggtitle("In-strength") + ylab("Interaction strength") + theme(legend.position = "none"),
                               plot_grid(Fig_ISvsDNA[[2]] + ggtitle("Out-strength") + ylab("Interaction strength") + theme(legend.position = "none"),
                                         NULL, legend01, nrow = 1, rel_widths = c(3.1,0.3,1.7)),
                               ncol = 1, rel_heights = c(3.4,1,0.3), labels = c("a", "b", NA))
# ggsave(sprintf("%s/PDF_ISvsDNA_merge.pdf", fig_folder1), Fig_ISvsDNA_merge,
#        width = 12, height = 12)

## IS v.s. species richness
Fig_ISvsRich_merge <- plot_grid(Fig_ISvsRich[[1]] + ggtitle("In-strength") + ylab("Interaction strength") + theme(legend.position = "none"),
                               plot_grid(Fig_ISvsRich[[2]] + ggtitle("Out-strength") + ylab("Interaction strength") + theme(legend.position = "none"),
                                         NULL, legend01, nrow = 1, rel_widths = c(1.15, 0.3,3.7)),
                               ncol = 1, rel_heights = c(2.5,1,0.3), labels = c("a", "b"))
# ggsave(sprintf("%s/PDF_ISvsRich_merge.pdf", fig_folder1), Fig_ISvsRich_merge,
#        width = 12, height = 10)

# <---------------------------------------------> #
#  Save figures
# <---------------------------------------------> #
# Main figures
## Figure 1
Fig_SiteDynamics <- plot_grid(Fig_SiteMap,
                              Fig_DNApatternCombined,
                              ncol = 1, rel_heights = c(1, 2.3),
                              labels = c("a", NA))
ggsave(file = sprintf("%s/Figure_01.pdf", fig_folder2),
       plot = Fig_SiteDynamics, width = 8, height = 12)
## Figure 2
quartz(file = sprintf("%s/Figure_02.jpg", fig_folder2),
       type = "jpg", width = 22, height = 14, dpi = 300)
Fig_NetworkAll; dev.off()
## Figure 3
ggsave(file = sprintf("%s/Figure_03.pdf", fig_folder2),
       plot = Fig_ISenv_Main, width = 12, height = 8, device = cairo_pdf)
## Figure 4
ggsave(file = sprintf("%s/Figure_04.pdf", fig_folder2),
       plot = Fig_ISvsTemp_merge, width = 14, height = 14)

# Supplementary figures
## Figure S1
ggsave(file = sprintf("%s/Figure_S01.pdf", fig_folder2),
       plot = Fig_qPCR, width = 12, height = 8)
## Figure S2
ggsave(file = sprintf("%s/Figure_S02.pdf", fig_folder2),
       plot = Fig_Pattern2, width = 12, height = 10)
## Figure S3
cairo_pdf(file = sprintf("%s/Figure_S03.pdf", fig_folder2),
       width = 10, height = 10)
Fig_CorPairs; dev.off()
## Figure S4
ggsave(file = sprintf("%s/Figure_S04.pdf", fig_folder2),
       plot = Fig_ISenv_SI, width = 12, height = 8, device = cairo_pdf)
## Figure S5
ggsave(file = sprintf("%s/Figure_S05.jpg", fig_folder2),
       plot = Fig_ISoverall1, width = 10, height = 5.5, dpi = 300)
## Figure S6
ggsave(file = sprintf("%s/Figure_S06.jpg", fig_folder2),
       plot = Fig_ISoverall2, width = 10, height = 5.5, dpi = 300)
## Figure S7
ggsave(file = sprintf("%s/Figure_S07.pdf", fig_folder2),
       plot = Fig_ISvsRich_merge, width = 14, height = 11)
## Figure S8
ggsave(file = sprintf("%s/Figure_S08.pdf", fig_folder2),
       plot = Fig_ISvsDNA_merge, width = 14, height = 14)


# ------------------------------------------- #
# Appendix
# ------------------------------------------- #
