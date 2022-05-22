####
#### Boso Peninsula project
#### Format figures
####

if(basename(getwd()) != "FigCode") setwd("FigCode")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(lubridate); packageVersion("lubridate") # 1.8.0, 2021.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(ggraph); packageVersion("ggraph") # 2.0.5, 2021.3.1
library(igraph); packageVersion("igraph") # 1.2.11, 2022.4.1
library(ggpubr); packageVersion("ggpubr") # 0.4.0, 2021.10.17
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
library(ggimage); packageVersion("ggimage") # 0.3.0, 2021.12.8
library(magick); packageVersion("magick") # 2.7.3, 2021.12.8
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"

# Prepare color palette
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.8.25
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
palette_custom <- get_palette(11)[c(1:4,7:11,6:5)]


# <---------------------------------------------> #
#  Load figure objects
# <---------------------------------------------> #
Fig_SiteMap <- image_read("0_Illustrations/SiteMap.png")
Fig_SiteMap <- ggdraw() + draw_image(Fig_SiteMap)
Fig_DNApattern01 <- readRDS(paste0(fig_folder1, "/Fig_DNApattern01.obj"))
Fig_DNApattern02 <- readRDS(paste0(fig_folder1, "/Fig_DNApattern02.obj"))
Fig_DNApatternCombined <- readRDS(paste0(fig_folder1, "/Fig_DNApatternCombined.obj"))
Fig_CorPairs <- readRDS(paste0(fig_folder1, "/Fig_CorPairs.obj"))
Fig_MedianTEall <- readRDS(paste0(fig_folder1, "/Fig_MedianTEall.obj"))
Fig_MedianTE_main <- readRDS(paste0(fig_folder1, "/FigData_MedianTE_main.obj"))
Fig_MedianTE_sub <- readRDS(paste0(fig_folder1, "/FigData_MedianTE_sub.obj"))
Fig_NetworkAll_raw <- readRDS(paste0(fig_folder1, "/Fig_NetworkAll.obj")) # Vector version
Fig_NetworkAll <- image_read("0_Illustrations/FishNetworkAll_pls_Fish.jpg") # Image version
Fig_NetworkAll <- ggdraw() + draw_image(Fig_NetworkAll)
Fig_NetworkProp <- readRDS(paste0(fig_folder1, "/Fig_NetworkProp.obj"))
Fig_TotalMedTE_env <- readRDS(paste0(fig_folder1, "/Fig_TotalMedTE_env.obj"))
Fig_TotalMedTE_env_GAM <- readRDS(paste0(fig_folder1, "/Fig_TotalMedTE_env_GAM.obj"))
Fig_NintTotalDNA <- readRDS(paste0(fig_folder1, "/Fig_NintTotalDNA.obj"))
FigData_MedianTEall <- readRDS(paste0(fig_folder1, "/FigData_MedianTEall.obj"))

# Load individual fish images
## Fish images will be in "ggfish" object.
source("F07_LoadFishImages.R")

# Site networks
net_n <- list(NULL)
for(site_i in 1:11){
  net_n_i <- readRDS(sprintf("%s/Fig_SiteNetworks/Fig_FishNetwork_St%02d.obj", fig_folder1, site_i)) +
    ggtitle(sprintf("Site%02d", site_i))
  net_n[[site_i]] <- net_n_i
}


# <---------------------------------------------> #
#  Compile figures
# <---------------------------------------------> #
net_n2 <- net_n
for(i in 1:length(net_n)) {
  net_n2[[i]] <- net_n2[[i]] + theme(legend.position = "none")
  net_n2[[i]]$layers[[3]] <- NULL
  net_n2[[i]]$layers[[3]] <- NULL
}

# Figure 1
Fig_SiteDynamics <- plot_grid(Fig_SiteMap,
                              Fig_DNApatternCombined,
                              ncol = 1, rel_heights = c(1, 2.5),
                              labels = c("a", NA))

# Figure 2
# Load network figures
net_legend_all <- get_legend(Fig_NetworkAll_raw + theme(legend.text = element_text(size = 10)))
net_legend1 <- ggpubr::as_ggplot(net_legend_all$grobs[[1]])
net_legend2 <- ggpubr::as_ggplot(net_legend_all$grobs[[2]])
net_legend3 <- ggpubr::as_ggplot(net_legend_all$grobs[[3]])
net_legend <- plot_grid(NULL, plot_grid(net_legend1, net_legend2, NULL, ncol = 3),
                        net_legend3, NULL, ncol = 1, rel_heights = c(1, 1.2, 3, 1))

Fig_NetworkAll2 <- plot_grid(Fig_NetworkAll, net_legend,
                             ncol = 2, rel_widths = c(2,1), labels = "a")
Fig_SiteNetworks <- plot_grid(net_n2[[1]], net_n2[[2]], net_n2[[3]], net_n2[[4]],
                              net_n2[[5]], net_n2[[6]], net_n2[[7]], net_n2[[8]],
                              net_n2[[9]], net_n2[[10]], net_n2[[11]],
                              ncol = 4, labels = letters[2:12])
Fig_NetworkAll3 <- plot_grid(Fig_NetworkAll2, Fig_SiteNetworks, ncol = 1,
                             rel_heights = c(0.9,1))

# Figure 3
Fig_MedianTE1 <- plot_grid(Fig_MedianTE_main[[1]] + ylab("Interaction strength (TE)"),
                           Fig_MedianTE_main[[2]] + ylab("Interaction strength (TE)") +
                             xlab(expression(atop(paste("Water temperature (", degree, "C)")))) ,
                           Fig_MedianTE_main[[3]] + ylab("Interaction strength (TE)"),
                           Fig_MedianTE_main[[4]] + ylab("Interaction strength (TE)"),
                           ncol = 2, labels = letters[1:4], align = "hv")
Fig_MedianTE_main[[2]] + coord_cartesian(ylim=c(0.02, 0.2))

# Figure 4
## Interaction strength + temperature
## Add fish images
Fig_TotalMedTE_sp_tmp <- Fig_TotalMedTE_env[[2]] +
  ylab("Interaction strength (TE)") +
  xlab(expression(paste("Median water temperature (", degree, "C)"))) +
  theme(legend.position = "none", strip.text = element_text(size = 8)) +
  ## Add fish images
  # 1st row
  patchwork::inset_element(ggfish$Canthigaster_rivulata,     left = 0.00, right = 0.05, bottom = 0.84, top = 0.89) +
  patchwork::inset_element(ggfish$Chaenogobius_annularis,    left = 0.30, right = 0.36, bottom = 0.84, top = 0.89) +
  patchwork::inset_element(ggfish$Dictyosoma_rubrimaculatum, left = 0.52, right = 0.57, bottom = 0.84, top = 0.89) +
  patchwork::inset_element(ggfish$Enchelycore_pardalis,      left = 0.74, right = 0.79, bottom = 0.84, top = 0.89) +
  patchwork::inset_element(ggfish$Engraulis_japonicus,       left = 0.945, right = 0.995, bottom = 0.84, top = 0.89) +
  # 2nd row
  patchwork::inset_element(ggfish$Entomacrodus_stellifer_stellifer,  left = 0.09, right = 0.15, bottom = 0.63, top = 0.68) +
  patchwork::inset_element(ggfish$Eviota_abax,                       left = 0.21, right = 0.27, bottom = 0.63, top = 0.68) +
  patchwork::inset_element(ggfish$Girella_lenonina,                  left = 0.52, right = 0.57, bottom = 0.63, top = 0.68) +
  patchwork::inset_element(ggfish$Hypoatherina_tsurugae,             left = 0.73, right = 0.78, bottom = 0.63, top = 0.68) +
  patchwork::inset_element(ggfish$Iso_flosmaris,                     left = 0.94, right = 0.99, bottom = 0.63, top = 0.68) +
  # 3rd row
  patchwork::inset_element(ggfish$Istiblennius_enosimae,   left = 0.09, right = 0.15, bottom = 0.42, top = 0.47) +
  patchwork::inset_element(ggfish$Microcanthus_strigatus,  left = 0.31, right = 0.37, bottom = 0.425, top = 0.475) +
  patchwork::inset_element(ggfish$Mugil_cephalus,          left = 0.52, right = 0.57, bottom = 0.42, top = 0.47) +
  patchwork::inset_element(ggfish$Oplegnathus_fasciatus,   left = 0.73, right = 0.789, bottom = 0.42, top = 0.47) +
  patchwork::inset_element(ggfish$Ostracion_immaclatus,    left = 0.95, right = 1.00, bottom = 0.42, top = 0.47) +
  # 4th row
  patchwork::inset_element(ggfish$Parablennius_yatabei,      left = 0.00, right = 0.05, bottom = 0.21, top = 0.26) +
  patchwork::inset_element(ggfish$Pempheris_schwenkii,       left = 0.31, right = 0.37, bottom = 0.21, top = 0.26) +
  patchwork::inset_element(ggfish$Prionurus_scalprum,        left = 0.52, right = 0.57, bottom = 0.21, top = 0.26) +
  patchwork::inset_element(ggfish$Pseudoblennius_percoides,  left = 0.73, right = 0.789, bottom = 0.21, top = 0.26) +
  patchwork::inset_element(ggfish$Pseudolabrus_eoethinus,    left = 0.95, right = 1.00, bottom = 0.21, top = 0.26) +
  # 5th row
  patchwork::inset_element(ggfish$Siganus_fuscescens,              left = 0.09, right = 0.15, bottom = 0.00, top = 0.05) +
  patchwork::inset_element(ggfish$Spratelloides_gracilis,          left = 0.31, right = 0.37, bottom = 0.00, top = 0.05) +
  patchwork::inset_element(ggfish$Stethojulis_interrupta_terina,   left = 0.52, right = 0.57, bottom = 0.00, top = 0.05) +
  patchwork::inset_element(ggfish$Thalassoma_cupido,               left = 0.64, right = 0.69, bottom = 0.00, top = 0.05) +
  patchwork::inset_element(ggfish$Trachurus_japonicus,             left = 0.95, right = 1.00, bottom = 0.00, top = 0.05) +
  NULL
## Save temporal image
quartz(file = sprintf("%s/Fig_TotalMedTE_sp_tmp.jpg", fig_folder2),
       type = "jpg", width = 12, height = 10, dpi = 600)
Fig_TotalMedTE_sp_tmp; dev.off()
Fig_TotalMedTE_sp_tmp <- image_read("0_FormattedFigs/Fig_TotalMedTE_sp_tmp.jpg") # Image version
Fig_TotalMedTE_sp_tmp <- ggdraw() + draw_image(Fig_TotalMedTE_sp_tmp)
system("rm 0_FormattedFigs/Fig_TotalMedTE_sp_tmp.jpg")

## Combine panels
Fig_TEenv <- plot_grid(Fig_TotalMedTE_sp_tmp,
                       plot_grid(NULL,
                                 Fig_TotalMedTE_env[[3]] +
                                   scale_y_continuous(breaks = seq(0,13,2)) +
                                   theme(legend.position = "right"),
                                 NULL, nrow = 1, rel_widths = c(0.02,1,0.4)),
                       ncol = 1, rel_heights = c(1,0.3),
                       labels = c("a", "b"))


# <---------------------------------------------> #
#  Supplementary figures
# <---------------------------------------------> #
label_func <- function(x) {
  ifelse(x == 0, "0", 
         parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
  )
}

# Japanese black seabream eDNA copies
Fig_qPCR <- Fig_DNApattern01[[1]] +
  scale_color_manual(name = "Site ID", values = palette_custom) +
  theme(panel.grid.major.x = element_line(linetype = 2, size = 0.4, color = "gray70")) +
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

# TE v.s. environmental variables
Fig_MedianTE2 <- plot_grid(Fig_MedianTE_sub[[1]] + ylab("Interaction strength (TE)"),
                           Fig_MedianTE_sub[[2]] + ylab("Interaction strength (TE)"),
                           Fig_MedianTE_sub[[3]] + ylab("Interaction strength (TE)"),
                           Fig_MedianTE_sub[[4]] + ylab("Interaction strength (TE)"),
                           ncol = 2, labels = letters[1:4], align = "hv")

## Combine panels
Fig_TEenv_GAM <- plot_grid(Fig_TotalMedTE_env_GAM[[1]],
                           plot_grid(NULL,
                                     Fig_TotalMedTE_env_GAM[[2]] +
                                       theme(legend.position = "right"),
                                     NULL, nrow = 1, rel_widths = c(0.02,1,0.7)),
                           ncol = 1, rel_heights = c(1,0.4),
                           labels = c("a", "b"))


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
#quartz(file = sprintf("%s/Figure_02.pdf", fig_folder2),
#       type = "pdf", width = 15, height = 22)
#Fig_NetworkAll3; dev.off()
quartz(file = sprintf("%s/Figure_02.jpg", fig_folder2),
       type = "jpg", width = 15, height = 22, dpi = 300)
Fig_NetworkAll3; dev.off()
## Figure 3
ggsave(file = sprintf("%s/Figure_03.pdf", fig_folder2),
       plot = Fig_MedianTE1, width = 10, height = 8)
## Figure 4
ggsave(file = sprintf("%s/Figure_04.pdf", fig_folder2),
       plot = Fig_TEenv, width = 12, height = 14)

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
       plot = Fig_TEenv_GAM, width = 18, height = 16)

