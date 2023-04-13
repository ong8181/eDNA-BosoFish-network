####
#### Boso Peninsula project
#### No.1: Visualize eDNA dynamics
####

#setwd("FigCode/")

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2022.11.10
library(lubridate); packageVersion("lubridate") # 1.8.0, 2021.11.10
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.10.13
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.3, 2022.11.10
library(ggsci); packageVersion("ggsci") # 2.9, 2020.10.3
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"
dir.create(fig_folder1)
dir.create(fig_folder2)

# Load original data
load("../02_VisualizePatternsOut/02_VisualizePatternsOut.RData")

# Prepare color palette
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
palette_custom <- get_palette(11)[c(1:4,7:11,6:5)]


#------------------------------------------#
# Visualize each fish dynamics
#------------------------------------------#
## Add temperature information
site_med_temp <- sample_df %>% group_by(site_code) %>%
  summarize(temp_med = median(water_temp, na.rm = T)) %>%
  pull(temp_med)
sample_df$water_temp_med <- rep(site_med_temp, 50)

dir.create(sprintf("%s/FigDNAconc", fig_folder1))
dir.create(sprintf("%s/FigSeqReads", fig_folder1))
if(F){
  for(taxa_n in 1:50){ # Top 50 fish species
    taxa_id <- taxa_orders_freq[taxa_n]
    g_taxa <- ggplot(sample_df, aes(x = date, y = log(asv_df_conv[,taxa_id]+0.5), color = site_code, group = site_code)) +
      geom_point() + geom_line() + facet_wrap(.~ site_code) + 
      geom_hline(yintercept = log(0.5), linetype = 2) +
      ggtitle(bquote(italic(.(tax_df[taxa_id, "scientific_name"])))) + 
      #ggtitle(tax_df[taxa_id, "scientific_name"]) + 
      scale_color_manual(values = palette_custom, name = "Site code") +
      xlab(NULL) + ylab("log(DNA conc. + 0.5) copies / ml water") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      NULL
    g_reads <- ggplot(sample_df, aes(x = date, y = log(asv_df[,taxa_id]+0.5), color = site_code, group = site_code)) +
      geom_point() + geom_line() + facet_wrap(.~ site_code) + 
      geom_hline(yintercept = log(0.5), linetype = 2) +
      ggtitle(bquote(italic(.(tax_df[taxa_id, "scientific_name"])))) + 
      scale_color_manual(values = palette_custom, name = "Site code") +
      xlab(NULL) + ylab("log(Sequence reads + 0.5)") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      NULL
    # Save figures
    quartz(file = sprintf("%s/FigDNAfreq/Top%03d_%s_conc.png",
                          fig_folder1, taxa_n, tax_df[taxa_id, "scientific_name"]),
           type = "png", width = 14, height = 6, dpi = 200)
    print(g_taxa); dev.off()
    quartz(file = sprintf("%s/FigSeqReads/Top%03d_%s_reads.png",
                          fig_folder1, taxa_n, tax_df[taxa_id, "scientific_name"]),
           type = "png", width = 14, height = 6, dpi = 200)
    print(g_reads); dev.off()
  }
}


###################################################################################

#------------------------------------------#
# Visualize key property dynamics
#------------------------------------------#
# 1. Visualize A.schlegelii DNA
#(p1 <- ggplot(sample_df, aes(x = date, y = qpcr_2ul * 50 / filtered_vol + 0.5, color = site_code, group = site_code)) +
    #geom_hline(yintercept = 0.6913333 * 50 / filtered_vol + 0.5, linetype = 1, color = "gray60") +
(p1 <- ggplot(sample_df, aes(x = date, y = qpcr_2ul + 0.5, color = site_code, group = site_code)) +
    geom_hline(yintercept = 0.6913333 + 0.5, linetype = 1, color = "gray60") +
    geom_point() + geom_line() + facet_wrap(.~ site_code) + 
    geom_hline(yintercept = 0.5, linetype = 2) +
    panel_border() +
    scale_color_manual(values = palette_custom, name = "Site code") +
    xlab(NULL) + ylab(expression(paste("DNA (+ 0.5) (copies / 2 ", mu, "l DNA extract)"))) +
    #xlab(NULL) + ylab("DNA (+ 0.5) copies / ml water") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12)) +
    scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
    ggtitle(expression(paste(italic("Acanthopagrus schlegelii"), " eDNA quantified by qPCR"))) +
    NULL)
# 2. qPCR vs Sequence reads
(p2 <- ggplot(NULL, aes(y = a_schlegelii_reads, x = sample_df$qpcr_2ul)) +
    geom_smooth(method = "lm") +
    geom_point(alpha = 0.3) + xlab("qPCR (copies/2 µl)") +
    ylab(expression(paste(italic("A.schlegelii"), " reads"))) +
    NULL)
# 3. Total DNA copies
(p3 <- ggplot(sample_df, aes(x = date, y = total_dna_per_ml + 0.5, color = site_code, group = site_code)) +
    geom_point() + geom_line() + facet_wrap(.~ site_code) + 
    scale_color_manual(values = palette_custom)+ xlab(NULL) + ylab("Total DNA copies (+ 0.5) / ml water") +
    panel_border() +
    geom_hline(yintercept = 0.5, linetype = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12)) +
    scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
    ggtitle("Total eDNA concentrations") +
    NULL)
(p3_2 <- ggplot(sample_df, aes(x = date, y = total_dna_per_ml + 0.5, color = site_code, group = site_code)) +
    geom_point(alpha = 0.5) + #geom_line() + 
    stat_smooth(se = FALSE, alpha = 0.5) +
    scale_color_manual(name = "Site ID", values = palette_custom) +
    xlab(NULL) +
    ylab("Total DNA copies (+ 0.5) / ml water") +
    geom_hline(yintercept = log(0.5), linetype = 2) +
    scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12)) +
    NULL)

all(names(rowSums(asv_df_conv > 0)) == rownames(sample_df))
sample_df$total_div <- rowSums(asv_df_conv > 0)
# 4. Total diversity
(p4 <- ggplot(sample_df, aes(x = date, y = total_div, color = site_code, group = site_code)) +
    geom_point() + geom_line() + facet_wrap(.~ site_code) + 
    scale_color_manual(values = palette_custom) + xlab(NULL) + ylab("Species richness") +
    geom_hline(yintercept = log(0.5), linetype = 2) +
    panel_border() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12)) +
    ggtitle("Species richness") +
    NULL)
(p4_2 <- ggplot(sample_df, aes(x = date, y = total_div, color = site_code, group = site_code)) +
    geom_point(alpha = 0.5) + #geom_line() + 
    stat_smooth(se = FALSE, alpha = 0.5) +
    #scale_color_gradient2(low = "royalblue", high = "red3", midpoint = 19) +
    scale_color_manual(name = "Site ID",
                       values = palette_custom)+ xlab(NULL) + ylab("Species richness") +
    geom_hline(yintercept = log(0.5), linetype = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12)) +
    NULL)



#------------------------------------------#
# Supporting figures
#------------------------------------------#
# Scattered plots
## Quant-IT (log)
dna_lm <- sample_df %>%
  smatr::sma("total_dna_per_ml + 0.5 ~ water_temp", data = ., log = "y")
dna_lm_slope <- dna_lm$coef[[1]][2,1]
dna_lm_intct <- dna_lm$coef[[1]][1,1]
div_lm <- sample_df %>%
  smatr::sma("total_div ~ water_temp", data = .)
div_lm_slope <- div_lm$coef[[1]][2,1]
div_lm_intct <- div_lm$coef[[1]][1,1]

(p5_1 <- ggplot(sample_df, aes(x = water_temp, y =  total_dna_per_ml + 0.5, color = site_code, group = site_code)) +
    geom_abline(intercept = dna_lm_intct, slope = dna_lm_slope, linetype = 2) +
    geom_point() +
    scale_color_manual(name = "Site ID", values = palette_custom) +
    scale_y_continuous(trans = "log10", labels = macam::label_10_to_power) +
    xlab(expression(paste("Water temperature (", degree, "C)"))) +
    ylab("Total DNA copies (+ 0.5) / ml water") +
    NULL)
(p5_2 <- ggplot(sample_df, aes(x = water_temp, y =  total_div, color = site_code, group = site_code)) +
    geom_abline(intercept = div_lm_intct, slope = div_lm_slope, linetype = 2) +
    geom_point() +
    scale_color_manual(name = "Site ID", values = palette_custom) +
    xlab(expression(paste("Water temperature (", degree, "C)"))) +
    ylab("Species richness") +
    NULL)
(p5_3 <- ggplot(sample_df, aes(x = date, y =  water_temp, color = site_code, group = site_code)) +
    geom_line() +
    scale_color_manual(name = "Site ID", values = palette_custom) +
    xlab(NULL) +
    ylab(expression(paste("Water temperature (", degree, "C)"))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12)) +
    NULL)
(p6 <- ggplot(data = NULL, aes(y = rowSums(asv_df_conv), x = ddf$conc38_interpolate)) +
    geom_point() + geom_smooth(method = "lm") + 
    scale_y_log10(limits = c(0.05, 1000)) + scale_x_log10(limits = c(0.05, 200)) +
    ylab("DNA conc corrected by qPCR") + xlab("DNA conc after 1st PCR (38-cycle)"))
summary(lm(rowSums(asv_df_conv) ~ ddf$conc38_interpolate))
summary(lm(ddf$conc38_interpolate ~ rowSums(asv_df_conv)))


#------------------------------------------#
# Save figures and objects
#------------------------------------------#
# Save PDF
ggsave(filename = sprintf("%s/PDF_Aschlegelii_DNA.pdf", fig_folder1),
       plot = p1, width = 14, height = 6)
ggsave(filename = sprintf("%s/PDF_qPCR_vs_Reads.pdf", fig_folder1),
       plot = p2, width = 5.5, height = 5)
ggsave(filename = sprintf("%s/PDF_Total_DNA.pdf", fig_folder1),
       plot = p3, width = 14, height = 6)
ggsave(filename = sprintf("%s/PDF_Total_Div.pdf", fig_folder1),
       plot = p4, width = 14, height = 6)
ggsave(filename = sprintf("%s/PDF_DNAconc38_vs_TapaStation.pdf", fig_folder1),
       plot = p6, width = 6, height = 6)
ggsave(filename = sprintf("%s/PDF_Total_DNA2.pdf", fig_folder1),
       plot = p3_2, width = 6, height = 4)
ggsave(filename = sprintf("%s/PDF_Total_Div2.pdf", fig_folder1),
       plot = p4_2, width = 6, height = 4)


# Save R Objects
p_all <- list(p1, p2, p3, p3_2, p4, p4_2)
p2_all <- list(p5_1, p5_2, p5_3, p6)
saveRDS(p_all, sprintf("%s/Fig_DNApattern01.obj", fig_folder1))
saveRDS(p2_all, sprintf("%s/Fig_DNApattern02.obj", fig_folder1))


# Pre-assembled for format
Fig_Pattern0 <- plot_grid(p3_2 +
                            theme(axis.text.x = element_text(angle = 0)) +
                            theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 0),
                                  panel.grid.major.x = element_line(linetype = 2, color = "gray70"),
                                  panel.grid.minor.x = element_line(linetype = 2, color = "gray70")),
                          p4_2 + theme(axis.text.x = element_text(angle = 0))  +
                            theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 0),
                                  panel.grid.major.x = element_line(linetype = 2, color = "gray70"),
                                  panel.grid.minor.x = element_line(linetype = 2, color = "gray70")),
                          ncol = 1, axis = "lrbt", align = "hv", labels = c("b","c"))
site_legend <- get_legend(p3_2)
Fig_Pattern <- plot_grid(Fig_Pattern0, site_legend, rel_widths = c(1,0.1))
saveRDS(Fig_Pattern, sprintf("%s/Fig_DNApatternCombined.obj", fig_folder1))


