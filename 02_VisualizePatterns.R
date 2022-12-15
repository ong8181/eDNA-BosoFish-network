####
#### Boso Peninsula project
#### No.2: Visualize fish-specific patterns
####

# Set random seeds (for reproduction)
output_folder02 <- "02_VisualizePatternsOut"
dir.create(output_folder02)
dir.create(sprintf("%s/FigDNAconc",output_folder02))
dir.create(sprintf("%s/FigDNAfreq",output_folder02))
dir.create(sprintf("%s/FigSeqReadsConc",output_folder02))
dir.create(sprintf("%s/FigSeqReadsFreq",output_folder02))

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.10.3
library(lubridate); packageVersion("lubridate") # 1.7.9, 2020.10.3
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.10.3
library(ggsci); packageVersion("ggsci") # 2.9, 2020.10.3
theme_set(theme_cowplot())

# Load workspace
load("01_CompileDataOut/01_CompileDataOut.RData")

# Check asv_df_onv
head(asv_df_conv)
all(rownames(asv_df_conv) == rownames(sample_df))

# Visualization
## Extract dominant 20 fish spee
taxa_orders_conc <- names(sort(colSums(asv_df_conv), decreasing = T))
taxa_orders_freq <- names(sort(colSums(asv_df_conv > 0), decreasing = T))

for(taxa_n in 1:50){
  taxa_id <- taxa_orders_freq[taxa_n]
  g_taxa <- ggplot(sample_df, aes(x = date, y = log(asv_df_conv[,taxa_id]+0.5), color = site_code, group = site_code)) +
    geom_point() + geom_line() + facet_wrap(.~ site_code) + 
    geom_hline(yintercept = log(0.5), linetype = 2) +
    ggtitle(tax_df[taxa_id, "scientific_name"]) + 
    scale_color_igv() + xlab(NULL) + ylab("log(DNA conc. + 0.5) copies / ml water") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  g_reads <- ggplot(sample_df, aes(x = date, y = log(asv_df[,taxa_id]+0.5), color = site_code, group = site_code)) +
    geom_point() + geom_line() + facet_wrap(.~ site_code) + 
    geom_hline(yintercept = log(0.5), linetype = 2) +
    ggtitle(tax_df[taxa_id, "scientific_name"]) + 
    scale_color_igv() + xlab(NULL) + ylab("log(Sequence reads + 0.5)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  # Save figures
  quartz(file = sprintf("%s/FigDNAfreq/Top%03d_%s_conc.png",
                        output_folder02,
                        taxa_n,
                        tax_df[taxa_id, "scientific_name"]),
         type = "png", width = 14, height = 6, dpi = 200)
  print(g_taxa); dev.off()
  
  quartz(file = sprintf("%s/FigSeqReads/Top%03d_%s_reads.png",
                        output_folder02,
                        taxa_n,
                        tax_df[taxa_id, "scientific_name"]),
         type = "png", width = 14, height = 6, dpi = 200)
  print(g_reads); dev.off()
}

# Check correspondence between DNA conc after 1st PCR and converted DNA concs
ddf <- dnaconc_df # shorten the object name
ddf$conc35_vol <- ddf$conc_pgul_35/ddf$water_vol_ml # Corrected DNA conc. by water volume
ddf$conc38_vol <- ddf$conc_pgul_38/ddf$water_vol_ml # Corrected DNA conc. by water volume
correct_valid_id <- !is.na(ddf$conc35_vol) & !is.na(ddf$conc38_vol)
cor_lm <- lm(ddf$conc38_vol[correct_valid_id] ~ ddf$conc35_vol[correct_valid_id])
summary(cor_lm)
g4 <- ggplot(ddf[correct_valid_id,], aes(x = log(conc35_vol), y = log(conc38_vol))) +
  geom_point() + geom_smooth(method = "lm") +
  xlab("Log(DNA conc after 35 cycle)") +
  ylab("Log(DNA conc after 38 cycle") + geom_abline(intercept = 0, slope = 1) +
  xlim(-6,6) + ylim(-2,4)

# Convert 35 cycle DNA conc to 38 cycles
ddf$conc35_interpolate <- ddf$conc35_vol
ddf$conc38_interpolate <- ddf$conc38_vol
ddf$conc35_interpolate[is.na(ddf$conc35_vol)] <- (ddf$conc38_vol[is.na(ddf$conc35_vol)] - cor_lm$coefficients[1])/ cor_lm$coefficients[2]
ddf$conc38_interpolate[is.na(ddf$conc38_vol)] <- ddf$conc35_vol[is.na(ddf$conc38_vol)] * cor_lm$coefficients[2] + cor_lm$coefficients[1]
g5 <- ggplot(ddf, aes(x = log(conc35_interpolate), y = log(conc38_interpolate))) +
  geom_point() + geom_smooth(method = "lm") +
  xlab("Log(DNA conc after 35 cycle)") +
  ylab("Log(DNA conc after 38 cycle") + geom_abline(intercept = 0, slope = 1) +
  xlim(-6,6) + ylim(-2,4)

# Compare qPCR-based correction and DNA conc after 1st PCR
g6 <- ggplot(data = NULL, aes(x = rowSums(asv_df_conv), y = ddf$conc35_interpolate)) +
  geom_point() + geom_smooth(method = "lm") + scale_x_log10() + scale_y_log10() +
  xlab("DNA conc corrected by qPCR") + ylab("DNA conc after 1st PCR (35-cycle)")
g7 <- ggplot(data = NULL, aes(y = rowSums(asv_df_conv), x = ddf$conc38_interpolate)) +
  geom_point() + geom_smooth(method = "lm") + 
  scale_y_log10(limits = c(0.05, 1000)) + scale_x_log10(limits = c(0.05, 200)) +
  ylab("DNA conc corrected by qPCR") + xlab("DNA conc after 1st PCR (38-cycle)")
summary(lm(rowSums(asv_df_conv) ~ ddf$conc38_interpolate))
summary(lm(ddf$conc38_interpolate ~ rowSums(asv_df_conv)))

# Save and output results
save.image(sprintf("%s/%s.RData", output_folder02, output_folder02))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_%s.txt",  output_folder02, substr(Sys.time(), 1, 10)))

