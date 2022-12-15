####
#### Boso Peninsula project
#### No.1: Compile data and correct DNA conc.
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder01 <- "01_CompileDataOut"
dir.create("00_SessionInfo")
dir.create(output_folder01)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.10.3
library(lubridate); packageVersion("lubridate") # 1.7.9, 2020.10.3
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.10.3
library(ggsci); packageVersion("ggsci") # 2.9, 2020.10.3
theme_set(theme_cowplot())

# Load original data
asv_df <- read.csv("ed_data/asv_table.csv", row.names = 1)
sample_df <- read.csv("ed_data/sample_sheet.csv", row.names = 1)
tax_df <- read.csv("ed_data/tax_sheet.csv", row.names = 1)
# Data sheet for sample check only
dnaconc_df <- read.csv("ed_data/0_for_preparation/dnaconc_sheet.csv", row.names = 1)
samplecheck <- read.csv("ed_data/0_for_preparation/sample_code_check.csv", row.names = 1)
taxacheck <- read.csv("ed_data/0_for_preparation/taxa_code_check.csv", row.names = 1)

# Convert date object
sample_df$date <- ymd(sample_df$date)

# Sample check
check01 <- all(rownames(samplecheck) == rownames(sample_df))
check02 <- all(paste(samplecheck$code2, samplecheck$code3) == paste(sample_df$time_code, sample_df$site_code))
check03 <- all(as.character(sample_df$site_name) == as.character(samplecheck$code4))
#samplecheck$code5
check04 <- all(rownames(taxacheck) == rownames(tax_df))
check05 <- all(taxacheck$Family == tax_df$family_w_n)
check06 <- all(taxacheck$Common.Name == tax_df$common_jp_name)
check07 <- all(taxacheck$Coastal.fish == tax_df$coastal_fish)

# Check row.names and col.names
## Sample ID
check08 <- all(rownames(asv_df) == rownames(sample_df))
check09 <- all(rownames(asv_df) == rownames(dnaconc_df))
check10 <- all(paste(sample_df$time_code, sample_df$site_code) == paste(dnaconc_df$time_code, dnaconc_df$site_code))
## ASV ID
check11 <- all(rownames(tax_df) == colnames(asv_df))

# Total check
all(check01, check02, check03, check04, check05, check06,
    check07, check08, check09, check10, check11)
if(all(check01, check02, check03, check04, check05, check06,
    check07, check08, check09, check10, check11)){
  rm(check01, check02, check03, check04, check05, check06,
     check07, check08, check09, check10, check11)
}
## All TRUE! => Data check OK


# Replace 0 A.schlegelii data with the minimum values
sample_df$qpcr_2ul_replaced <- sample_df$qpcr_2ul
min_qpcr <- min(sample_df$qpcr_2ul[sample_df$qpcr_2ul > 0])
sample_df$qpcr_2ul_replaced[sample_df$qpcr_2ul < min_qpcr] <- min_qpcr
## Zero DNA copies detected from 16.7% of total samples (92/550). (83.3% positive)

# Assign A.schlegelii read count to sample df
a_schlegelii_reads <- as.numeric(unlist(asv_df[tax_df$scientific_name == "Acanthopagrus schlegelii"]))
sample_df$stdfish_reads_replaced <- sample_df$stdfish_reads <- a_schlegelii_reads
## Replace 0 A.schlegelii reads with minimum reads
min_reads <- min(a_schlegelii_reads[a_schlegelii_reads > 0])
sample_df$stdfish_reads_replaced[sample_df$stdfish_reads < min_reads] <- min_reads
## Zero reads detected from 8.4% of total samples (46/550). (91.6% positive)

# Add stdfish positive or negative
sample_df$stdfish_positive <- sample_df$stdfish_reads > 0 & sample_df$qpcr_2ul > 0

# Calculate conversion factor
sample_df$conv_factor <- sample_df$qpcr_2ul_replaced / sample_df$stdfish_reads_replaced
hist(log(sample_df$conv_factor)); which.min(sample_df$conv_factor)
sample_df$total_reads <- rowSums(asv_df)
sample_df$total_dna_estimated <- rowSums(asv_df) * sample_df$conv_factor
sample_df$total_dna_per_ml <- sample_df$total_dna_estimated * (100/2) / sample_df$filtered_vol

# Visualization
## Visualize A.schlegelii DNA
g1 <- ggplot(sample_df, aes(x = date, y = log(qpcr_2ul + 0.5), color = site_code, group = site_code)) +
  geom_point() + geom_line() + facet_wrap(.~ site_code) + 
  geom_hline(yintercept = log(0.6913333 + 0.5), linetype = 1, color = "red3") +
  geom_hline(yintercept = log(0.5), linetype = 2) +
  scale_color_igv() + xlab(NULL) + ylab("log(DNA conc. + 0.5) copies / 2 µl") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12))

pdf(sprintf("%s/Aschlegelii_DNA.pdf", output_folder01), width = 14, height = 6)
g1; dev.off()

g2 <- ggplot(NULL, aes(y = a_schlegelii_reads, x = sample_df$qpcr_2ul)) +
  geom_point(alpha = 0.3) + xlab("qPCR (copies/2 µl)") + ylab("A.schlegelii reads") +
  geom_smooth(method = "lm")
pdf(sprintf("%s/qPCR_vs_Reads.pdf", output_folder01), width = 5.5, height = 5)
g2; dev.off()

g3 <- ggplot(sample_df, aes(x = date, y = log(total_dna_per_ml + 0.5), color = site_code, group = site_code)) +
  geom_point() + geom_line() + facet_wrap(.~ site_code) + 
  scale_color_igv() + xlab(NULL) + ylab("log(Total DNA copies + 0.5) / ml water") +
  geom_hline(yintercept = log(0.5), linetype = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12))

pdf(sprintf("%s/Total_DNA.pdf", output_folder01), width = 14, height = 6)
g3; dev.off()


# Save converted DNA table
asv_df_conv <- asv_df * sample_df$conv_factor * (100/2) / sample_df$filtered_vol

# Output all data
write.csv(asv_df, "data/asv_table.csv")
write.csv(asv_df_conv, "data/asv_table_conv.csv")
write.csv(sample_df, "data/sample_sheet.csv")
write.csv(tax_df, "data/tax_sheet.csv")

# Save and output results
save.image(sprintf("%s/%s.RData", output_folder01, output_folder01))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_%s.txt",  output_folder01, substr(Sys.time(), 1, 10)))

