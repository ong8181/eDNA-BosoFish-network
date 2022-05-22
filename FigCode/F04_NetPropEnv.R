####
#### Boso Peninsula project
#### No. 4 Figure: Environmental variables and network patterns
####

# Load workspace
load("../07_EnvironmentNetworkOut/07_EnvironmentNetworkOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.8.25
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(mgcv); packageVersion("mgcv") # 1.8.38, 2021.11.10
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

# Prepare output folders
fig_folder1 <- "0_RawFigs"
fig_folder2 <- "0_FormattedFigs"

# Prepare color palette
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
get_palette2 <- colorRampPalette(brewer.pal(8, "Set1"))
palette_custom <- get_palette(11)[c(1:4,7:11,6:5)]
palette_custom3 <- get_palette2(28)


# <---------------------------------------------> #
#  Visualization
# <---------------------------------------------> #
## Prepare facet labels
new_facet_lab <- c("Latitude", "Linear distance", "Richness",
                   "Richness (dominant)", "Mean salinity", "Max. temperature",
                   "Mean temperature", "Median temperature", "Min. temperature",
                   "Mean tide", "Total eDNA", "Mean wave (m)")
names(new_facet_lab) <- sort(unique(effect_te_long$name)[1:12])
## One point == one species
## (= one species has several interactions (TE values), and mean value of TE is used)
# Effect TE patterns
(e3 <- ggplot(effect_te_long %>% filter(name != "n_int"),
              aes(x = value, y = te_med, color = name)) +
    geom_point(alpha = 0.5) +
    stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
    facet_wrap(.~ name, ncol = 3, scales = "free",
               labeller = labeller(name = new_facet_lab)) +
    scale_color_manual(values = get_palette(12),
                       labels = new_facet_lab) +
    scale_y_log10() +
    xlab("Environmental variables") + ylab("Effect TE") +
    panel_border() + ggtitle("Median TE") +
    NULL)
# Causal TE patterns
(c3 <- ggplot(cause_te_long %>% filter(name != "n_int") %>% droplevels(),
              aes(x = value, y = te_med, color = name)) +
    geom_point(alpha = 0.5) +
    stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
    facet_wrap(.~ name, ncol = 3, scales = "free",
               labeller = labeller(name = new_facet_lab)) +
    scale_color_manual(values = get_palette(12),
                       labels = new_facet_lab) +
    scale_y_log10() +
    xlab("Environmental variables") + ylab("Causal TE") +
    panel_border() + ggtitle("Median TE") +
    NULL)
## Total TE patterns
(t3 <- ggplot(sum_te_long %>% filter(name != "n_int" & name != "n_int_cause" & name != "n_int_effect"),
              aes(x = value, y = te_med, color = name)) +
    geom_point(alpha = 0.5) +
    stat_smooth(color = "gray50", method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
    facet_wrap(.~ name, ncol = 3, scales = "free",
               labeller = labeller(name = new_facet_lab)) +
    scale_color_manual(values = get_palette(12),
                       label = new_facet_lab) +
    scale_y_log10() +
    xlab("Environmental variables") + ylab("Total TE") +
    panel_border() + ggtitle("Median TE") +
    NULL)


## Separate panels for the manuscript
#name_cond_main <- c("lat", "temp_med", "richness_dom", "total_edna")
#name_cond_sub <- c("temp_max", "temp_min", "salinity_mean", "richness")
# Add site ID to sum_te_long
sum_te_long$site_id <- sum_te_long$site_sp %>%
  str_sub(start = 1, end = 4) %>% as.factor()
sum_te_long$sp_id <- sum_te_long$site_sp %>%
  str_sub(start = 6) %>% as.factor()
## Latitude ------------------------------
gamm1 <- gamm(te_med ~ s(value),
             family = Gamma(link = "log"), random = list(sp_id = ~1),
             data = sum_te_long %>% filter(name == "lat"))
summary(gamm1$gam) #; summary(gamm1$lme)
summary(gamm1$gam)$s.table
x_min1 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "lat") %>% select(value))))
x_max1 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "lat") %>% select(value))))
new_d1 <- data.frame(value = seq(x_min1, x_max1, by = 0.001))
new_d1$pred <- exp(predict.gam(gamm1$gam, newdata = new_d1))
(s1 <- sum_te_long %>% filter(name == "lat") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_line(data = new_d1, aes(x = value, y = pred), color = "gray50", size = 1) +
    geom_point(alpha = 0.5, color = "gray20") +
    scale_y_log10() + ylab("Total TE") + xlab(expression(paste("Latitude (", degree, "N)"))) + 
    NULL) # Significant, P = 2.680223e-06


## Median water temperature ------------------------------
gamm2 <- gamm(te_med ~ s(value),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "temp_med"))
summary(gamm2$gam) #; summary(gamm1$lme)
summary(gamm2$gam)$s.table
x_min2 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "temp_med") %>% select(value))))
x_max2 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "temp_med") %>% select(value))))
new_d2 <- data.frame(value = seq(x_min2, x_max2, by = 0.01))
new_d2$pred <- exp(predict.gam(gamm2$gam, newdata = new_d2))
(s2 <- sum_te_long %>% filter(name == "temp_med") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_line(data = new_d2, aes(x = value, y = pred), color = "gray50", size = 1) +
    geom_point(alpha = 0.5, color = "red2") +
    scale_y_log10() + ylab("Total TE") +
    xlab(expression(atop(paste("Water temperature (", degree, "C)"), "(Median value)"))) + 
    NULL) # Significant, P = <2e-16


## Richness ------------------------------
gamm3 <- gamm(te_med ~ s(value, k = 9),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "richness"))
summary(gamm3$gam) #; summary(gamm1$lme)
summary(gamm3$gam)$s.table
x_min3 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "richness") %>% select(value))))
x_max3 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "richness") %>% select(value))))
new_d3 <- data.frame(value = seq(x_min3, x_max3, by = 0.01))
new_d3$pred <- exp(predict.gam(gamm3$gam, newdata = new_d3))
(s3 <- sum_te_long %>% filter(name == "richness") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_line(data = new_d3, aes(x = value, y = pred), color = "gray50", size = 1) +
    geom_point(alpha = 0.5, color = "deepskyblue3") +
    scale_y_log10() + ylab("Total TE") + xlab("Species richness") + 
    NULL) # Significant, P = 0.00304012


## Total eDNA ------------------------------
gamm4 <- gamm(te_med ~ s(value),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "total_edna"))
summary(gamm4$gam) #; summary(gamm1$lme)
summary(gamm4$gam)$s.table
(s4 <- sum_te_long %>% filter(name == "total_edna") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_point(alpha = 0.5, color = "royalblue") +
    scale_y_log10() + ylab("Total TE") + xlab("Total eDNA (copies/ml water)") + 
    NULL) # Marginally significant, P = 0.05777008


## Minimum water temperature ------------------------------
gamm5 <- gamm(te_med ~ s(value, k = 9),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "temp_min"))
summary(gamm5$gam) #; summary(gamm1$lme)
summary(gamm5$gam)$s.table
x_min5 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "temp_min") %>% select(value))))
x_max5 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "temp_min") %>% select(value))))
new_d5 <- data.frame(value = seq(x_min5, x_max5, by = 0.01))
new_d5$pred <- exp(predict.gam(gamm5$gam, newdata = new_d5))
(s5 <- sum_te_long %>% filter(name == "temp_min") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_point(alpha = 0.5, color = "skyblue") +
    geom_line(data = new_d5, aes(x = value, y = pred), color = "gray50", size = 1) +
    scale_y_log10() + ylab("Total TE") +
    xlab(expression(atop(paste("Water temperature (", degree, "C)"), "(Minimum value)"))) + 
    NULL) # Significant, P = 0


## Maximum water temperature ------------------------------
gamm6 <- gamm(te_med ~ s(value, k = 8),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "temp_max"))
summary(gamm6$gam) #; summary(gamm1$lme)
summary(gamm6$gam)$s.table
x_min6 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "temp_max") %>% select(value))))
x_max6 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "temp_max") %>% select(value))))
new_d6 <- data.frame(value = seq(x_min6, x_max6, by = 0.01))
new_d6$pred <- exp(predict.gam(gamm6$gam, newdata = new_d6))
(s6 <- sum_te_long %>% filter(name == "temp_max") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_point(alpha = 0.5, color = "red4") +
    geom_line(data = new_d6, aes(x = value, y = pred), color = "gray50", size = 1) +
    scale_y_log10() + ylab("Total TE") +
    xlab(expression(atop(paste("Water temperature (", degree, "C)"), "(Maximum value)"))) + 
    NULL) # Significant, P = 0.002385199


## Maximum water temperature ------------------------------
gamm7 <- gamm(te_med ~ s(value),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "salinity_mean"))
summary(gamm7$gam) #; summary(gamm1$lme)
summary(gamm7$gam)$s.table
x_min7 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "salinity_mean") %>% select(value))))
x_max7 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "salinity_mean") %>% select(value))))
new_d7 <- data.frame(value = seq(x_min7, x_max7, by = 0.01))
new_d7$pred <- exp(predict.gam(gamm7$gam, newdata = new_d7))
(s7 <- sum_te_long %>% filter(name == "salinity_mean") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_line(data = new_d7, aes(x = value, y = pred), color = "gray50", size = 1) +
    geom_point(alpha = 0.5, color = "royalblue") +
    scale_y_log10() + ylab("Total TE") + xlab(expression("Mean salinity (\u2030)")) + 
    NULL) # Significant, P = 0.002677969


## Richness (dom) temperature ------------------------------
gamm8 <- gamm(te_med ~ s(value, k = 9),
              family = Gamma(link = "log"), random = list(sp_id = ~1),
              data = sum_te_long %>% filter(name == "richness_dom"))
summary(gamm8$gam) #; summary(gamm1$lme)
summary(gamm8$gam)$s.table
x_min8 <- min(as.numeric(unlist(sum_te_long %>% filter(name == "richness_dom") %>% select(value))))
x_max8 <- max(as.numeric(unlist(sum_te_long %>% filter(name == "richness_dom") %>% select(value))))
new_d8 <- data.frame(value = seq(x_min8, x_max8, by = 0.01))
new_d8$pred <- exp(predict.gam(gamm8$gam, newdata = new_d8))
(s8 <- sum_te_long %>% filter(name == "richness_dom") %>%
    ggplot(aes(x = value, y = te_med)) +
    geom_line(data = new_d8, aes(x = value, y = pred), color = "gray50", size = 1) +
    geom_point(alpha = 0.5, color = "deepskyblue3") +
    scale_y_log10() + ylab("Total TE") + xlab("Species richness\n(rare species excluded)") + 
    NULL) # Significant, P = 0


# <---------------------------------------------> #
#  The number of interactions and abundance
# <---------------------------------------------> #
sitesp_match_id <- match(paste0(sum_te_df$sp, "_", sum_te_df$site), colnames(asv_df_sitesp))
sum_te_df$site_abundance <- colSums(asv_df_sitesp[sitesp_match_id])
n_int_df <- sum_te_df %>%
  select(site, site_sp, site_abundance, n_int, n_int_cause, n_int_effect) %>%
  pivot_longer(cols = -c(site, site_sp, site_abundance), names_to = "name")

## Label power numbers
label_func <- function(x) {
  ifelse(x == 0, "0", 
         parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
  )
}
## Prepare facet labels
new_facet_lab2 <- c("N of interactions", "N of giving interactions", "N of received interactions")
names(new_facet_lab2) <- sort(unique(n_int_df$name))


(n1 <- ggplot(n_int_df, aes(x = site_abundance, y = value, color = site)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "lm", se = F, alpha = 0.75) +
    scale_color_manual(values = palette_custom, name = "Site ID") +
    facet_wrap(. ~ name, labeller = labeller(name = new_facet_lab2)) +
    scale_x_continuous(trans = "log10", labels = label_func) +
    scale_y_log10() +
    panel_border() +
    xlab("total eDNA abundance of each species per site (copies / ml)") +
    ylab("N of interactions") +
    NULL)


# <---------------------------------------------> #
#                    Save results                 #
# <---------------------------------------------> #
# Save figures
ggsave(sprintf("%s/PDF_MedianEffectTE.pdf", fig_folder1), plot = e3, width = 12, height = 10)
ggsave(sprintf("%s/PDF_MedianCauseTE.pdf", fig_folder1), plot = c3, width = 12, height = 10)
ggsave(sprintf("%s/PDF_MedianTotalTE.pdf", fig_folder1), plot = t3, width = 12, height = 10)
ggsave(sprintf("%s/PDF_NintTotalDNA.pdf", fig_folder1), plot = n1, width = 10, height = 6)

# Save objects
saveRDS(list(e3, c3, t3), sprintf("%s/Fig_MedianTEall.obj", fig_folder1))
saveRDS(n1, sprintf("%s/Fig_NintTotalDNA.obj", fig_folder1))
saveRDS(list(effect_te_long, cause_te_long, sum_te_long), sprintf("%s/FigData_MedianTEall.obj", fig_folder1))
saveRDS(list(s1, s2, s3, s4), sprintf("%s/FigData_MedianTE_main.obj", fig_folder1))
saveRDS(list(s5, s6, s7, s8), sprintf("%s/FigData_MedianTE_sub.obj", fig_folder1))


