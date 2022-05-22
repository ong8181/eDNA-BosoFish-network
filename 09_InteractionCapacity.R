####
#### Boso Peninsula project
#### No. 9 Test interaction capacity hypothesis
####

# Load workspace
load("07_EnvironmentNetworkOut/07_EnvironmentNetworkOut.RData")

# Load tidyverse
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.8.25
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.8.25
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.3.1
library(plotly); packageVersion("plotly") # 4.10.0, 2021.11.9
library(mgcv); packageVersion("mgcv") # 1.8.38, 2021.11.9
library(ggsci); packageVersion("ggsci") # 2.9, 2021.8.26
options(tibble.print_min = 20)
options(tibble.width = Inf)
theme_set(theme_cowplot())

# Set random seeds (for reproduction)
od_name <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od_name, end = -3), "Out")); rm(od_name)
dir.create(output_folder)


# <------------------------------------------------------> #
#  Calculate site-specific interaction strength (= TE)
# <------------------------------------------------------> #
site_te <- uic_same_st %>% group_by(site) %>%
  summarize(te_mean = mean(te, na.rm = T),
            te_med = median(te, na.rm = T),
            te_max = max(te, na.rm = T),
            te_min = min(te, na.rm = T)) %>%
  bind_cols(site_summary_df)


# <---------------------------------------------> #
#  Visualization
# <---------------------------------------------> #
s1 <- ggplot(site_te, aes(x = richness, y = te_mean)) +
  stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = F, size = 0.2) +
  geom_point(alpha = 1) +
  xlab("Community diversity") +
  ylab("Mean TE") +
  NULL
summary(lm(te_mean ~ temp_mean + richness + total_edna, data = site_te))
rmd1 <- lm(richness ~ temp_mean + total_edna + te_mean, data = site_te)
rmd2 <- lm(richness ~ temp_mean + total_edna, data = site_te)
rmd3 <- lm(richness ~ temp_mean + te_mean, data = site_te)
rmd4 <- lm(richness ~ total_edna + te_mean, data = site_te)
rmd5 <- lm(richness ~ temp_mean, data = site_te)
rmd6 <- lm(richness ~ te_mean, data = site_te)
rmd7 <- lm(richness ~ total_edna, data = site_te)
AIC(rmd1, rmd2, rmd3, rmd4, rmd5, rmd6, rmd7)
BIC(rmd1, rmd2, rmd3, rmd4, rmd5, rmd6, rmd7)


## 3D scatter plot
library(akima)
library(rgl)
model1 <- gam(te_mean ~ te(temp_mean) + te(richness), data = site_te)
x_seq <- seq(min(site_te$temp_mean, na.rm=TRUE), max(site_te$temp_mean, na.rm=TRUE), length=50)
y_seq <- seq(min(site_te$richness, na.rm=TRUE), max(site_te$richness, na.rm=TRUE), length=50)
predfun <- function(x, y){
  newdat <- data.frame(temp_mean = x, richness = y)
  predict(model1, newdata=newdat)
}
fit <- outer(x_seq, y_seq, Vectorize(predfun))

plot_ly() %>% 
  add_markers(x = ~site_te$temp_mean, y = ~site_te$richness, z = ~site_te$te_mean) %>% 
  add_surface(x = ~x_seq, y = ~y_seq, z = t(fit), alpha = 0.5) %>%
  layout(scene = list(xaxis = list(title = "Temperature"),
                      yaxis = list(title = "Richness"),
                      zaxis = list(title = "Mean TE")))


# <---------------------------------------------> #
# Save results
# <---------------------------------------------> #
# Save workspace and objects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder, substr(Sys.time(), 1, 10)))

