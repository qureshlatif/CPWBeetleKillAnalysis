library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
#setwd("/home/RMBO.LOCAL/rstudio03")
setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

samp.ha <- sum(area.band) * 0.0001

## Spruce-fir stratum ##
stratum <- "SF" #Set to LP or SF
mod <- loadObject("mod_RESQ_outbreak_HZdist_SF_reduced")

# Plot outbreak relationships #
Cov <- eval(as.name(str_c("Cov.", stratum)))
x.pctdead <- seq(0, 1, length.out = 20)
x.YSO <- seq(0, max(Cov[, "YSO"], na.rm = T))
plt.tbl <- expand.grid(x.pctdead, x.YSO) %>%
  rename(pctdead = Var1, YSO = Var2) %>%
  mutate(pctdead.z = (pctdead - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)) %>%
  mutate(YSO.z = (YSO - mean(Cov[, "YSO"], na.rm = T)) / sd(Cov[, "YSO"], na.rm = T)) %>%
  mutate(pctdead = pctdead * 100)

pred.N <- matrix(NA, nrow = nrow(plt.tbl), ncol = mod$mcmc.info$n.samples)
for(j in 1:nrow(pred.N))
  pred.N[j, ] <- (exp(mod$sims.list$beta0.mean +
                       mod$sims.list$bl.pdead * plt.tbl$pctdead.z[j] +
                       mod$sims.list$bl.YSO * plt.tbl$YSO.z[j] +
                       mod$sims.list$bl.pdXYSO * plt.tbl$pctdead.z[j] * plt.tbl$YSO.z[j]) /
    samp.ha) * 100
plt.tbl$N.md <- apply(pred.N, 1, median)
plt.tbl$N.95lo <- apply(pred.N, 1, function(x) quantile(x, prob = 0.025, type = 8))
plt.tbl$N.95hi <- apply(pred.N, 1, function(x) quantile(x, prob = 0.975, type = 8))

dat.plt <- plt.tbl %>% filter(YSO %in% c(min(YSO), max(YSO))) %>%
   mutate(YSO = as.factor(YSO))
p.pdead <- ggplot(dat.plt, aes(x = pctdead, y = N.md)) +
  geom_ribbon(aes(ymin = N.95lo, ymax = N.95hi, fill = YSO), alpha = 0.3) +
  geom_line(aes(color = YSO), size = 2) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  xlab("Dead conifer (%)") + ylab(NULL) +
  labs(fill = "Years\nsince\noutbreak") +
  labs(color = "Years\nsince\noutbreak")

dat.plt <- plt.tbl %>% filter(pctdead %in% c(min(pctdead), max(pctdead))) %>%
  mutate(pctdead = as.factor(pctdead))
p.YSO <- ggplot(dat.plt, aes(x = YSO, y = N.md)) +
  geom_ribbon(aes(ymin = N.95lo, ymax = N.95hi, fill = pctdead), alpha = 0.3) +
  geom_line(aes(color = pctdead), size = 2) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  xlab("Years since outbreak") + ylab(NULL) +
  labs(fill = "Dead\nconifer\n(%)") +
  labs(color = "Dead\nconifer\n(%)")

p.heat <- ggplot(plt.tbl, aes(YSO, pctdead)) +
  geom_tile(aes(fill = N.md), color = "white") +
  scale_fill_gradient(low = "#0072B2", high = "#D55E00") +
  ylab("Dead conifer (%)") +
  xlab("Years since outbreak") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(fill = "Median\ndensity\nestimate\n(per 100 ha)")

p.SF <- ggdraw() +
  draw_plot(p.pdead, x = 0, y = 0, width = 0.333333, height = 1) +
  draw_plot(p.YSO, x = 0.333333, y = 0, width = 0.333333, height = 1) +
  draw_plot(p.heat, x = 0.666667, y = 0, width = 0.333333, height = 1)

## Lodgepole pine stratum ##
stratum <- "LP" #Set to LP or SF
mod <- loadObject("mod_RESQ_outbreak_HZdist_LP_reduced")

# Plot outbreak relationships #
Cov <- eval(as.name(str_c("Cov.", stratum)))
x.pctdead <- seq(0, 1, length.out = 20)
x.YSO <- seq(0, max(Cov[, "YSO"], na.rm = T))
plt.tbl <- expand.grid(x.pctdead, x.YSO) %>%
  rename(pctdead = Var1, YSO = Var2) %>%
  mutate(pctdead.z = (pctdead - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)) %>%
  mutate(YSO.z = (YSO - mean(Cov[, "YSO"], na.rm = T)) / sd(Cov[, "YSO"], na.rm = T)) %>%
  mutate(pctdead = pctdead * 100)

pred.N <- matrix(NA, nrow = nrow(plt.tbl), ncol = mod$mcmc.info$n.samples)
for(j in 1:nrow(pred.N))
  pred.N[j, ] <- (exp(mod$sims.list$beta0.mean +
                        mod$sims.list$bl.pdead * plt.tbl$pctdead.z[j] +
                        mod$sims.list$bl.YSO * plt.tbl$YSO.z[j] +
                        mod$sims.list$bl.pdXYSO * plt.tbl$pctdead.z[j] * plt.tbl$YSO.z[j]) /
                    samp.ha) * 100
plt.tbl$N.md <- apply(pred.N, 1, median)
plt.tbl$N.95lo <- apply(pred.N, 1, function(x) quantile(x, prob = 0.025, type = 8))
plt.tbl$N.95hi <- apply(pred.N, 1, function(x) quantile(x, prob = 0.975, type = 8))

dat.plt <- plt.tbl %>% filter(YSO %in% c(min(YSO), max(YSO))) %>%
  mutate(YSO = as.factor(YSO))
p.pdead <- ggplot(dat.plt, aes(x = pctdead, y = N.md)) +
  geom_ribbon(aes(ymin = N.95lo, ymax = N.95hi, fill = YSO), alpha = 0.3) +
  geom_line(aes(color = YSO), size = 2) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  xlab("Dead conifer (%)") + ylab(NULL) +
  labs(fill = "Years\nsince\noutbreak") +
  labs(color = "Years\nsince\noutbreak")

dat.plt <- plt.tbl %>% filter(pctdead %in% c(min(pctdead), max(pctdead))) %>%
  mutate(pctdead = as.factor(pctdead))
p.YSO <- ggplot(dat.plt, aes(x = YSO, y = N.md)) +
  geom_ribbon(aes(ymin = N.95lo, ymax = N.95hi, fill = pctdead), alpha = 0.3) +
  geom_line(aes(color = pctdead), size = 2) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  xlab("Years since outbreak") + ylab(NULL) +
  labs(fill = "Dead\nconifer\n(%)") +
  labs(color = "Dead\nconifer\n(%)")

p.heat <- ggplot(plt.tbl, aes(YSO, pctdead)) +
  geom_tile(aes(fill = N.md), color = "white") +
  scale_fill_gradient(low = "#0072B2", high = "#D55E00") +
  ylab("Dead conifer (%)") +
  xlab("Years since outbreak") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(fill = "Median\ndensity\nestimate\n(per 100 ha)")

p.LP <- ggdraw() +
  draw_plot(p.pdead, x = 0, y = 0, width = 0.333333, height = 1) +
  draw_plot(p.YSO, x = 0.333333, y = 0, width = 0.333333, height = 1) +
  draw_plot(p.heat, x = 0.666667, y = 0, width = 0.333333, height = 1)

## Put em together ##
p <- ggdraw() +
  draw_plot(p.LP, x = 0.05, y = 0.5, width = 0.95, height = 0.45) +
  draw_plot(p.SF, x = 0.05, y = 0, width = 0.95, height = 0.45) +
  draw_plot_label(c("Density per 100 ha", "Lodgepole pine", "Spruce-fir"),
                  x = c(0.02, 0.5, 0.5), y = c(0.4, 0.98, 0.48),
                  angle = c(90, 0, 0), hjust = c(0, 0, 0))

save_plot(str_c("figure_RESQ_outbreak.tiff"), p, ncol = 3, nrow = 2, dpi = 200)
