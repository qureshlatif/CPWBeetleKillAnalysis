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

maxD <- 460 # upper limit for y-axis
maxYSO <- 15

## Spruce-fir stratum ##
stratum <- "SF" #Set to LP or SF
mod <- loadObject("mod_RESQ_outbreak_HZdist_SF")

# Plot outbreak relationships #
Cov <- eval(as.name(str_c("Cov.", stratum)))
x.pctdead <- seq(0, 1, length.out = 20)
x.YSO <- seq(-1, maxYSO)
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
                       mod$sims.list$bl.YSO2 * (plt.tbl$YSO.z[j]^2) +
                       mod$sims.list$bl.pdXYSO * plt.tbl$pctdead.z[j] * plt.tbl$YSO.z[j]) /
    samp.ha) * 100
plt.tbl$N.md <- apply(pred.N, 1, median)
plt.tbl$N.95lo <- apply(pred.N, 1, function(x) quantile(x, prob = 0.05, type = 8))
plt.tbl$N.95hi <- apply(pred.N, 1, function(x) quantile(x, prob = 0.95, type = 8))

dat.plt <- plt.tbl %>% filter(pctdead %in% c(min(pctdead), max(pctdead))) %>%
  mutate(pctdead = as.factor(pctdead))
p.SF <- ggplot(dat.plt, aes(x = YSO, y = N.md)) +
  geom_ribbon(aes(ymin = N.95lo, ymax = N.95hi, fill = pctdead), alpha = 0.3) +
  ylim(0, maxD) +
  geom_line(aes(color = pctdead), size = 2) +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  annotate("text", x = -0.8, y = (maxD - 50), label = "No outbreak", angle = 90) +
  scale_color_manual(values = c("#009E73", "#D55E00")) +
  scale_fill_manual(values = c("#009E73", "#D55E00")) +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "DCon") +
  labs(color = "DCon")

## Lodgepole pine stratum ##
stratum <- "LP" #Set to LP or SF
mod <- loadObject("mod_RESQ_outbreak_HZdist_LP")

# Plot outbreak relationships #
Cov <- eval(as.name(str_c("Cov.", stratum)))
x.pctdead <- seq(0, 1, length.out = 20)
x.YSO <- seq(-1, maxYSO)
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
                        mod$sims.list$bl.YSO2 * (plt.tbl$YSO.z[j] ^2) +
                        mod$sims.list$bl.pdXYSO * plt.tbl$pctdead.z[j] * plt.tbl$YSO.z[j]) /
                    samp.ha) * 100
plt.tbl$N.md <- apply(pred.N, 1, median)
plt.tbl$N.95lo <- apply(pred.N, 1, function(x) quantile(x, prob = 0.05, type = 8))
plt.tbl$N.95hi <- apply(pred.N, 1, function(x) quantile(x, prob = 0.95, type = 8))

dat.plt <- plt.tbl %>% filter(pctdead %in% c(min(pctdead), max(pctdead))) %>%
  mutate(pctdead = as.factor(pctdead))
p.LP <- ggplot(dat.plt, aes(x = YSO, y = N.md)) +
  ylim(0, maxD) +
  geom_ribbon(aes(ymin = N.95lo, ymax = N.95hi, fill = pctdead), alpha = 0.3) +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  annotate("text", x = -0.8, y = (maxD - 50), label = "No outbreak", angle = 90) +
  geom_line(aes(color = pctdead), size = 2) +
  scale_color_manual(values = c("#009E73", "#D55E00")) +
  scale_fill_manual(values = c("#009E73", "#D55E00")) +
  xlab(NULL) + ylab(NULL) +
  guides(fill = F, color = F)

## Put em together ##
p <- ggdraw() +
  draw_plot(p.LP, x = 0.05, y = 0.05, width = 0.45, height = 0.9) +
  draw_plot(p.SF, x = 0.5, y = 0.05, width = 0.5, height = 0.9) +
  draw_plot_label(c("Density per 100 ha", "Years since outbreak", "Lodgepole pine", "Spruce-fir"),
                  x = c(0.02, 0.4, 0.22, 0.68), y = c(0.3, 0.06, 0.98, 0.98),
                  angle = c(90, 0, 0, 0), hjust = c(0, 0, 0, 0))

save_plot(str_c("figure_RESQ_outbreak_MS.tiff"), p, ncol = 2.5, nrow = 1.25, dpi = 200)
