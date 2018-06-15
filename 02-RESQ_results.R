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

#___ Inputs ___#
stratum <- "LP" # LP or SF for lodgepole or spruce-fire, respectively
det.model <- "HN" # Either HN for half-normal or HZ for hazard rate
#______________#

mod <- loadObject(str_c("mod_RESQ_allcovs_", det.model, "dist_", stratum))

# Parameter summary table #
sum.table <- mod$summary %>% tbl_df %>%
  mutate(parameter = dimnames(mod$summary)[[1]])
if(det.model == "HZ")
  sum.table <- sum.table %>%
  filter(parameter %in% c("beta0.mean", "beta0.sd", "bl.pdead", "bl.pdead2",
                          "bl.outbrk", "bl.YSO", "bl.YSO2", "bl.pdXYSO", "bl.TWIP",
                          "a0", "b", "a.Time", "a.Time2", "a.DOY", "a.DOY2",
                          "a.ccov", "a.shcov"))
if(det.model == "HN")
  sum.table <- sum.table %>%
  filter(parameter %in% c("beta0.mean", "beta0.sd", "bl.pdead", "bl.pdead2",
                          "bl.outbrk", "bl.YSO", "bl.YSO2", "bl.pdXYSO", "bl.TWIP",
                          "bt.0", "bt.Time", "bt.Time2", "bt.DOY", "bt.DOY2",
                          "bt.ccov", "bt.shcov"))
sum.table <- sum.table %>%
  mutate(estimate = str_c(round(`50%`, digits = 3), "(", round(`2.5%`, digits = 3), ",", round(`97.5%`, digits = 3), ")")) %>%
  select(parameter, estimate, Rhat, n.eff, overlap0, f) %>%
  mutate(n.eff = n.eff %>% as.integer)

write.csv(sum.table, str_c("Param_summ_RESQ_allcovs_", det.model,"dist", stratum, ".csv"), row.names = F)

# Plot outbreak relationships #
Cov <- eval(as.name(str_c("Cov.", stratum)))
x.pctdead <- seq(0, 1, length.out = 20)
#z.pctdead <- (x.pctdead - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)
x.YSO <- seq(0, max(Cov[, "YSO"], na.rm = T))
#z.YSO <- (x.YSO - mean(Cov[, "YSO"], na.rm = T)) / sd(Cov[, "YSO"], na.rm = T)
plt.tbl <- expand.grid(x.pctdead, x.YSO) %>%
  rename(pctdead = Var1, YSO = Var2) %>%
  mutate(pctdead.z = (pctdead - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)) %>%
  mutate(YSO.z = (YSO - mean(Cov[, "YSO"], na.rm = T)) / sd(Cov[, "YSO"], na.rm = T)) %>%
  mutate(Outbreak = 1)
plt.tbl <- plt.tbl %>% rbind(c(mean(Cov[-ind, "DeadConif"], na.rm = T),
                               -2,
                               ((mean(Cov[-ind, "DeadConif"], na.rm = T) - mean(Cov[, "DeadConif"], na.rm = T)) / sd(Cov[, "DeadConif"], na.rm = T)),
                               0, 0))

pred.N <- matrix(NA, nrow = nrow(plt.tbl), ncol = mod$mcmc.info$n.samples)
for(j in 1:nrow(pred.N))
  pred.N[j, ] <- exp(mod$sims.list$beta0.mean +
                       mod$sims.list$bl.pdead * plt.tbl$pctdead.z[j] +
                       mod$sims.list$bl.pdead2 * (plt.tbl$pctdead.z[j]^2) +
                       mod$sims.list$bl.outbrk * plt.tbl$Outbreak[j] +
                       mod$sims.list$bl.YSO * plt.tbl$YSO.z[j] +
                       mod$sims.list$bl.YSO2 * (plt.tbl$YSO[j]^2) +
                       mod$sims.list$bl.pdXYSO * plt.tbl$pctdead.z[j] * plt.tbl$YSO.z[j])
plt.tbl$N.md <- apply(pred.N, 1, median)
plt.tbl$N.95lo <- apply(pred.N, 1, function(x) quantile(x, prob = 0.025, type = 8))
plt.tbl$N.95hi <- apply(pred.N, 1, function(x) quantile(x, prob = 0.975, type = 8))
