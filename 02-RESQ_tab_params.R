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
stratum <- "SF" #Set to LP or SF
mod <- loadObject("mod_RESQ_outbreak_HZdist_SF_reduced")
pars <- c("beta0.mean", "beta0.sd", # Parameters of interest
          "bd.TWIP", "bd.heatload", "bd.TWI", "bd.RDens", "bd.WILD",
          "bl.pdead", "bl.YSO", "bl.YSO2", "bl.pdXYSO",
          "bl.RCovAS", "bl.RCovES", "bl.RCovPine",
          "a0", "b", "a.Time", "a.Time2", "a.DOY", "a.DOY2",
          "a.pdead", "a.YSO", "a.YSO2", "a.pdXYSO",
          "bt.0", "bt.Time", "bt.Time2", "bt.DOY",
          "bt.DOY2", "bt.ccov", "bt.shcov")
tab.out <- "Param_summ_RESQ.csv"
samp.ha <- sum(area.band) * 0.0001
#______________#

# Parameter summary table #
sum.table <- mod$summary %>% tbl_df %>%
  mutate(parameter = dimnames(mod$summary)[[1]]) %>%
  filter(parameter %in% pars) %>%
  mutate(estimate = str_c(round(`50%`, digits = 3), "(", round(`2.5%`, digits = 3), ",", round(`97.5%`, digits = 3), ")")) %>%
  select(parameter, estimate, Rhat, n.eff, overlap0, f) %>%
  mutate(n.eff = n.eff %>% as.integer)

write.csv(sum.table, tab.out, row.names = F)
