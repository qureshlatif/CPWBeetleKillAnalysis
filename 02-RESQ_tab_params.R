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
mod <- loadObject("mod_RESQ_habitat_HZdist_SF")
pars <- c("beta0.mean", "beta0.sd", # Parameters of interest
          "bd.TWIP", "bd.heatload", "bd.TWI", "bd.RDens", "bd.WILD",
          "bl.pdead", "bl.YSO", "bl.YSO2", "bl.pdXYSO",
          "bl.ccov", "bl.RCovAS", "bl.RCovES", "bl.RCovPine",
          "bl.shcov", "bl.RSC_Con", "bl.GHerb", "bl.Gwoody", "bl.GDD",
          "a0", "b", "a.Time", "a.Time2", "a.DOY", "a.DOY2",
          "a.pdead", "a.YSO", "a.YSO2", "a.pdXYSO",
          "a.ccov", "a.shcov",
          "bt.0", "bt.Time", "bt.Time2", "bt.DOY",
          "bt.DOY2", "bt.ccov", "bt.shcov")
tab.out <- "Param_summ_RESQ.csv"
#______________#

# Parameter summary table #
sum.table <- mod$summary %>% tbl_df %>%
  mutate(parameter = dimnames(mod$summary)[[1]]) %>%
  filter(parameter %in% pars) %>%
  mutate(estimate = str_c(round(`50%`, digits = 3), "(", round(`2.5%`, digits = 3), ",", round(`97.5%`, digits = 3), ")")) %>%
  select(parameter, estimate, Rhat, n.eff, overlap0, f) %>%
  mutate(n.eff = n.eff %>% as.integer)

write.csv(sum.table, tab.out, row.names = F)

### For MS ###
mods <- c("mod_RESQ_outbreak_HZdist_LP", "mod_RESQ_outbreak_HZdist_SF",
          "mod_RESQ_habitat_HZdist_LP", "mod_RESQ_habitat_HZdist_SF")
mod.nams <- c("LP_outbrk", "SF_outbrk", "LP_hab", "SF_hab")
sum.tab <- matrix(NA, nrow = length(pars), ncol = length(mods),
                  dimnames = list(pars, mod.nams))

for(j in 1:length(mods)) {
  mod <- loadObject(mods[j])
  for(i in 1:length(pars)) {
    if(any(names(mod$sims.list) == pars[i])) {
      p <- mod$sims.list[[pars[i]]]
      sum.tab[i, j] <- str_c(median(p) %>% round(digits = 2),
                             "(",
                             quantile(p, prob = 0.05, type = 8) %>% round(digits = 2),
                             ",",
                             quantile(p, prob = 0.95, type = 8) %>% round(digits = 2),
                             ")")
    }
  }}

write.csv(sum.tab, tab.out, row.names = T)
