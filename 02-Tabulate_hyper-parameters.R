library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#___________________ Covariate relations only ______________________#
stratum <- c("LP", "LP", "SF", "SF") # Select LP or SF
mod.type <- c("outbreak", "habitat", "outbreak", "habitat")
mods <- c("mod_LPcommunity_outbreak_reduced2", "mod_LPcommunity_habitat_reduced",
          "mod_SFcommunity_outbreak_reduced2", "mod_SFcommunity_habitat_reduced")
params.mn <- c("Betab.PctDead", "Betab.YSO", "Betab.YSO2", "Betab.PctDdXYSO",
               "Betab.CanCov", "Betab.RCovAS", "Betab.RCovES",
               "Betab.RCovPine", "Betab.ShCov", "Betab.RSC_Con",
               "Betab.GHerb", "Betab.Gwoody", "Betab.GDD")
params.sd <- str_c("sigma.", params.mn)
cols <- expand.grid(str_c(stratum, mod.type, sep = "."), c("mean", "SD")) %>%
  select(Var2, Var1) %>%
  mutate(Var3 = str_c(Var2, Var1, sep = ".")) %>%
  pull(Var3)
cols <- cols[c(1, 5, 2, 6, 3, 7, 4, 8)]
out <- matrix(NA, nrow = length(params.mn), ncol = length(cols),
              dimnames = list(params.mn, cols))

for(m in 1:length(mods)) {
  mod <- loadObject(mods[m])
  for(p in which(params.mn %in% names(mod$sims.list))) {
    out[p, (m*2 - 1)] <- str_c(
      median(mod$sims.list[[params.mn[p]]]) %>% round(digits = 2),
      " (",
      quantile(mod$sims.list[[params.mn[p]]],
               prob = 0.05, type = 8) %>% round(digits = 2),
      ",",
      quantile(mod$sims.list[[params.mn[p]]],
               prob = 0.95, type = 8) %>% round(digits = 2),
      ")")
    out[p, m*2] <- str_c(
      median(mod$sims.list[[params.sd[p]]]) %>% round(digits = 2),
      " (",
      quantile(mod$sims.list[[params.sd[p]]],
               prob = 0.05, type = 8) %>% round(digits = 2),
      ",",
      quantile(mod$sims.list[[params.sd[p]]],
               prob = 0.95, type = 8) %>% round(digits = 2),
      ")")
  }
}

out %>% write.csv("Hyper_parameters.csv", row.names = T)
#__________________________________________________________________#

#____________________________ All _________________________________#
mods <- c("mod_LPcommunity_outbreak_reduced2", "mod_LPcommunity_habitat_reduced",
          "mod_SFcommunity_outbreak_reduced2", "mod_SFcommunity_habitat_reduced")
params <- c("omega", "rho.ab", "rho.bd", "alpha0", "sigma.a0",
            "beta0", "sigma.b0", "delta0", "sigma.d0",
            "Betab.PctDead", "sigma.Betab.PctDead", "Betab.YSO", "sigma.Betab.YSO",
            "Betab.YSO2", "sigma.Betab.YSO2", "Betab.PctDdXYSO", "sigma.Betab.PctDdXYSO",
            "Betab.CanCov", "sigma.Betab.CanCov", "Betab.RCovAS", "sigma.Betab.RCovAS",
            "Betab.RCovES", "sigma.Betab.RCovES", "Betab.RCovPine", "sigma.Betab.RCovPine",
            "Betab.ShCov", "sigma.Betab.ShCov", "Betab.RSC_Con", "sigma.Betab.RSC_Con",
            "Betab.GHerb", "sigma.Betab.GHerb", "Betab.Gwoody", "sigma.Betab.Gwoody",
            "Betab.GDD", "sigma.Betab.GDD", "Betaa.Time", "sigma.Betaa.Time",
            "Betaa.Time2", "sigma.Betaa.Time2", "Betaa.DOY", "sigma.Betaa.DOY",
            "Betaa.DOY2", "sigma.Betaa.DOY2", "Betaa.PctDead", "sigma.Betaa.PctDead",
            "Betaa.YSO", "sigma.Betaa.YSO", "Betaa.CCov", "sigma.Betaa.CCov",
            "Betaa.SHCov", "sigma.Betaa.SHCov")
out <- matrix(NA, nrow = length(params), ncol = length(mods),
              dimnames = list(params, mods))

for(m in 1:length(mods)) {
  mod <- loadObject(mods[m])
  for(p in which(params %in% names(mod$sims.list)))
    out[p, m] <- str_c(
      median(mod$sims.list[[params[p]]]) %>% round(digits = 2),
      " (",
      quantile(mod$sims.list[[params[p]]],
               prob = 0.05, type = 8) %>% round(digits = 2),
      ",",
      quantile(mod$sims.list[[params[p]]],
               prob = 0.95, type = 8) %>% round(digits = 2),
      ")")
}

out %>% write.csv("Hyper_parameters.csv", row.names = T)
#__________________________________________________________________#