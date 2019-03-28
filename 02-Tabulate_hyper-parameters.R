library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

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
