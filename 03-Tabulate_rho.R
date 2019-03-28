require(dplyr)
require(stringr)
require(ggplot2)
require(cowplot)
require(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

mods <- c("mod_LPcommunity_outbreak_reduced2", "mod_LPcommunity_habitat_reduced", "mod_SFcommunity_outbreak_reduced2", "mod_SFcommunity_habitat_reduced")
rows <- c("LPoutbreak", "LPhabitat", "SFoutbreak", "SFhabitat")
cols <- c("rho.ab", "rho.bd")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

for(i in 1:length(rows)) {
  mod <- loadObject(mods[i])
  out[i, "rho.ab"] <- str_c(mod$sims.list$rho.ab %>% median %>% round(digits = 2),
                            " (",
                            mod$sims.list$rho.ab %>% quantile(prob = 0.05, type = 8) %>% round(digits = 2),
                            ",",
                            mod$sims.list$rho.ab %>% quantile(prob = 0.95, type = 8) %>% round(digits = 2),
                            ")")
  out[i, "rho.bd"] <- str_c(mod$sims.list$rho.bd %>% median %>% round(digits = 2),
                            " (",
                            mod$sims.list$rho.bd %>% quantile(prob = 0.05, type = 8) %>% round(digits = 2),
                            ",",
                            mod$sims.list$rho.bd %>% quantile(prob = 0.95, type = 8) %>% round(digits = 2),
                            ")")
}
out %>% write.csv("rho.csv", row.names = T)
