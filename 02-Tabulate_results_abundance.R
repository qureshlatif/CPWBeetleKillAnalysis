library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled_abundance.RData")

#____ Convergence ____#
cols <- c("Spp", "strata", "min_n.eff", "max_Rhat")
out <- matrix(NA, nrow = length(Spp_LP) + length(Spp_SF), ncol = length(cols))
dimnames(out) <- list(NULL, cols)

i <- 1
for(strt in (strata %>% str_sub(-2, -1))) {
  SPP <- str_c("Spp_", strt) %>% as.name %>% eval
  for(spp in SPP) {
    mod <- loadObject(str_c("abund_models/mod_", spp,"_abundance_outbreak_", strt))
    out[i, "Spp"] <- spp
    out[i, "strata"] <- strt
    out[i, "min_n.eff"] <- min(mod$summary[, "n.eff"])
    out[i, "max_Rhat"] <- max(mod$summary[, "Rhat"])
    i <- i + 1
    }
  }
rm(i)

write.csv(out, "abund_models/Convergence_summary.csv", row.names = F)
#____________________#

#_____Tabulate parameter estimates_____#
conv_sum <- out %>% data.frame(stringsAsFactors = F) %>%
  mutate(min_n.eff = min_n.eff %>% as.integer) %>%
  mutate(max_Rhat = max_Rhat %>% as.numeric)
Spp_LP_keep <- conv_sum %>%
  filter(strata == "LP" & min_n.eff >= 30) %>%
  pull(Spp)
Spp_SF_keep <- conv_sum %>%
  filter(strata == "SF" & min_n.eff >= 30) %>%
  pull(Spp)

cols <- c("Spp", "strata",
          "beta0.mean", "beta0.sd",
          "bl.pdead", "bl.YSO", "bl.YSO2", "bl.pdXYSO",
          "bl.RCovAS", "bl.RCovES", "bl.RCovPine",
          ##___ detection parameters ___##
          "a0", "a.Time", "a.Time2",
          "a.DOY", "a.DOY2",
          "a.pdead", "a.YSO",
          "b") # GOF
out <- matrix(NA, nrow = length(Spp_LP_keep) + length(Spp_SF_keep), ncol = length(cols))
dimnames(out) <- list(NULL, cols)

i <- 1
for(strt in (strata %>% str_sub(-2, -1))) {
  SPP <- str_c("Spp_", strt, "_keep") %>% as.name %>% eval
  for(spp in SPP) {
    mod <- loadObject(str_c("abund_models/mod_", spp,"_abundance_outbreak_", strt))
    out[i, "Spp"] <- spp
    out[i, "strata"] <- strt
    vars <- cols[-c(1:2)]
    vars <- cols[which(cols %in% names(mod$sims.list))]
    for(var in vars){
      if((mod$sims.list[[var]] %>% quantile(prob = 0.05, type = 8)) > 0 |
         (mod$sims.list[[var]] %>% quantile(prob = 0.95, type = 8)) < 0) {
        out[i, var] <- str_c(mod$summary[var, "50%"] %>% round(digits = 2),
                           "(",
                           mod$sims.list[[var]] %>%
                             quantile(prob = 0.05, type = 8) %>%
                             round(digits = 2),
                           ",",
                           mod$sims.list[[var]] %>%
                             quantile(prob = 0.95, type = 8) %>%
                             round(digits = 2),
                           ")*")} else {
                             out[i, var] <- str_c(mod$summary[var, "50%"] %>% round(digits = 2),
                                                  "(",
                                                  mod$sims.list[[var]] %>%
                                                    quantile(prob = 0.05, type = 8) %>%
                                                    round(digits = 2),
                                                  ",",
                                                  mod$sims.list[[var]] %>%
                                                    quantile(prob = 0.95, type = 8) %>%
                                                    round(digits = 2),
                                                  ")")
                           }
    }
    i <- i + 1
  }
}
rm(i)

write.csv(out, "Abundance_model_estimates.csv", row.names = F)


