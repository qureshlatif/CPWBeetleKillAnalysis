library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

# Habitat models #
mod_LP <- loadObject("mod_RESQ_habitat_HZdist_LP")
mod_SF <- loadObject("mod_RESQ_habitat_HZdist_SF")
samp.ha <- sum(area.band) * 0.0001

dat.mnN.LP <- data.frame(avg.md = (median(mod_LP$sims.list$beta0.mean) %>%
                                 exp) / samp.ha * 100) %>%
  mutate(avg.lo = (quantile((mod_LP$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.05, type = 8)),
         avg.hi = (quantile((mod_LP$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.95, type = 8)),
         lo.md = (median(mod_LP$sims.list$beta0.mean -
                        1.96*mod_LP$sims.list$beta0.sd) %>%
                 exp) / samp.ha * 100,
         hi.md = (median(mod_LP$sims.list$beta0.mean +
                        1.96*mod_LP$sims.list$beta0.sd) %>%
                 exp) / samp.ha * 100)
dat.mnN.SF <- data.frame(avg.md = (median(mod_SF$sims.list$beta0.mean) %>%
                                     exp) / samp.ha * 100) %>%
  mutate(avg.lo = (quantile((mod_SF$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.05, type = 8)),
         avg.hi = (quantile((mod_SF$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.95, type = 8)),
         lo.md = (median(mod_SF$sims.list$beta0.mean -
                           1.96*mod_SF$sims.list$beta0.sd) %>%
                    exp) / samp.ha * 100,
         hi.md = (median(mod_SF$sims.list$beta0.mean +
                           1.96*mod_SF$sims.list$beta0.sd) %>%
                    exp) / samp.ha * 100)
dat.mnN.LP
dat.mnN.SF

# Habitat models #
mod_LP <- loadObject("mod_RESQ_outbreak_HZdist_LP")
mod_SF <- loadObject("mod_RESQ_outbreak_HZdist_SF")
samp.ha <- sum(area.band) * 0.0001

dat.mnN.LP <- data.frame(avg.md = (median(mod_LP$sims.list$beta0.mean) %>%
                                     exp) / samp.ha * 100) %>%
  mutate(avg.lo = (quantile((mod_LP$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.05, type = 8)),
         avg.hi = (quantile((mod_LP$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.95, type = 8)),
         lo.md = (median(mod_LP$sims.list$beta0.mean -
                           1.96*mod_LP$sims.list$beta0.sd) %>%
                    exp) / samp.ha * 100,
         hi.md = (median(mod_LP$sims.list$beta0.mean +
                           1.96*mod_LP$sims.list$beta0.sd) %>%
                    exp) / samp.ha * 100)
dat.mnN.SF <- data.frame(avg.md = (median(mod_SF$sims.list$beta0.mean) %>%
                                     exp) / samp.ha * 100) %>%
  mutate(avg.lo = (quantile((mod_SF$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.05, type = 8)),
         avg.hi = (quantile((mod_SF$sims.list$beta0.mean %>% exp) / samp.ha * 100,
                            prob = 0.95, type = 8)),
         lo.md = (median(mod_SF$sims.list$beta0.mean -
                           1.96*mod_SF$sims.list$beta0.sd) %>%
                    exp) / samp.ha * 100,
         hi.md = (median(mod_SF$sims.list$beta0.mean +
                           1.96*mod_SF$sims.list$beta0.sd) %>%
                    exp) / samp.ha * 100)
dat.mnN.LP
dat.mnN.SF
