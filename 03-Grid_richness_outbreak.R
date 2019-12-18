library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#______Lodgepole pine stratum_______#
stratum <- "LP"
mod <- loadObject("mod_LPcommunity_outbreak_reduced2")

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
SPR <- mod$sims.list$SR.grid

# Diff between outbreak and non-outbreak grids #
SPR %>% apply(1, (function(x) tapply(x, outbreak_grids, mean))) %>%
  apply(1, function(x) str_c(median(x) %>% round(digits = 2),
                             "(",
                             quantile(x, type = 8, prob = 0.05) %>% round(digits = 2),
                             ",",
                             quantile(x, type = 8, prob = 0.95) %>% round(digits = 2),
                             ")"))

x <- SPR %>% apply(1, (function(x) tapply(x, outbreak_grids, mean)))
x <- x[2, ] - x[1, ]
str_c(median(x) %>% round(digits = 2),
      "(",
      quantile(x, type = 8, prob = 0.05) %>% round(digits = 2),
      ",",
      quantile(x, type = 8, prob = 0.95) %>% round(digits = 2),
      ")")

#______Spruce-fir stratum_______#
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced2")

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
SPR <- mod$sims.list$SR.grid

# Diff between outbreak and non-outbreak grids #
SPR %>% apply(1, (function(x) tapply(x, outbreak_grids, mean))) %>%
  apply(1, function(x) str_c(median(x) %>% round(digits = 2),
                             "(",
                             quantile(x, type = 8, prob = 0.05) %>% round(digits = 2),
                             ",",
                             quantile(x, type = 8, prob = 0.95) %>% round(digits = 2),
                             ")"))

x <- SPR %>% apply(1, (function(x) tapply(x, outbreak_grids, mean)))
x <- x[2, ] - x[1, ]
str_c(median(x) %>% round(digits = 2),
      "(",
      quantile(x, type = 8, prob = 0.05) %>% round(digits = 2),
      ",",
      quantile(x, type = 8, prob = 0.95) %>% round(digits = 2),
      ")")
