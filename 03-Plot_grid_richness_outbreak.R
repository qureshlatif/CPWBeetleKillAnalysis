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

YSO.b <- Cov[, "YSO"]
YSO.b[which(is.na(YSO.b))] <- -1
YSO.d <- tapply(YSO.b, ID, sd) / mean(YSO.b)

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[is.na(PctDead.b)] <- mean(PctDead.b, na.rm = T)
PctDead.d <- tapply(PctDead.b, ID, sd) / mean(PctDead.b)

x.het <- YSO.d + PctDead.d

SPR <- mod$sims.list$SR.grid
dat.SR <- data.frame(x = x.het, outbreak_grids) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

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

# Relation with heterogeneity #
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = 10)
B0 <- B1 <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  dat = data.frame(Y = SPR[i, ], x = x.het)
  m <- suppressWarnings(glm(Y ~ x, data = dat, family = "poisson"))
  mn <- m$coefficients
  vc <- vcov(m)
  cfs <- MASS::mvrnorm(1, mn, vc)
  B0[i] <- cfs["(Intercept)"]
  B1[i] <- cfs["x"]
  X <- cbind(1, seq(min(x.het), max(x.het), length.out = 10)) %>% as.matrix
  Y[i, ] <- exp(X %*% cfs)
}
dat.pred <- data.frame(x.het = seq(min(x.het), max(x.het), length.out = 10)) %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))

B0 <- str_c(median(B0) %>% round(digits = 2),
                 " (",
                 quantile(B0, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B0, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")
B1 <- str_c(median(B1) %>% round(digits = 2),
                 " (",
                 quantile(B1, prob = 0.025, type = 8) %>% round(digits = 2),
                 ",",
                 quantile(B1, prob = 0.975, type = 8) %>% round(digits = 2),
                 ")")

p.LP.het <- ggplot(data = dat.SR, aes(x = x, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = x.het, ymin = Y.lo, ymax = Y.hi),
             alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x.het, y = Y.md),
            size = 1.5, inherit.aes = F) + 
  labs(x = NULL, y = NULL)

#______Spruce-fir stratum_______#
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced2")

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks

YSO.b <- Cov[, "YSO"]
YSO.b[which(is.na(YSO.b))] <- -1
YSO.d <- tapply(YSO.b, ID, sd) / mean(YSO.b)

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[is.na(PctDead.b)] <- mean(PctDead.b, na.rm = T)
PctDead.d <- tapply(PctDead.b, ID, sd) / mean(PctDead.b)

x.het <- YSO.d + PctDead.d

SPR <- mod$sims.list$SR.grid
dat.SR <- data.frame(x = x.het, outbreak_grids) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

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

# Relation with heterogeneity #
Y <- matrix(NA, nrow = dim(SPR)[1], ncol = 10)
B0 <- B1 <- numeric(length = dim(SPR)[1])
for(i in 1:dim(SPR)[1]) {
  dat = data.frame(Y = SPR[i, ], x = x.het)
  m <- suppressWarnings(glm(Y ~ x, data = dat, family = "poisson"))
  mn <- m$coefficients
  vc <- vcov(m)
  cfs <- MASS::mvrnorm(1, mn, vc)
  B0[i] <- cfs["(Intercept)"]
  B1[i] <- cfs["x"]
  X <- cbind(1, seq(min(x.het), max(x.het), length.out = 10)) %>% as.matrix
  Y[i, ] <- exp(X %*% cfs)
}
dat.pred <- data.frame(x.het = seq(min(x.het), max(x.het), length.out = 10)) %>%
  mutate(Y.md = apply(Y, 2, median),
         Y.lo = apply(Y, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y, 2, function(x) quantile(x, prob = 0.975, type = 8)))

B0 <- str_c(median(B0) %>% round(digits = 2),
            " (",
            quantile(B0, prob = 0.025, type = 8) %>% round(digits = 2),
            ",",
            quantile(B0, prob = 0.975, type = 8) %>% round(digits = 2),
            ")")
B1 <- str_c(median(B1) %>% round(digits = 2),
            " (",
            quantile(B1, prob = 0.025, type = 8) %>% round(digits = 2),
            ",",
            quantile(B1, prob = 0.975, type = 8) %>% round(digits = 2),
            ")")

p.SF <- ggplot(data = dat.SR, aes(x = x, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = x.het, ymin = Y.lo, ymax = Y.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x.het, y = Y.md),
            size = 1.5, inherit.aes = F) + 
  labs(x = NULL, y = NULL)

# #_____ Put them all together _____#
# p <- ggdraw() + 
#   draw_plot(p.DCon.LP, x = 0.05, y = 0.525, width = 0.475, height = 0.425) +
#   draw_plot(p.YSO.LP, x = 0.525, y = 0.525, width = 0.475, height = 0.425) +
#   draw_plot(p.DCon.SF, x = 0.05, y = 0.05, width = 0.475, height = 0.425) +
#   draw_plot(p.YSO.SF, x = 0.525, y = 0.05, width = 0.475, height = 0.425) +
#   draw_plot_label(c("Bird species richness", "Dead conifer", "Years since outbreak", "Lodgepole", "Spruce-fir"),
#                   x = c(0, 0.25, 0.7, 0.5, 0.5), y = c(0.4, 0.05, 0.05, 0.98, 0.5),
#                   size = c(20, 17, 17, 17, 17), angle = c(90, 0, 0, 0, 0),
#                   hjust = c(0, 0, 0, 0, 0))
# 
# save_plot("Plot_richness_outbreak.tiff", p, ncol = 2.5, nrow = 2.5, dpi = 200)
