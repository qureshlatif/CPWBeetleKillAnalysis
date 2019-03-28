library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#______Lodgepole pine stratum_______#
stratum <- "LP"
mod <- loadObject("mod_LPcommunity_outbreak_reduced2")
maxyso <- 12 # Set to 12 for LP and 9 for SF

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks

YSO.b <- Cov[, "YSO"]
YSO.b[which(!outbreak_grids[ID])] <- -1
YSO.check <- YSO.b %>% # This is for screening out late-outbreak PctDead values
  replace(which(!outbreak_grids[ID]), -1) %>%
  replace(., which(is.na(.)), tapply(., ID, median, na.rm = T)[ID][which(is.na(.))] %>% round)
x.mn <- mean(YSO.b, na.rm = T); x.sd <- sd(YSO.b, na.rm = T)
YSO.bz <- (YSO.b - x.mn) / x.sd # Re-scale point values
YSO.bz[which(is.na(YSO.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(YSO.bz[which(!is.na(YSO.bz) & outbreak_grids[ID])])
YSO.b <- (YSO.bz * x.sd) + x.mn

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(YSO.check > maxyso)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
x.mn <- mean(PctDead.b, na.rm = T); x.sd <- sd(PctDead.b, na.rm = T)
PctDead.bz <- (PctDead.b - x.mn) / x.sd # Re-scale point values
PctDead.bz[which(is.na(PctDead.bz) & !outbreak_grids[ID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & !outbreak_grids[ID])])
PctDead.bz[which(is.na(PctDead.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & outbreak_grids[ID])])
PctDead.b <- (PctDead.bz * x.sd) + x.mn

SPR <- mod$sims.list$SR.point
dat.SR <- data.frame(PctDead = PctDead.b * 100, YSO = YSO.b + runif(length(YSO.b), -0.4, 0.4)) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

dat.pred <- read.csv(str_c("Spp_richness_pred_DCon_cache_", stratum, ".csv"), header = T)
p.DCon.LP <- ggplot(data = dat.SR, aes(x = PctDead, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = x.pd, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x.pd, y = pred.md),
            size = 1.5, inherit.aes = F) + 
  labs(x = NULL, y = NULL)

dat.pred <- read.csv(str_c("Spp_richness_pred_YSO_cache_", stratum, ".csv"), header = T)
p.YSO.LP <- ggplot(data = dat.SR, aes(x = YSO, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = x.yso, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x.yso, y = pred.md),
            size = 1.5, inherit.aes = F) + 
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  annotate("text", x = -1, y = 16, label = "No outbreak", angle = 90) +
  scale_x_continuous(breaks = seq(0, 18, 2)) + 
  labs(x = NULL, y = NULL)

#______Spruce-fir stratum_______#
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced2")
maxyso <- 9 # Set to 12 for LP and 9 for SF

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks

YSO.b <- Cov[, "YSO"]
YSO.b[which(!outbreak_grids[ID])] <- -1
YSO.check <- YSO.b %>% # This is for screening out late-outbreak PctDead values
  replace(which(!outbreak_grids[ID]), -1) %>%
  replace(., which(is.na(.)), tapply(., ID, median, na.rm = T)[ID][which(is.na(.))] %>% round)
x.mn <- mean(YSO.b, na.rm = T); x.sd <- sd(YSO.b, na.rm = T)
YSO.bz <- (YSO.b - x.mn) / x.sd # Re-scale point values
YSO.bz[which(is.na(YSO.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(YSO.bz[which(!is.na(YSO.bz) & outbreak_grids[ID])])
YSO.b <- (YSO.bz * x.sd) + x.mn

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(YSO.check > maxyso)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
x.mn <- mean(PctDead.b, na.rm = T); x.sd <- sd(PctDead.b, na.rm = T)
PctDead.bz <- (PctDead.b - x.mn) / x.sd # Re-scale point values
PctDead.bz[which(is.na(PctDead.bz) & !outbreak_grids[ID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & !outbreak_grids[ID])])
PctDead.bz[which(is.na(PctDead.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & outbreak_grids[ID])])
PctDead.b <- (PctDead.bz * x.sd) + x.mn

SPR <- mod$sims.list$SR.point
dat.SR <- data.frame(PctDead = PctDead.b * 100, YSO = YSO.b + runif(length(YSO.b), -0.4, 0.4)) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

dat.pred <- read.csv(str_c("Spp_richness_pred_DCon_cache_", stratum, ".csv"), header = T)
p.DCon.SF <- ggplot(data = dat.SR, aes(x = PctDead, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = x.pd, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x.pd, y = pred.md),
            size = 1.5, inherit.aes = F) + 
  labs(x = NULL, y = NULL)

dat.pred <- read.csv(str_c("Spp_richness_pred_YSO_cache_", stratum, ".csv"), header = T)
p.YSO.SF <- ggplot(data = dat.SR, aes(x = YSO, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred, aes(x = x.yso, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x.yso, y = pred.md),
            size = 1.5, inherit.aes = F) + 
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  annotate("text", x = -1, y = 16, label = "No outbreak", angle = 90) +
  scale_x_continuous(breaks = seq(0, 18, 2)) + 
  labs(x = NULL, y = NULL)

#_____ Put them all together _____#
p <- ggdraw() + 
  draw_plot(p.DCon.LP, x = 0.05, y = 0.525, width = 0.475, height = 0.425) +
  draw_plot(p.YSO.LP, x = 0.525, y = 0.525, width = 0.475, height = 0.425) +
  draw_plot(p.DCon.SF, x = 0.05, y = 0.05, width = 0.475, height = 0.425) +
  draw_plot(p.YSO.SF, x = 0.525, y = 0.05, width = 0.475, height = 0.425) +
  draw_plot_label(c("Bird species richness", "Dead conifer", "Years since outbreak", "Lodgepole", "Spruce-fir"),
                  x = c(0, 0.25, 0.7, 0.5, 0.5), y = c(0.4, 0.05, 0.05, 0.98, 0.5),
                  size = c(20, 17, 17, 17, 17), angle = c(90, 0, 0, 0, 0),
                  hjust = c(0, 0, 0, 0, 0))

save_plot("Plot_richness_outbreak.tiff", p, ncol = 2.5, nrow = 2.5, dpi = 200)
