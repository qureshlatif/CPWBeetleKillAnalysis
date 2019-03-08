library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(lme4)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#______Lodgepole pine stratum_______#
stratum <- "LP"
mod <- loadObject("mod_LPcommunity_outbreak_reduced3")
maxyso <- 12 # Set to 12 for LP and 9 for SF

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(Cov[, "YSO"] > maxyso)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
x.mn <- mean(PctDead.b, na.rm = T); x.sd <- sd(PctDead.b, na.rm = T)
PctDead.bz <- (PctDead.b - x.mn) / x.sd # Re-scale point values
PctDead.bz[which(is.na(PctDead.bz) & !outbreak_grids[ID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & !outbreak_grids[ID])])
PctDead.bz[which(is.na(PctDead.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & outbreak_grids[ID])])
PctDead.b <- (PctDead.bz * x.sd) + x.mn
dat.pred.PDd <- data.frame(PctDead.x = seq(min(PctDead.b, na.rm = T), max(PctDead.b, na.rm = T), length.out = 20)) %>%
  mutate(PctDead.z = (PctDead.x - mean(PctDead.b, na.rm = T)) / sd(PctDead.b, na.rm = T),
         YSO.z = 0)
YSO.b <- Cov[, "YSO"]
x.mn <- mean(YSO.b, na.rm = T); x.sd <- sd(YSO.b, na.rm = T)
YSO.bz <- (YSO.b - x.mn) / x.sd # Re-scale point values
YSO.bz[which(!outbreak_grids[ID])] <- min(YSO.bz, na.rm = T) # Insert zero for non-outbreak grids
YSO.bz[which(is.na(YSO.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(YSO.bz[which(!is.na(YSO.bz) & outbreak_grids[ID])])
YSO.b <- (YSO.bz * x.sd) + x.mn
dat.pred.YSO <- data.frame(YSO.x = seq(min(YSO.b, na.rm = T), max(YSO.b, na.rm = T), length.out = 20)) %>%
  mutate(PctDead.z = 0,
         YSO.z = (YSO.x - mean(YSO.b, na.rm = T)) / sd(YSO.b, na.rm = T))
rm(x.mn, x.sd)

SPR <- mod$sims.list$SR.point
dat.pred.PDd <- read.csv(str_c("Spp_richness_datPDd_pred_cache_", stratum, ".csv"), header = T)
dat.pred.YSO <- read.csv(str_c("Spp_richness_datYSO_pred_cache_", stratum, ".csv"), header = T)
B0_pnt <- loadObject(str_c("B0_pnt_N_outbreak_cache_", stratum))
B_PDd <- loadObject(str_c("B_PDd_pnt_N_outbreak_cache_", stratum))
B_YSO <- loadObject(str_c("B_YSO_pnt_N_outbreak_cache_", stratum))
B_YSO2 <- loadObject(str_c("B_YSO2_pnt_N_outbreak_cache_", stratum))
B0_pntSD <- loadObject(str_c("B0_pntSD_pnt_N_outbreak_cache_", stratum))

dat.SR <- data.frame(PctDead = PctDead.b, YSO = YSO.b) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

p.PctDead.LP <- ggplot(data = dat.SR %>% filter(!is.na(PctDead)), aes(x = PctDead, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.PDd, aes(x = PctDead.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred.PDd, aes(x = PctDead.x, y = Y.md), size = 1.5, color = "blue", inherit.aes = F) + 
  labs(x = NULL, y = NULL)

p.YSO.LP <- ggplot(data = dat.SR %>% filter(!is.na(YSO)), aes(x = YSO, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.YSO, aes(x = YSO.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred.YSO, aes(x = YSO.x, y = Y.md), size = 1.5, color = "blue", inherit.aes = F) + 
  scale_x_continuous(breaks = seq(0, 18, by = 2)) +
  #scale_x_continuous(breaks = seq(min(YSO.b, na.rm = T), max(YSO.b, na.rm = T), by = 2)) +
  labs(x = NULL, y = NULL)

#______Spruce-fir stratum_______#
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced3")
maxyso <- 9 # Set to 12 for LP and 9 for SF

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(Cov[, "YSO"] > maxyso)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
x.mn <- mean(PctDead.b, na.rm = T); x.sd <- sd(PctDead.b, na.rm = T)
PctDead.bz <- (PctDead.b - x.mn) / x.sd # Re-scale point values
PctDead.bz[which(is.na(PctDead.bz) & !outbreak_grids[ID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & !outbreak_grids[ID])])
PctDead.bz[which(is.na(PctDead.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & outbreak_grids[ID])])
PctDead.b <- (PctDead.bz * x.sd) + x.mn
dat.pred.PDd <- data.frame(PctDead.x = seq(min(PctDead.b, na.rm = T), max(PctDead.b, na.rm = T), length.out = 20)) %>%
  mutate(PctDead.z = (PctDead.x - mean(PctDead.b, na.rm = T)) / sd(PctDead.b, na.rm = T),
         YSO.z = 0)
YSO.b <- Cov[, "YSO"]
x.mn <- mean(YSO.b, na.rm = T); x.sd <- sd(YSO.b, na.rm = T)
YSO.bz <- (YSO.b - x.mn) / x.sd # Re-scale point values
YSO.bz[which(!outbreak_grids[ID])] <- min(YSO.bz, na.rm = T) # Insert zero for non-outbreak grids
YSO.bz[which(is.na(YSO.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(YSO.bz[which(!is.na(YSO.bz) & outbreak_grids[ID])])
YSO.b <- (YSO.bz * x.sd) + x.mn
dat.pred.YSO <- data.frame(YSO.x = seq(min(YSO.b, na.rm = T), max(YSO.b, na.rm = T), length.out = 20)) %>%
  mutate(PctDead.z = 0,
         YSO.z = (YSO.x - mean(YSO.b, na.rm = T)) / sd(YSO.b, na.rm = T))
rm(x.mn, x.sd)

SPR <- mod$sims.list$SR.point
dat.pred.PDd <- read.csv(str_c("Spp_richness_datPDd_pred_cache_", stratum, ".csv"), header = T)
dat.pred.YSO <- read.csv(str_c("Spp_richness_datYSO_pred_cache_", stratum, ".csv"), header = T)
B0_pnt <- loadObject(str_c("B0_pnt_N_outbreak_cache_", stratum))
B_PDd <- loadObject(str_c("B_PDd_pnt_N_outbreak_cache_", stratum))
B_YSO <- loadObject(str_c("B_YSO_pnt_N_outbreak_cache_", stratum))
B_YSO2 <- loadObject(str_c("B_YSO2_pnt_N_outbreak_cache_", stratum))
B0_pntSD <- loadObject(str_c("B0_pntSD_pnt_N_outbreak_cache_", stratum))

dat.SR <- data.frame(PctDead = PctDead.b, YSO = YSO.b) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

p.PctDead.SF <- ggplot(data = dat.SR %>% filter(!is.na(PctDead)), aes(x = PctDead, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.PDd, aes(x = PctDead.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred.PDd, aes(x = PctDead.x, y = Y.md), size = 1.5, color = "blue", inherit.aes = F) + 
  labs(x = NULL, y = NULL)

p.YSO.SF <- ggplot(data = dat.SR %>% filter(!is.na(YSO)), aes(x = YSO, y = Y)) + 
  geom_point(alpha = 0.3) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.3) +
  geom_ribbon(data = dat.pred.YSO, aes(x = YSO.x, ymin = Y.lo, ymax = Y.hi), alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred.YSO, aes(x = YSO.x, y = Y.md), size = 1.5, color = "blue", inherit.aes = F) + 
  scale_x_continuous(breaks = seq(min(YSO.b, na.rm = T), max(YSO.b, na.rm = T), by = 2)) +
  labs(x = NULL, y = NULL)

#_____ Put them all together _____#
p <- ggdraw() + 
  draw_plot(p.PctDead.LP, x = 0.05, y = 0.525, width = 0.475, height = 0.475) +
  draw_plot(p.YSO.LP, x = 0.525, y = 0.525, width = 0.475, height = 0.475) +
  draw_plot(p.PctDead.SF, x = 0.05, y = 0.05, width = 0.475, height = 0.475) +
  draw_plot(p.YSO.SF, x = 0.525, y = 0.05, width = 0.475, height = 0.475) +
  draw_plot_label(c("N[spruce-fir]", "N[lodgepole]", "Percent~dead~conifer", "Years~since~outbreak"),
                  x = c(0, 0, 0.2, 0.65), y = c(0.25, 0.7, 0.05, 0.05),
                  size = c(20, 20, 17, 17), angle = c(90, 90, 0, 0),
                  hjust = c(0, 0, 0, 0), parse = T)

save_plot("Plot_richness_outbreak.tiff", p, ncol = 2, nrow = 2, dpi = 200)
