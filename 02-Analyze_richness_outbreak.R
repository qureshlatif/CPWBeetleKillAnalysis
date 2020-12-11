library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#_______Script inputs_______#
stratum <- "LP"
mod <- loadObject("mod_LPcommunity_outbreak_reduced2")
maxyso <- 12 # Set to 12 for LP and 9 for SF
#___________________________#
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval

#### Plot species richness ####
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
x.yso <- -1:18
z.yso <- (x.yso - x.mn) / x.sd

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(YSO.check > maxyso)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
x.mn <- mean(PctDead.b, na.rm = T); x.sd <- sd(PctDead.b, na.rm = T)
PctDead.bz <- (PctDead.b - x.mn) / x.sd # Re-scale point values
PctDead.bz[which(is.na(PctDead.bz) & !outbreak_grids[ID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & !outbreak_grids[ID])])
PctDead.bz[which(is.na(PctDead.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.bz[which(!is.na(PctDead.bz) & outbreak_grids[ID])])
PctDead.b <- (PctDead.bz * x.sd) + x.mn
x.pd <- seq(0, 100, by = 10)
z.pd <- ((x.pd / 100) - x.mn) / x.mn

dat.pred.DCon <- data.frame(x.pd = x.pd, z.pd = z.pd)
dat.pred.YSO <- data.frame(x.yso = x.yso, z.yso = z.yso)

rm(x.mn, x.sd, x.yso, z.yso, x.pd, z.pd)

omega <- mod$sims.list$omega # Probability of species occurring in the super community
psi <- expit(mod$sims.list$d0)
b0 <- mod$sims.list$b0
b.pdead <- mod$sims.list$bb.pdead
b.yso <- mod$sims.list$bb.YSO
b.yso2 <- mod$sims.list$bb.YSO2

SR.pred.DCon <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred.DCon))
for(i in 1:dim(SR.pred.DCon)[2]) {
  theta <- expit(b0 + b.pdead*dat.pred.DCon$z.pd[i])
  SR.pred.DCon[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred.DCon <- dat.pred.DCon %>%
  mutate(pred.md = apply(SR.pred.DCon, 2, median),
         pred.lo = apply(SR.pred.DCon, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred.DCon, 2, function(x) quantile(x, prob = 0.95, type = 8)))
write.csv(dat.pred.DCon, str_c("Spp_richness_pred_DCon_cache_", stratum, ".csv"), row.names = F)

SR.pred.YSO <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred.YSO))
for(i in 1:dim(SR.pred.YSO)[2]) {
  theta <- expit(b0 + b.yso*dat.pred.YSO$z.yso[i] + b.yso2*(dat.pred.YSO$z.yso[i]^2))
  SR.pred.YSO[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred.YSO <- dat.pred.YSO %>%
  mutate(pred.md = apply(SR.pred.YSO, 2, median),
         pred.lo = apply(SR.pred.YSO, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred.YSO, 2, function(x) quantile(x, prob = 0.95, type = 8)))
write.csv(dat.pred.YSO, str_c("Spp_richness_pred_YSO_cache_", stratum, ".csv"), row.names = F)
