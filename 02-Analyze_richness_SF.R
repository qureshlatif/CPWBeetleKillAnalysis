library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(lme4)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#_______Script inputs_______#
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced3")
maxyso <- 9 # Set to 12 for LP and 9 for SF
#___________________________#
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval

#### Plot species richness ####
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
  mutate(PctDead.z = (PctDead.x - x.mn) / x.sd,
         YSO.z = 0)

YSO.b <- Cov[, "YSO"]
x.mn <- mean(YSO.b, na.rm = T); x.sd <- sd(YSO.b, na.rm = T)
YSO.bz <- (YSO.b - x.mn) / x.sd # Re-scale point values
YSO.bz[which(!outbreak_grids[ID])] <- min(YSO.bz, na.rm = T) # Insert zero for non-outbreak grids
YSO.bz[which(is.na(YSO.bz) & outbreak_grids[ID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(YSO.bz[which(!is.na(YSO.bz) & outbreak_grids[ID])])
YSO.b <- (YSO.bz * x.sd) + x.mn
dat.pred.YSO <- data.frame(YSO.x = seq(min(YSO.b, na.rm = T), max(YSO.b, na.rm = T), by = 2) %>% round) %>%
  mutate(PctDead.z = 0,
         YSO.z = (YSO.x - x.mn) / x.sd)

rm(x.mn, x.sd)

SPR <- mod$sims.list$SR.point
# Derive posterior samples for plotting spp richness trend #
Y.PDd <- matrix(NA, nrow = dim(SPR)[1], ncol = nrow(dat.pred.PDd))
Y.YSO <- matrix(NA, nrow = dim(SPR)[1], ncol = nrow(dat.pred.YSO))
B0_pnt <- B_PDd <- B_YSO <- B_YSO2 <- B0_pntSD <- rep(NA, length = dim(SPR)[1])
options(warn=1) # Change to warn=2 to investigate convergence warnings.
for(i in 1:dim(SPR)[1]) {
  dat <- data.frame(y = SPR[i,],
                    PctDead.z = PctDead.bz,
                    YSO.z = YSO.bz,
                    ID = ID)
  m <- glmer(y ~ PctDead.z + YSO.z + I(YSO.z^2) + (1|ID), data = dat, family = "poisson")
  B0_pnt[i] <- summary(m)$coefficients["(Intercept)", "Estimate"]
  B_PDd[i] <- summary(m)$coefficients["PctDead.z", "Estimate"]
  B_YSO[i] <- summary(m)$coefficients["YSO.z", "Estimate"]
  B_YSO2[i] <- summary(m)$coefficients["I(YSO.z^2)", "Estimate"]
  B0_pntSD[i] <- as.data.frame(VarCorr(m))$sdcor
  Y.PDd[i, ] <- predict(m, dat.pred.PDd, re.form = NA, type = "response")
  Y.YSO[i, ] <- predict(m, dat.pred.YSO, re.form = NA, type = "response")
}
B0_pnt <- str_c(median(B0_pnt, na.rm = T) %>% round(digits = 5),
                " (",
                quantile(B0_pnt, prob = 0.05, type = 8, na.rm = T) %>% round(digits = 5),
                ",",
                quantile(B0_pnt, prob = 0.95, type = 8, na.rm = T) %>% round(digits = 5),
                ")")
B_PDd <- str_c(median(B_PDd, na.rm = T) %>% round(digits = 5),
               " (",
               quantile(B_PDd, prob = 0.05, type = 8, na.rm = T) %>% round(digits = 5),
               ",",
               quantile(B_PDd, prob = 0.95, type = 8, na.rm = T) %>% round(digits = 5),
               ")")
B_YSO <- str_c(median(B_YSO, na.rm = T) %>% round(digits = 5),
               " (",
               quantile(B_YSO, prob = 0.05, type = 8, na.rm = T) %>% round(digits = 5),
               ",",
               quantile(B_YSO, prob = 0.95, type = 8, na.rm = T) %>% round(digits = 5),
               ")")
B_YSO2 <- str_c(median(B_YSO2, na.rm = T) %>% round(digits = 5),
                " (",
                quantile(B_YSO2, prob = 0.05, type = 8, na.rm = T) %>% round(digits = 5),
                ",",
                quantile(B_YSO2, prob = 0.95, type = 8, na.rm = T) %>% round(digits = 5),
                ")")
B0_pntSD <- str_c(median(B0_pntSD, na.rm = T) %>% round(digits = 5),
                  " (",
                  quantile(B0_pntSD, prob = 0.05, type = 8, na.rm = T) %>% round(digits = 5),
                  ",",
                  quantile(B0_pntSD, prob = 0.95, type = 8, na.rm = T) %>% round(digits = 5),
                  ")")

dat.pred.PDd <- dat.pred.PDd %>%
  mutate(Y.md = apply(Y.PDd, 2, median),
         Y.lo = apply(Y.PDd, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y.PDd, 2, function(x) quantile(x, prob = 0.975, type = 8)))
dat.pred.YSO <- dat.pred.YSO %>%
  mutate(Y.md = apply(Y.YSO, 2, median),
         Y.lo = apply(Y.YSO, 2, function(x) quantile(x, prob = 0.025, type = 8)),
         Y.hi = apply(Y.YSO, 2, function(x) quantile(x, prob = 0.975, type = 8)))
write.csv(dat.pred.PDd, str_c("Spp_richness_datPDd_pred_cache_", stratum, ".csv"), row.names = F)
write.csv(dat.pred.YSO, str_c("Spp_richness_datYSO_pred_cache_", stratum, ".csv"), row.names = F)
saveObject(B0_pnt, str_c("B0_pnt_N_outbreak_cache_", stratum))
saveObject(B_PDd, str_c("B_PDd_pnt_N_outbreak_cache_", stratum))
saveObject(B_YSO, str_c("B_YSO_pnt_N_outbreak_cache_", stratum))
saveObject(B_YSO2, str_c("B_YSO2_pnt_N_outbreak_cache_", stratum))
saveObject(B0_pntSD, str_c("B0_pntSD_pnt_N_outbreak_cache_", stratum))
rm(i, m, Y.PDd, Y.YSO)
