library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#______Lodgepole pine stratum_______#
stratum <- "LP"
mod <- loadObject("mod_LPcommunity_habitat_reduced")

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor
miny <- 0; maxy <- 20

# Compile covariates #
ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.d <- tapply(ccov.b, ID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, ID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- ccov.d[ID[is.na(ccov.b)]]
rm(ccov.d, ccov.sd, ccov.b.missing)

RCovAS.b <- Cov[, "RCOV_AS"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovAS.d <- tapply(RCovAS.b, ID, mean, na.rm = T) # Grid-level values
RCovAS.b.missing <- is.na(RCovAS.b) %>% as.integer # Index missing values to be imputed
RCovAS.sd <- tapply(RCovAS.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovAS.sd==0))
  RCovAS.sd[which(RCovAS.sd==0)] <- min(RCovAS.sd[which(RCovAS.sd>0)])
RCovAS.lower <- min(RCovAS.b, na.rm = T) # Lower bound for imputing missing values
RCovAS.b[which(is.na(RCovAS.b))] <- RCovAS.d[ID][which(is.na(RCovAS.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RCovAS.d, RCovAS.sd, RCovAS.b.missing, RCovAS.lower)

RCovES.b <- Cov[, "RCOV_ES"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovES.d <- tapply(RCovES.b, ID, mean, na.rm = T) # Grid-level values
RCovES.b.missing <- is.na(RCovES.b) %>% as.integer # Index missing values to be imputed
RCovES.sd <- tapply(RCovES.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovES.sd==0))
  RCovES.sd[which(RCovES.sd==0)] <- min(RCovES.sd[which(RCovES.sd>0)])
RCovES.lower <- min(RCovES.b, na.rm = T) # Lower bound for imputing missing values
RCovES.b[which(is.na(RCovES.b))] <- RCovES.d[ID][which(is.na(RCovES.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RCovES.d, RCovES.sd, RCovES.b.missing, RCovES.lower)

RCovPine.b <- Cov[, "RCOV_Pine"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovPine.d <- tapply(RCovPine.b, ID, mean, na.rm = T) # Grid-level values
RCovPine.b.missing <- is.na(RCovPine.b) %>% as.integer # Index missing values to be imputed
RCovPine.sd <- tapply(RCovPine.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovPine.sd==0))
  RCovPine.sd[which(RCovPine.sd==0)] <- min(RCovPine.sd[which(RCovPine.sd>0)])
RCovPine.lower <- min(RCovPine.b, na.rm = T) # Lower bound for imputing missing values
RCovPine.b[which(is.na(RCovPine.b))] <- RCovPine.d[ID][which(is.na(RCovPine.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RCovPine.d, RCovPine.sd, RCovPine.b.missing, RCovPine.lower)

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.d <- tapply(shcov.b, ID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, ID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- shcov.d[ID[is.na(shcov.b)]]
rm(shcov.d, shcov.sd, shcov.b.missing)

RSC_Con.b <- Cov[, "RCShrb_UC"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RSC_Con.d <- tapply(RSC_Con.b, ID, mean, na.rm = T) # Grid-level values
RSC_Con.b.missing <- is.na(RSC_Con.b) %>% as.integer # Index missing values to be imputed
RSC_Con.sd <- tapply(RSC_Con.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RSC_Con.sd==0))
  RSC_Con.sd[which(RSC_Con.sd==0)] <- min(RSC_Con.sd[which(RSC_Con.sd>0)])
RSC_Con.lower <- min(RSC_Con.b, na.rm = T) # Lower bound for imputing missing values
RSC_Con.b[which(is.na(RSC_Con.b))] <- RSC_Con.d[ID][which(is.na(RSC_Con.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RSC_Con.d, RSC_Con.sd, RSC_Con.b.missing, RSC_Con.lower)

GHerb.b <- Cov[, "HerbCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GHerb.d <- tapply(GHerb.b, ID, mean, na.rm = T) # Grid-level values
GHerb.b.missing <- is.na(GHerb.b) %>% as.integer # Index missing values to be imputed
GHerb.sd <- tapply(GHerb.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GHerb.sd==0))
  GHerb.sd[which(GHerb.sd==0)] <- min(GHerb.sd[which(GHerb.sd>0)])
GHerb.lower <- min(GHerb.b, na.rm = T) # Lower bound for imputing missing values
GHerb.b[which(is.na(GHerb.b))] <- GHerb.d[ID][which(is.na(GHerb.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(GHerb.d, GHerb.sd, GHerb.b.missing, GHerb.lower)

Gwoody.b <- Cov[, "WoodyCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
Gwoody.d <- tapply(Gwoody.b, ID, mean, na.rm = T) # Grid-level values
Gwoody.b.missing <- is.na(Gwoody.b) %>% as.integer # Index missing values to be imputed
Gwoody.sd <- tapply(Gwoody.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(Gwoody.sd==0))
  Gwoody.sd[which(Gwoody.sd==0)] <- min(Gwoody.sd[which(Gwoody.sd>0)])
Gwoody.lower <- min(Gwoody.b, na.rm = T) # Lower bound for imputing missing values
Gwoody.b[which(is.na(Gwoody.b))] <- Gwoody.d[ID][which(is.na(Gwoody.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(Gwoody.d, Gwoody.sd, Gwoody.b.missing, Gwoody.lower)

GDD.b <- Cov[, "DDCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GDD.d <- tapply(GDD.b, ID, mean, na.rm = T) # Grid-level values
GDD.b.missing <- is.na(GDD.b) %>% as.integer # Index missing values to be imputed
GDD.sd <- tapply(GDD.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GDD.sd==0))
  GDD.sd[which(GDD.sd==0)] <- min(GDD.sd[which(GDD.sd>0)])
GDD.lower <- min(GDD.b, na.rm = T) # Lower bound for imputing missing values
GDD.b[which(is.na(GDD.b))] <- GDD.d[ID][which(is.na(GDD.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(GDD.d, GDD.sd, GDD.b.missing, GDD.lower)

# Tabulate partially observed richness #
SPR <- mod$sims.list$SR.point
dat.SR <- Cov %>% data.frame() %>%
  select(CanCov:DDCov) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

## Component plots ##
mns <- dat.SR %>%
  select(CanCov:DDCov) %>%
  as.matrix %>%
  apply(2, mean, na.rm = T)
sds <- dat.SR %>%
  select(CanCov:DDCov) %>%
  as.matrix %>%
  apply(2, sd, na.rm = T)
names(mns) <- names(sds) <- c("CanCov", "Aspen", "Spruce",
                              "Pine", "ShrubCov", "ConShrb",
                              "Herb", "Woody", "DeadDown")
omega <- mod$sims.list$omega
psi <- expit(mod$sims.list$d0)
b0 <- mod$sims.list$b0

vnam <- "Pine"
v <- RCovPine.b
b <- mod$sims.list$bb.RCovPine
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(RCOV_Pine)), aes(x = RCOV_Pine, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "ConShrb"
v <- RSC_Con.b
b <- mod$sims.list$bb.RSC_Con
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(RCOV_Pine)), aes(x = RCOV_Pine, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "Herb"
v <- GHerb.b
b <- mod$sims.list$bb.GHerb
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(HerbCov)), aes(x = HerbCov, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "Woody"
v <- Gwoody.b
b <- mod$sims.list$bb.Gwoody
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(WoodyCov)), aes(x = WoodyCov, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

#______Spruce-fir stratum_______#
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_habitat_reduced")

Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
ID <- Cov[, "gridIndex"] %>% as.factor

# Compile covariates #
ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.d <- tapply(ccov.b, ID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, ID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- ccov.d[ID[is.na(ccov.b)]]
rm(ccov.d, ccov.sd, ccov.b.missing)

RCovAS.b <- Cov[, "RCOV_AS"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovAS.d <- tapply(RCovAS.b, ID, mean, na.rm = T) # Grid-level values
RCovAS.b.missing <- is.na(RCovAS.b) %>% as.integer # Index missing values to be imputed
RCovAS.sd <- tapply(RCovAS.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovAS.sd==0))
  RCovAS.sd[which(RCovAS.sd==0)] <- min(RCovAS.sd[which(RCovAS.sd>0)])
RCovAS.lower <- min(RCovAS.b, na.rm = T) # Lower bound for imputing missing values
RCovAS.b[which(is.na(RCovAS.b))] <- RCovAS.d[ID][which(is.na(RCovAS.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RCovAS.d, RCovAS.sd, RCovAS.b.missing, RCovAS.lower)

RCovES.b <- Cov[, "RCOV_ES"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovES.d <- tapply(RCovES.b, ID, mean, na.rm = T) # Grid-level values
RCovES.b.missing <- is.na(RCovES.b) %>% as.integer # Index missing values to be imputed
RCovES.sd <- tapply(RCovES.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovES.sd==0))
  RCovES.sd[which(RCovES.sd==0)] <- min(RCovES.sd[which(RCovES.sd>0)])
RCovES.lower <- min(RCovES.b, na.rm = T) # Lower bound for imputing missing values
RCovES.b[which(is.na(RCovES.b))] <- RCovES.d[ID][which(is.na(RCovES.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RCovES.d, RCovES.sd, RCovES.b.missing, RCovES.lower)

RCovPine.b <- Cov[, "RCOV_Pine"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovPine.d <- tapply(RCovPine.b, ID, mean, na.rm = T) # Grid-level values
RCovPine.b.missing <- is.na(RCovPine.b) %>% as.integer # Index missing values to be imputed
RCovPine.sd <- tapply(RCovPine.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovPine.sd==0))
  RCovPine.sd[which(RCovPine.sd==0)] <- min(RCovPine.sd[which(RCovPine.sd>0)])
RCovPine.lower <- min(RCovPine.b, na.rm = T) # Lower bound for imputing missing values
RCovPine.b[which(is.na(RCovPine.b))] <- RCovPine.d[ID][which(is.na(RCovPine.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RCovPine.d, RCovPine.sd, RCovPine.b.missing, RCovPine.lower)

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.d <- tapply(shcov.b, ID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, ID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- shcov.d[ID[is.na(shcov.b)]]
rm(shcov.d, shcov.sd, shcov.b.missing)

RSC_Con.b <- Cov[, "RCShrb_UC"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RSC_Con.d <- tapply(RSC_Con.b, ID, mean, na.rm = T) # Grid-level values
RSC_Con.b.missing <- is.na(RSC_Con.b) %>% as.integer # Index missing values to be imputed
RSC_Con.sd <- tapply(RSC_Con.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RSC_Con.sd==0))
  RSC_Con.sd[which(RSC_Con.sd==0)] <- min(RSC_Con.sd[which(RSC_Con.sd>0)])
RSC_Con.lower <- min(RSC_Con.b, na.rm = T) # Lower bound for imputing missing values
RSC_Con.b[which(is.na(RSC_Con.b))] <- RSC_Con.d[ID][which(is.na(RSC_Con.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(RSC_Con.d, RSC_Con.sd, RSC_Con.b.missing, RSC_Con.lower)

GHerb.b <- Cov[, "HerbCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GHerb.d <- tapply(GHerb.b, ID, mean, na.rm = T) # Grid-level values
GHerb.b.missing <- is.na(GHerb.b) %>% as.integer # Index missing values to be imputed
GHerb.sd <- tapply(GHerb.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GHerb.sd==0))
  GHerb.sd[which(GHerb.sd==0)] <- min(GHerb.sd[which(GHerb.sd>0)])
GHerb.lower <- min(GHerb.b, na.rm = T) # Lower bound for imputing missing values
GHerb.b[which(is.na(GHerb.b))] <- GHerb.d[ID][which(is.na(GHerb.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(GHerb.d, GHerb.sd, GHerb.b.missing, GHerb.lower)

Gwoody.b <- Cov[, "WoodyCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
Gwoody.d <- tapply(Gwoody.b, ID, mean, na.rm = T) # Grid-level values
Gwoody.b.missing <- is.na(Gwoody.b) %>% as.integer # Index missing values to be imputed
Gwoody.sd <- tapply(Gwoody.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(Gwoody.sd==0))
  Gwoody.sd[which(Gwoody.sd==0)] <- min(Gwoody.sd[which(Gwoody.sd>0)])
Gwoody.lower <- min(Gwoody.b, na.rm = T) # Lower bound for imputing missing values
Gwoody.b[which(is.na(Gwoody.b))] <- Gwoody.d[ID][which(is.na(Gwoody.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(Gwoody.d, Gwoody.sd, Gwoody.b.missing, Gwoody.lower)

GDD.b <- Cov[, "DDCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GDD.d <- tapply(GDD.b, ID, mean, na.rm = T) # Grid-level values
GDD.b.missing <- is.na(GDD.b) %>% as.integer # Index missing values to be imputed
GDD.sd <- tapply(GDD.b, ID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GDD.sd==0))
  GDD.sd[which(GDD.sd==0)] <- min(GDD.sd[which(GDD.sd>0)])
GDD.lower <- min(GDD.b, na.rm = T) # Lower bound for imputing missing values
GDD.b[which(is.na(GDD.b))] <- GDD.d[ID][which(is.na(GDD.b))] # Insert means for imputing PctDead for points in non-outbreak grids
rm(GDD.d, GDD.sd, GDD.b.missing, GDD.lower)

# Tabulate partially observed richness #
SPR <- mod$sims.list$SR.point
dat.SR <- Cov %>% data.frame() %>%
  select(CanCov:DDCov) %>%
  mutate(Y = apply(SPR, 2, median) %>% as.numeric,
         Y.lo = apply(SPR, 2,function(x) quantile(x,prob=0.025,type=8)) %>%
           as.numeric,
         Y.hi = apply(SPR, 2, function(x) quantile(x,prob=0.975,type=8)) %>%
           as.numeric)

## Component plots ##
mns <- dat.SR %>%
  select(CanCov:DDCov) %>%
  as.matrix %>%
  apply(2, mean, na.rm = T)
sds <- dat.SR %>%
  select(CanCov:DDCov) %>%
  as.matrix %>%
  apply(2, sd, na.rm = T)
names(mns) <- names(sds) <- c("CanCov", "Aspen", "Spruce",
                              "Pine", "ShrubCov", "ConShrb",
                              "Herb", "Woody", "DeadDown")
omega <- mod$sims.list$omega
psi <- expit(mod$sims.list$d0)
b0 <- mod$sims.list$b0

vnam <- "Pine"
v <- RCovPine.b
b <- mod$sims.list$bb.RCovPine
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(RCOV_Pine)), aes(x = RCOV_Pine, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "ShrubCov"
v <- shcov.b
b <- mod$sims.list$bb.ShCov
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(shrub_cover)), aes(x = shrub_cover, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "ConShrb"
v <- RSC_Con.b
b <- mod$sims.list$bb.RSC_Con
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(RCShrb_UC)), aes(x = RCShrb_UC, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "Herb"
v <- GHerb.b
b <- mod$sims.list$bb.GHerb
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(HerbCov)), aes(x = HerbCov, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

vnam <- "DeadDown"
v <- GDD.b
b <- mod$sims.list$bb.GDD
dat.pred <- data.frame(z = seq(min(v, na.rm = T), max(v, na.rm = T), length.out = 10)) %>%
  mutate(x = z * sds[vnam] + mns[vnam])
SR.pred <- matrix(NA, nrow = dim(b0)[1], ncol = nrow(dat.pred))
for(i in 1:dim(SR.pred)[2]) {
  theta <- expit(b0 + b*dat.pred$z[i])
  SR.pred[, i] <- omega * apply(psi * theta, 1, sum)
}
dat.pred <- dat.pred %>%
  mutate(pred.md = apply(SR.pred, 2, median),
         pred.lo = apply(SR.pred, 2, function(x) quantile(x, prob = 0.05, type = 8)),
         pred.hi = apply(SR.pred, 2, function(x) quantile(x, prob = 0.95, type = 8)))
p <- ggplot(data = dat.SR %>% filter(!is.na(DDCov)), aes(x = DDCov, y = Y)) + 
  geom_point(alpha = 0.1) + 
  geom_errorbar(aes(ymin = Y.lo, ymax = Y.hi), width = 0, alpha = 0.1) +
  geom_ribbon(data = dat.pred, aes(x = x, ymin = pred.lo, ymax = pred.hi),
              alpha = 0.3, inherit.aes = F) +
  geom_line(data = dat.pred, aes(x = x, y = pred.md),
            size = 1.5, inherit.aes = F) +
  labs(x = vnam, y = NULL) +
  ylim(miny, maxy)
assign(str_c("p.", vnam, ".", stratum), p)

#_____ Put them all together _____#
p <- ggdraw() + 
  draw_plot(p.Pine.LP, x = 0.03, y = 0.5, width = 0.2425, height = 0.45) +
  draw_plot(p.ConShrb.LP, x = 0.2725, y = 0.5, width = 0.2425, height = 0.45) +
  draw_plot(p.Herb.LP, x = 0.5150, y = 0.5, width = 0.2425, height = 0.45) +
  draw_plot(p.Woody.LP, x = 0.7575, y = 0.5, width = 0.2425, height = 0.45) +
  draw_plot(p.Pine.SF, x = 0.03, y = 0, width = 0.194, height = 0.45) +
  draw_plot(p.ShrubCov.SF, x = 0.224, y = 0, width = 0.194, height = 0.45) +
  draw_plot(p.ConShrb.SF, x = 0.418, y = 0, width = 0.194, height = 0.45) +
  draw_plot(p.Herb.SF, x = 0.612, y = 0, width = 0.194, height = 0.45) +
  draw_plot(p.DeadDown.SF, x = 0.806, y = 0, width = 0.194, height = 0.45) +
  draw_plot_label(c("Bird species richness", "Lodgepole", "Spruce-fir"),
                  x = c(0, 0.48, 0.48), y = c(0.37, 0.98, 0.49), size = c(20, 20, 20),
                  angle = c(90, 0, 0), hjust = c(0, 0, 0))

save_plot("Plot_richness_habitat.tiff", p, ncol = 4, nrow = 2.5, dpi = 200)
