library(jagsUI)
library(stringr)
library(dplyr)

setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
load("Data_compiled.RData")

## Lodgepole pine stratum ##
# Detection data #
Y <- Y.LP
TPeriod <- T.LP
gridID <- Cov.LP[, "gridIndex"]
n.grid <- max(gridID)
n.point <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Covariates #
PctDead.b <- Cov.LP[, "DeadConif"] # Point-level values
PctDead.b[which(Cov.LP[, "YSO"] > 12)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
PctDead.d <- tapply(PctDead.b, gridID, mean, na.rm = T) # Grid-level values
PctDead.b[which(is.na(PctDead.b) & is.na(Cov.LP[, "YSO"]))] <- # Insert means for imputing PctDead outside ADS polygons
  mean(PctDead.b[which(!is.na(PctDead.b) & is.na(Cov.LP[, "YSO"]))])
PctDead.b[which(is.na(PctDead.b) & !is.na(Cov.LP[, "YSO"]))] <- # Insert means for imputing PctDead inside ADS polygons
  mean(PctDead.b[which(!is.na(PctDead.b) & !is.na(Cov.LP[, "YSO"]))])
PctDead.sd <- sd(PctDead.b[which(is.na(Cov.LP[, "YSO"]))]) %>% rep(length(PctDead.b)) # SDs for imputing missing values outside ADS polygons
PctDead.sd[which(!is.na(Cov.LP[, "YSO"]))] <- sd(PctDead.b[which(!is.na(Cov.LP[, "YSO"]))]) # SDs for imputing missing values inside ADS polygons
PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values
PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed

YSO.b <- Cov.LP[, "YSO"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YSO.d <- tapply(YSO.b, gridID, mean, na.rm = T) # Grid-level values
YSO.b[is.na(YSO.b)] <- 0
YSO.d[is.na(YSO.d)] <- 0

Outbrk.b <- Cov.LP[, "YSO"] %>% (function(x) as.integer(!is.na(x)))
Outbrk.d <- tapply(Outbrk.b, gridID, max)

TWIP.b <- Cov.LP[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
TWIP.d <- tapply(TWIP.b, gridID, mean, na.rm = T) # Grid-level values
TWIP.sd <- tapply(TWIP.b, gridID, sd, na.rm = T) # Grid-specific SDs for imputing missing values
TWIP.b.missing <- is.na(TWIP.b) %>% as.integer # Index missing values to be imputed
TWIP.b[is.na(TWIP.b)] <- 0

DOY.b <- Cov.LP[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov.LP[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
Time.b[point.list.LP == "LP-057-06"] <- 
  mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
# Imputed missing time value based on survey times for point likely completed just before and just after.
  # Inferred which points after viewing the other points and their associated times in ArcMap.

ccov.b <- Cov.LP[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.means <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- 0

shcov.b <- Cov.LP[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.means <- tapply(shcov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- 0

data <- list("Y", "TPeriod", "gridID", "n.grid", "n.point", "n.spp", "PctDead.b",
             "PctDead.d", "PctDead.sd", "PctDead.lower", "PctDead.b.missing",
             "YSO.b", "YSO.d", "Outbrk.b", "Outbrk.d",
             "TWIP.b", "TWIP.d", "TWIP.sd", "TWIP.b.missing",
             "DOY.b", "Time.b",
             "ccov.b", "ccov.means", "ccov.sd", "ccov.b.missing",
             "shcov.b", "shcov.means", "shcov.sd", "shcov.b.missing")

# Assemble the initial values for JAGS.
u.init <- Y
dimnames(u.init) <- NULL
z.init <- matrix(NA, nrow = n.grid, ncol = n.spp)
for(sp in 1:n.spp) {
  z.init[, sp] <- tapply(u.init[, sp], gridID, max)
}
w.init <- apply(Y, 2, max) %>% as.integer

inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1),
       tvar.sigma.d0 = rnorm(1), tvar.Betad.PctDead = rnorm(1), tvar.Betad.PctDead2 = rnorm(1),
       tvar.Betad.Outbrk = rnorm(1), tvar.Betad.YSO = rnorm(1), tvar.Betad.YSO2 = rnorm(1),
       tvar.Betad.PctDdXYSO = rnorm(1), tvar.Betad.TWIP = rnorm(1), tvar.Betab.PctDead = rnorm(1),
       tvar.Betab.PctDead2 = rnorm(1), tvar.Betab.Outbrk = rnorm(1), tvar.Betab.YSO = rnorm(1),
       tvar.Betab.YSO2 = rnorm(1), tvar.Betab.PctDdXYSO = rnorm(1), tvar.Betab.TWIP = rnorm(1),
       tvar.Betaa.Time = rnorm(1), tvar.Betaa.Time2 = rnorm(1), tvar.Betaa.DOY = rnorm(1),
       tvar.Betaa.DOY2 = rnorm(1), tvar.Betaa.CCov = rnorm(1), tvar.Betaa.SHCov = rnorm(1))

# Assemble the parameters vector for JAGS (What we want to track).
parameters <- c("d0", "b0", "a0",
                "bd.pdead", "bd.pdead2", "bd.outbrk", "bd.YSO", "bd.YSO2", "bd.pdXYSO", "bd.TWIP",
                "bb.pdead", "bb.pdead2", "bb.outbrk", "bb.YSO", "bb.YSO2", "bb.pdXYSO", "bb.TWIP",
                "ba.Time", "ba.Time2", "ba.DOY", "ba.DOY2", "ba.ccov", "ba.shcov",
                "omega", "rho.ab", "rho.bd")

# MCMC values.  Change to what works well for you.
nc <- 6
nb <- 1000
ni <- 90000
nt <- 100

st.time <- Sys.time()
out <- jagsUI(data, inits, parameters.to.save = parameters, "jags_outbreak.txt", n.thin=nt, n.chains=nc,
              n.burnin=nb, n.iter=ni, parallel=TRUE)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)
