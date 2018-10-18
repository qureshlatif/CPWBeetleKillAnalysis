library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#________ Script inputs________#
stratum <- "LP"
maxYSOForPD <- 12 # Set to 12 for LP and 9 for SF
model.file <- "CPWBeetleKillAnalysis/model_outbreak_LP_reduced.jags"

# MCMC values
nc <- 3 # number of chains
nb <- 10000 # burn in
ni <- 100000 # number of iterations
nt <- 10 # thinning

save.out <- "mod_SFcommunity_outbreak"
#______________________________#


# Data objects to send to JAGS
data <- list("Y", "TPeriod", "gridID", "n.grid", "n.point", "n.spp", "PctDead.b",
             "PctDead.d", "PctDead.sd", "PctDead.lower", "PctDead.b.missing",
             "YSO.b", "YSO.mins", "YSO.maxs", "YSO.missing", "YSO.d",
             "heatload.d", "TWI.d", "TWIP.d", "RDens.d", "WILD.d",
             "RCovAS.d", "RCovAS.b", "RCovAS.sd", "RCovAS.b.missing", "RCovAS.lower",
             "RCovES.d", "RCovES.b", "RCovES.sd", "RCovES.b.missing", "RCovES.lower",
             "RCovPine.d", "RCovPine.b", "RCovPine.sd", "RCovPine.b.missing", "RCovPine.lower",
             "DOY.b", "Time.b")

# Stuff to save from JAGS
parameters <- c("omega", "rho.ab", "rho.bd",
                "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
                
                "Betad.PctDead", "sigma.Betad.PctDead", "Betad.YSO", "sigma.Betad.YSO",
                "Betad.heatload", "sigma.Betad.heatload", "Betad.TWI", "sigma.Betad.TWI", "Betad.TWIP", "sigma.Betad.TWIP",
                "Betad.RDens", "sigma.Betad.RDens", "Betad.WILD", "sigma.Betad.WILD",
                "Betad.RCovAS", "sigma.Betad.RCovAS", "Betad.RCovES", "sigma.Betad.RCovES", "Betad.RCovPine", "sigma.Betad.RCovPine",
                "Betab.PctDead", "sigma.Betab.PctDead",
                "Betab.YSO", "sigma.Betab.YSO", "Betab.YSO2", "sigma.Betab.YSO2", "Betab.PctDdXYSO", "sigma.Betab.PctDdXYSO",
                "Betab.RCovAS", "sigma.Betab.RCovAS", "Betab.RCovES", "sigma.Betab.RCovES", "Betab.RCovPine", "sigma.Betab.RCovPine",
                "Betaa.Time", "sigma.Betaa.Time", "Betaa.Time2", "sigma.Betaa.Time2",
                "Betaa.DOY", "sigma.Betaa.DOY", "Betaa.DOY2", "sigma.Betaa.DOY2",
                "Betaa.PctDead", "sigma.Betaa.PctDead",
                "Betaa.YSO", "sigma.Betaa.YSO", "Betaa.YSO2", "sigma.Betaa.YSO2", "Betaa.PctDdXYSO", "sigma.Betaa.PctDdXYSO",
                
                "d0", "b0", "a0",
                "bd.pdead", "bd.YSO", "bd.heatload", "bd.TWI", "bd.TWIP",
                "bd.RDens", "bd.WILD",
                "bd.RCovAS", "bd.RCovES", "bd.RCovPine",
                "bb.RCovAS", "bb.RCovES", "bb.RCovPine",
                "bb.pdead", "bb.YSO", "bb.YSO2", "bb.pdXYSO",
                "bb.RCovAS", "bb.RCovES", "bb.RCovPine",
                "ba.Time", "ba.Time2", "ba.DOY", "ba.DOY2",
                "ba.pdead", "ba.YSO", "ba.YSO2", "ba.pdXYSO",
                
                "SR.grid", "SR.point")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1), tvar.sigma.d0 = rnorm(1),
       tvar.Betad.PctDead = rnorm(1), tvar.Betad.YSO = rnorm(1),
       tvar.Betad.heatload = rnorm(1), tvar.Betad.TWI = rnorm(1), tvar.Betad.TWIP = rnorm(1),
       tvar.Betad.RDens = rnorm(1), tvar.Betad.WILD = rnorm(1),
       tvar.Betad.RCovAS = rnorm(1), tvar.Betad.RCovES = rnorm(1), tvar.Betad.RCovPine = rnorm(1),
       tvar.Betab.PctDead = rnorm(1), tvar.Betab.YSO = rnorm(1),
       tvar.Betab.YSO2 = rnorm(1), tvar.Betab.PctDdXYSO = rnorm(1),
       tvar.Betab.RCovAS = rnorm(1), tvar.Betab.RCovES = rnorm(1), tvar.Betab.RCovPine = rnorm(1),
       tvar.Betaa.Time = rnorm(1), tvar.Betaa.Time2 = rnorm(1),
       tvar.Betaa.DOY = rnorm(1), tvar.Betaa.DOY2 = rnorm(1),
       tvar.Betaa.PctDead = rnorm(1), tvar.Betaa.YSO = rnorm(1),
       tvar.Betaa.YSO2 = rnorm(1), tvar.Betaa.PctDdXYSO = rnorm(1))

# Detection data #
Y <- str_c("Y.", stratum) %>% as.name %>% eval
TPeriod <- str_c("T.", stratum) %>% as.name %>% eval
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
gridID <- Cov[, "gridIndex"]
n.grid <- max(gridID)
n.point <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Covariates #
outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks

PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(Cov[, "YSO"] > maxYSOForPD)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
PctDead.d <- tapply(PctDead.b, gridID, mean, na.rm = T) # Grid-level values
PctDead.b[which(is.na(PctDead.b) & !outbreak_grids[gridID])] <- # Insert means for imputing PctDead for points in non-outbreak grids
  mean(PctDead.b[which(!is.na(PctDead.b) & !outbreak_grids[gridID])])
PctDead.sd <- sd(PctDead.b[which(!outbreak_grids[gridID])]) %>% rep(length(PctDead.b)) # SDs for imputing missing values  for non-outbreak grids
PctDead.b[which(is.na(PctDead.b) & outbreak_grids[gridID])] <- # Insert means for imputing PctDead for points within outbreak grids
  mean(PctDead.b[which(!is.na(PctDead.b) & outbreak_grids[gridID])])
PctDead.sd[which(outbreak_grids[gridID])] <- sd(PctDead.b[which(outbreak_grids[gridID])]) # SDs for imputing missing values inside ADS polygons
PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values
PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed

YSO.b <- Cov[, "YSO"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YSO.d <- tapply(YSO.b, gridID, mean, na.rm = T) # Grid-level values
YSO.mins <- tapply(YSO.b, gridID, min, na.rm = T) %>% as.numeric # Grid-level values
YSO.maxs <- tapply(YSO.b, gridID, max, na.rm = T) %>% as.numeric  # Grid-level values
YSO.mins[which(YSO.mins == Inf)] <- -1
YSO.maxs[which(YSO.maxs == -Inf)] <- 1
ind <- which(YSO.mins >= YSO.maxs)
YSO.mins[ind] <- YSO.mins[ind] - 0.01
YSO.maxs[ind] <- YSO.maxs[ind] + 0.01
rm(ind)
YSO.missing <- (is.na(YSO.b) & outbreak_grids[gridID]) %>% as.integer
YSO.b[is.na(YSO.b)] <- 0
YSO.d[is.na(YSO.d)] <- 0

#TWIP.b <- Cov[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
heatload.d <- Cov[, "heatload"]  %>% # Grid-level values only
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWI.d <- Cov[, "TWI"]  %>% # Grid-level values only
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWIP.d <- Cov[, "TWIP"]  %>% # Grid-level values only
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

RDens.d <- Cov[, "Rd_dens1km"]  %>% # Grid-level values only
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

WILD.d <- Cov[, "WILD"]  %>% # Grid-level values only
  tapply(gridID, function(x) (mean(x) >= 0.5)*1)

RCovAS.b <- Cov[, "RCOV_AS"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovAS.d <- tapply(RCovAS.b, gridID, mean, na.rm = T) # Grid-level values
RCovAS.b.missing <- is.na(RCovAS.b) %>% as.integer # Index missing values to be imputed
RCovAS.sd <- tapply(RCovAS.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovAS.sd==0))
  RCovAS.sd[which(RCovAS.sd==0)] <- min(RCovAS.sd[which(RCovAS.sd>0)])
RCovAS.lower <- min(RCovAS.b, na.rm = T) # Lower bound for imputing missing values
RCovAS.b[which(is.na(RCovAS.b))] <- RCovAS.d[gridID][which(is.na(RCovAS.b))] # Insert means for imputing PctDead for points in non-outbreak grids

RCovES.b <- Cov[, "RCOV_ES"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovES.d <- tapply(RCovES.b, gridID, mean, na.rm = T) # Grid-level values
RCovES.b.missing <- is.na(RCovES.b) %>% as.integer # Index missing values to be imputed
RCovES.sd <- tapply(RCovES.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovES.sd==0))
  RCovES.sd[which(RCovES.sd==0)] <- min(RCovES.sd[which(RCovES.sd>0)])
RCovES.lower <- min(RCovES.b, na.rm = T) # Lower bound for imputing missing values
RCovES.b[which(is.na(RCovES.b))] <- RCovES.d[gridID][which(is.na(RCovES.b))] # Insert means for imputing PctDead for points in non-outbreak grids

RCovPine.b <- Cov[, "RCOV_Pine"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RCovPine.d <- tapply(RCovPine.b, gridID, mean, na.rm = T) # Grid-level values
RCovPine.b.missing <- is.na(RCovPine.b) %>% as.integer # Index missing values to be imputed
RCovPine.sd <- tapply(RCovPine.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RCovPine.sd==0))
  RCovPine.sd[which(RCovPine.sd==0)] <- min(RCovPine.sd[which(RCovPine.sd>0)])
RCovPine.lower <- min(RCovPine.b, na.rm = T) # Lower bound for imputing missing values
RCovPine.b[which(is.na(RCovPine.b))] <- RCovPine.d[gridID][which(is.na(RCovPine.b))] # Insert means for imputing PctDead for points in non-outbreak grids

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
if(stratum == "LP")
  Time.b[point.list.LP == "LP-057-06"] <-
    mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
      # Impute missing time value in LP data based on survey times for point likely completed just before and just after.
        # Inferred which points after viewing the other points and their associated times in ArcMap.

# Assemble the initial values for JAGS.
u.init <- Y
dimnames(u.init) <- NULL
z.init <- matrix(NA, nrow = n.grid, ncol = n.spp)
for(sp in 1:n.spp) {
  z.init[, sp] <- tapply(u.init[, sp], gridID, max)
}
w.init <- apply(Y, 2, max) %>% as.integer

# Fit model
st.time <- Sys.time()
out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc,
              n.burnin=nb, n.iter=ni, parallel=TRUE)
end.time <- Sys.time()
run.time <- end.time - st.time
run.time
rm(st.time,end.time)

# Check the basics
max(out$summary[,"Rhat"])
#sort(out$summary[,"Rhat"], decreasing = T)[1:100]

min(out$summary[,"n.eff"])
#sort(out$summary[,"n.eff"])[1:50]

# Save output
library(R.utils)
saveObject(out, save.out)
