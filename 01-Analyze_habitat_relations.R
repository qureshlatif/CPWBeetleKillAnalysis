library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#### Script inputs ####
stratum <- "LP"
maxYSOForPD <- 12 # Set to 12 for LP and 9 for SF
model.file <- "CPWBeetleKillAnalysis/model_habitat_LP_global.jags"
RESQ.model <- "mod_RESQ_outbreak_HZdist_LP_reduced"

# Data objects to send to JAGS
data <- list("Y", "TPeriod", "gridID", "n.grid", "n.point", "n.spp",
             "PctDead.d", "PctDead.b", "PctDead.sd", "PctDead.lower", "PctDead.b.missing",
             "TWIP.d", "TWI.d", "heatload.d", "RDens.d", "WILD.d",
             "RCovAS.d", "RCovAS.b", "RCovAS.sd", "RCovAS.b.missing", "RCovAS.lower",
             "RCovES.d", "RCovES.b", "RCovES.sd", "RCovES.b.missing", "RCovES.lower",
             "RCovPine.d", "RCovPine.b", "RCovPine.sd", "RCovPine.b.missing", "RCovPine.lower",
             "ccov.b", "ccov.d", "ccov.sd", "ccov.b.missing",
             "shcov.b", "shcov.d", "shcov.sd", "shcov.b.missing",
             "RSC_Con.d", "RSC_Con.b", "RSC_Con.sd", "RSC_Con.b.missing", "RSC_Con.lower",
             "GHerb.d", "GHerb.b", "GHerb.sd", "GHerb.b.missing", "GHerb.lower",
             "Gwoody.d", "Gwoody.b", "Gwoody.sd", "Gwoody.b.missing", "Gwoody.lower",
             "GDD.d", "GDD.b", "GDD.sd", "GDD.b.missing", "GDD.lower",
             "RESQ.d", "RESQ.b", "RESQ.wts", "RESQ.b.simp", 
             "DOY.b", "Time.b")

# Stuff to save from JAGS
parameters <- c("omega", "rho.ab", "rho.bd",
                "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
                
                "Betad.PctDead", "sigma.Betad.PctDead", "Betad.YSO", "sigma.Betad.YSO",
                "Betad.TWIP", "sigma.Betad.TWIP", "Betad.RDens", "sigma.Betad.RDens", "Betad.WILD", "sigma.Betad.WILD",
                "Betad.RCovAS", "sigma.Betad.RCovAS", "Betad.RCovES", "sigma.Betad.RCovES", "Betad.RCovPine", "sigma.Betad.RCovPine",
                "Betad.CanCov", "sigma.Betad.CanCov", "Betad.ShCov", "sigma.Betad.ShCov", "Betad.RSC_Con", "sigma.Betad.RSC_Con",
                "Betad.GHerb", "sigma.Betad.GHerb", "Betad.Gwoody", "sigma.Betad.Gwoody", "Betad.GDD", "sigma.Betad.GDD",
                "Betad.RESQ", "sigma.Betad.RESQ",
                
                "Betab.PctDead", "sigma.Betab.PctDead",
                "Betab.RCovAS", "sigma.Betab.RCovAS", "Betab.RCovES", "sigma.Betab.RCovES", "Betab.RCovPine", "sigma.Betab.RCovPine",
                "Betab.CanCov", "sigma.Betab.CanCov", "Betab.ShCov", "sigma.Betab.ShCov", "Betab.RSC_Con", "sigma.Betab.RSC_Con",
                "Betab.GHerb", "sigma.Betab.GHerb", "Betab.Gwoody", "sigma.Betab.Gwoody", "Betab.GDD", "sigma.Betab.GDD",
                "Betab.RESQ", "sigma.Betab.RESQ",
                
                "Betaa.Time", "sigma.Betaa.Time", "Betaa.Time2", "sigma.Betaa.Time2",
                "Betaa.DOY", "sigma.Betaa.DOY", "Betaa.DOY2", "sigma.Betaa.DOY2",
                "Betaa.CCov", "sigma.Betaa.CCov",
                "Betaa.SHCov", "sigma.Betaa.SHCov",
                
                "d0", "b0", "a0",
                "bd.pdead", "bd.TWIP", "bd.RDens", "bd.WILD",
                "bd.RCovAS", "bd.RCovES", "bd.RCovPine",
                "bd.CanCov", "bd.ShCov", "bd.RSC_Con",
                "bd.GHerb", "bd.Gwoody", "bd.GDD", "bd.RESQ",
                
                "bb.RCovAS", "bb.RCovES", "bb.RCovPine",
                "bb.pdead", "bb.CanCov", "bb.ShCov", "bb.RSC_Con",
                "bb.GHerb", "bb.Gwoody", "bb.GDD", "bb.RESQ",
                
                "ba.Time", "ba.Time2", "ba.DOY", "ba.DOY2", "ba.ccov", "ba.shcov",
                
                "SR.grid", "SR.point")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1), tvar.sigma.d0 = rnorm(1),
       tvar.Betad.PctDead = rnorm(1), tvar.Betad.TWIP = rnorm(1), tvar.Betad.TWI = rnorm(1), tvar.Betad.heatload = rnorm(1),
       tvar.Betad.RDens = rnorm(1), tvar.Betad.WILD = rnorm(1),
       tvar.Betad.RCovAS = rnorm(1), tvar.Betad.RCovES = rnorm(1), tvar.Betad.RCovPine = rnorm(1),
       tvar.Betad.CanCov = rnorm(1), tvar.Betad.ShCov = rnorm(1), tvar.Betad.RSC_Con = rnorm(1),
       tvar.Betad.GHerb = rnorm(1), tvar.Betad.Gwoody = rnorm(1), tvar.Betad.GDD = rnorm(1), tvar.Betad.RESQ = rnorm(1),
       
       tvar.Betab.PctDead = rnorm(1), tvar.Betab.RCovAS = rnorm(1), tvar.Betab.RCovES = rnorm(1), tvar.Betab.RCovPine = rnorm(1),
       tvar.Betab.CanCov = rnorm(1), tvar.Betab.ShCov = rnorm(1), tvar.Betab.RSC_Con = rnorm(1),
       tvar.Betab.GHerb = rnorm(1), tvar.Betab.Gwoody = rnorm(1), tvar.Betab.GDD = rnorm(1), tvar.Betab.RESQ = rnorm(1),
       
       tvar.Betaa.Time = rnorm(1), tvar.Betaa.Time2 = rnorm(1),
       tvar.Betaa.DOY = rnorm(1), tvar.Betaa.DOY2 = rnorm(1),
       tvar.Betaa.CCov = rnorm(1), tvar.Betaa.SHCov = rnorm(1))

# MCMC values
nc <- 2 #3 # number of chains
nb <- 5000 #10000 # burn in
ni <- 10000 #100000 # number of iterations
nt <- 1 #10 # thinning

save.out <- "mod_LPcommunity_habitat"
##########################

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

#TWIP.b <- Cov[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
TWIP.d <- Cov[, "TWIP"]  %>% # Grid-level values only
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

TWI.d <- Cov[, "TWI"]  %>% # Grid-level values only
  tapply(gridID, mean) %>%
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))

heatload.d <- Cov[, "heatload"]  %>% # Grid-level values only
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

ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.d <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- ccov.d[gridID[is.na(ccov.b)]]

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.d <- tapply(shcov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- shcov.d[gridID[is.na(shcov.b)]]

RSC_Con.b <- Cov[, "RCShrb_UC"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
RSC_Con.d <- tapply(RSC_Con.b, gridID, mean, na.rm = T) # Grid-level values
RSC_Con.b.missing <- is.na(RSC_Con.b) %>% as.integer # Index missing values to be imputed
RSC_Con.sd <- tapply(RSC_Con.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(RSC_Con.sd==0))
  RSC_Con.sd[which(RSC_Con.sd==0)] <- min(RSC_Con.sd[which(RSC_Con.sd>0)])
RSC_Con.lower <- min(RSC_Con.b, na.rm = T) # Lower bound for imputing missing values
RSC_Con.b[which(is.na(RSC_Con.b))] <- RSC_Con.d[gridID][which(is.na(RSC_Con.b))] # Insert means for imputing PctDead for points in non-outbreak grids

GHerb.b <- Cov[, "HerbCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GHerb.d <- tapply(GHerb.b, gridID, mean, na.rm = T) # Grid-level values
GHerb.b.missing <- is.na(GHerb.b) %>% as.integer # Index missing values to be imputed
GHerb.sd <- tapply(GHerb.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GHerb.sd==0))
  GHerb.sd[which(GHerb.sd==0)] <- min(GHerb.sd[which(GHerb.sd>0)])
GHerb.lower <- min(GHerb.b, na.rm = T) # Lower bound for imputing missing values
GHerb.b[which(is.na(GHerb.b))] <- GHerb.d[gridID][which(is.na(GHerb.b))] # Insert means for imputing PctDead for points in non-outbreak grids

Gwoody.b <- Cov[, "WoodyCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
Gwoody.d <- tapply(Gwoody.b, gridID, mean, na.rm = T) # Grid-level values
Gwoody.b.missing <- is.na(Gwoody.b) %>% as.integer # Index missing values to be imputed
Gwoody.sd <- tapply(Gwoody.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(Gwoody.sd==0))
  Gwoody.sd[which(Gwoody.sd==0)] <- min(Gwoody.sd[which(Gwoody.sd>0)])
Gwoody.lower <- min(Gwoody.b, na.rm = T) # Lower bound for imputing missing values
Gwoody.b[which(is.na(Gwoody.b))] <- Gwoody.d[gridID][which(is.na(Gwoody.b))] # Insert means for imputing PctDead for points in non-outbreak grids

GDD.b <- Cov[, "DDCov"] %>% # Point-level values
  (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
GDD.d <- tapply(GDD.b, gridID, mean, na.rm = T) # Grid-level values
GDD.b.missing <- is.na(GDD.b) %>% as.integer # Index missing values to be imputed
GDD.sd <- tapply(GDD.b, gridID, sd, na.rm = T) # SDs for imputing missing values  for non-outbreak grids
if(any(GDD.sd==0))
  GDD.sd[which(GDD.sd==0)] <- min(GDD.sd[which(GDD.sd>0)])
GDD.lower <- min(GDD.b, na.rm = T) # Lower bound for imputing missing values
GDD.b[which(is.na(GDD.b))] <- GDD.d[gridID][which(is.na(GDD.b))] # Insert means for imputing PctDead for points in non-outbreak grids

library(R.utils)
mod <- loadObject(RESQ.model)
RESQ.b <- mod$sims.list$N
#cl.size.samps <- mod$sims.list$cl.size
#RESQ.b <- apply(N, c(1, 2), function(x) sum(base::sample(cl.size.samps, x)))
rm(mod) #N, cl.size.samps
RESQ.d <- matrix(NA, nrow = dim(RESQ.b)[1], ncol = n.grid)
for(i in 1:nrow(RESQ.d)) RESQ.d[i, ] <- tapply(RESQ.b[i, ], gridID, sum) # Probably can make this faster with series of applys or mapply or something
# z-score
RESQ.b <- RESQ.b %>%
  (function(x) (x - mean(x)) / sd(x))
RESQ.d <- RESQ.d %>%
  (function(x) (x - mean(x)) / sd(x))
RESQ.wts <- rep(1, dim(RESQ.b)[1])
RESQ.b.simp <- apply(RESQ.b, 2, mean)

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
