library(jagsUI)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#### Script inputs ####
stratum <- "LP"
model.file <- "model_outbreak_2ndOrderAtPnt.jags"

# Data objects to send to JAGS
data <- list("Y", "TPeriod", "gridID", "n.grid", "n.point", "n.spp", "PctDead.b",
             "PctDead.d", "PctDead.sd", "PctDead.lower", "PctDead.b.missing",
             "YSO.b", "YSO.d", "Outbrk.b", "Outbrk.d",
             "TWIP.b", "TWIP.d", "TWIP.sd", "TWIP.b.missing",
             "DOY.b", "Time.b",
             "ccov.b", "ccov.means", "ccov.sd", "ccov.b.missing",
             "shcov.b", "shcov.means", "shcov.sd", "shcov.b.missing")

# Stuff to save from JAGS
parameters <- c("omega", "rho.ab", "rho.bd",
                "alpha0", "sigma.a0", "beta0", "sigma.b0", "delta0", "sigma.d0",
                
                "Betad.PctDead", "sigma.Betad.PctDead", "Betad.Outbrk", "sigma.Betad.Outbrk", "Betad.YSO", "sigma.Betad.YSO",
                "Betad.TWIP", "sigma.Betad.TWIP",
                "Betab.PctDead", "sigma.Betab.PctDead", "Betab.PctDead2", "sigma.Betab.PctDead2", "Betab.Outbrk", "sigma.Betab.Outbrk",
                "Betab.YSO", "sigma.Betab.YSO", "Betab.YSO2", "sigma.Betab.YSO2", "Betab.PctDdXYSO", "sigma.Betab.PctDdXYSO",
                "Betab.TWIP", "sigma.Betab.TWIP",
                "Betaa.Time", "sigma.Betaa.Time", "Betaa.DOY", "sigma.Betaa.DOY", "Betaa.CCov", "sigma.Betaa.CCov",
                "Betaa.SHCov", "sigma.Betaa.SHCov",
                
                "d0", "b0", "a0",
                "bd.pdead", "bd.outbrk", "bd.YSO", "bd.TWIP",
                "bb.pdead", "bb.outbrk", "bb.YSO", "bb.TWIP",
                "ba.Time", "ba.DOY", "ba.ccov", "ba.shcov")

# Function for setting initial values in JAGS
inits <- function()
  list(z=z.init, u=u.init, w=w.init, tvar.sigma.a0 = rnorm(1), tvar.sigma.b0 = rnorm(1), tvar.sigma.d0 = rnorm(1),
       tvar.Betad.PctDead = rnorm(1), tvar.Betad.Outbrk = rnorm(1), tvar.Betad.YSO = rnorm(1), tvar.Betad.TWIP = rnorm(1),
       tvar.Betab.PctDead = rnorm(1), tvar.Betab.PctDead2 = rnorm(1), tvar.Betab.Outbrk = rnorm(1), tvar.Betab.YSO = rnorm(1),
       tvar.Betab.YSO2 = rnorm(1), tvar.Betab.PctDdXYSO = rnorm(1), tvar.Betab.TWIP = rnorm(1), tvar.Betaa.Time = rnorm(1),
       tvar.Betaa.DOY = rnorm(1), tvar.Betaa.CCov = rnorm(1), tvar.Betaa.SHCov = rnorm(1))

# MCMC values
nc <- 3 # number of chains
nb <- 1000 # burn in
ni <- 15000 # number of iterations
nt <- 10 # thinning

save.out <- "mod_LPcommunity_no2ndorder"
##########################

## Lodgepole pine stratum ##
# Detection data #
Y <- str_c("Y.", stratum) %>% as.name %>% eval
TPeriod <- str_c("T.", stratum) %>% as.name %>% eval
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval
gridID <- Cov[, "gridIndex"]
n.grid <- max(gridID)
n.point <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Covariates #
PctDead.b <- Cov[, "DeadConif"] # Point-level values
PctDead.b[which(Cov[, "YSO"] > 12)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
PctDead.d <- tapply(PctDead.b, gridID, mean, na.rm = T) # Grid-level values
PctDead.b[which(is.na(PctDead.b) & is.na(Cov[, "YSO"]))] <- # Insert means for imputing PctDead outside ADS polygons
  mean(PctDead.b[which(!is.na(PctDead.b) & is.na(Cov[, "YSO"]))])
PctDead.b[which(is.na(PctDead.b) & !is.na(Cov[, "YSO"]))] <- # Insert means for imputing PctDead inside ADS polygons
  mean(PctDead.b[which(!is.na(PctDead.b) & !is.na(Cov[, "YSO"]))])
PctDead.sd <- sd(PctDead.b[which(is.na(Cov[, "YSO"]))]) %>% rep(length(PctDead.b)) # SDs for imputing missing values outside ADS polygons
PctDead.sd[which(!is.na(Cov[, "YSO"]))] <- sd(PctDead.b[which(!is.na(Cov[, "YSO"]))]) # SDs for imputing missing values inside ADS polygons
PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values
PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed

YSO.b <- Cov[, "YSO"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
YSO.d <- tapply(YSO.b, gridID, mean, na.rm = T) # Grid-level values
YSO.b[is.na(YSO.b)] <- 0
YSO.d[is.na(YSO.d)] <- 0

Outbrk.b <- Cov[, "YSO"] %>% (function(x) as.integer(!is.na(x)))
Outbrk.d <- tapply(Outbrk.b, gridID, max)

TWIP.b <- Cov[, "TWIP"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
TWIP.d <- tapply(TWIP.b, gridID, mean, na.rm = T) # Grid-level values
TWIP.sd <- tapply(TWIP.b, gridID, sd, na.rm = T) # Grid-specific SDs for imputing missing values
TWIP.b.missing <- is.na(TWIP.b) %>% as.integer # Index missing values to be imputed
TWIP.b[is.na(TWIP.b)] <- 0

DOY.b <- Cov[, "DayOfYear"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values

Time.b <- Cov[, "Time"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
if(stratum == "LP")
  Time.b[point.list.LP == "LP-057-06"] <-
    mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
      # Impute missing time value in LP data based on survey times for point likely completed just before and just after.
        # Inferred which points after viewing the other points and their associated times in ArcMap.

ccov.b <- Cov[, "CanCov"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
ccov.means <- tapply(ccov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
ccov.sd <- tapply(ccov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
ccov.sd[which(ccov.sd == 0)] <- min(ccov.sd[which(ccov.sd > 0)]) # Zeros won't work for SD!
ccov.b.missing <- is.na(ccov.b) %>% as.integer # Index missing values to be imputed
ccov.b[is.na(ccov.b)] <- 0

shcov.b <- Cov[, "shrub_cover"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
shcov.means <- tapply(shcov.b, gridID, mean, na.rm = T) # Grid-level means for imputing missing values
shcov.sd <- tapply(shcov.b, gridID, sd, na.rm = T) # Grid-level SDs for imputing missing values
shcov.b.missing <- is.na(shcov.b) %>% as.integer # Index missing values to be imputed
shcov.b[is.na(shcov.b)] <- 0

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
