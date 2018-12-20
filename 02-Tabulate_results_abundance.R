library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled_abundance.RData")

#_____ Inputs ______#
SPP <- c("ATTW", "MOCH", "RBNU")
#___________________#

#____ Convergence ____#
cols <- c("strata", "min_n.eff", "max_Rhat")
out <- matrix(NA, nrow = length(SPP), ncol = length(cols))
dimnames(out) <- list(SPP, cols)

for(spp in SPP) {
  mod <- loadObject(str_c("abund_models/mod_", spp,"_abundance_outbreak"))
  out[spp, "strata"] <- "SF"
  out[spp, "min_n.eff"] <- min(mod$summary[, "n.eff"])
  out[spp, "max_Rhat"] <- max(mod$summary[, "Rhat"])
}

write.csv(out, "abund_models/Convergence_summary.csv", row.names = T)
#____________________#

parameters <- c("beta0.mean", "beta0.sd", "N.mean", "p.mean", # Assemble the parameters vector for JAGS (What we want to track).
                "beta0", "N", "cl.size",
                "bd.heatload", "bd.TWI",
                "bl.pdead", "bl.YSO", "bl.YSO2", "bl.pdXYSO",
                "bl.RCovAS", "bl.RCovES", "bl.RCovPine",
                ##___ Hazard rate parameters ___##
                "a0", "a.Time", "a.Time2",
                "a.DOY", "a.DOY2",
                "a.pdead", "a.YSO",
                "b",
                ##______________________________##
                "lambda", "a", # Needed for WAIC
                "test") # GOF

inits <- function() # Setting these based on posterior distribution from an initial successful run.
  list(N = Y, beta0.mean = rnorm(1, 1, 0.1), beta0.sd = rnorm(1, 0.66, 0.07),
       a0 = rnorm(1, 3.5, 0.5), b = rnorm(1, 3.16, 0.5)) # for hazard rate model

for(strt in str_sub(strata, -2, -1)) {
  SPP <- eval(as.name(str_c("Spp_", strt)))
  ifelse(strt == "LP", maxYSOForPD <- 12, maxYSOForPD <- 9) # Set to 12 for LP and 9 for SF
  model.file <- str_c(modfile.loc, "model_abundance_outbreak_", strt, ".jags")
  
  # Covariates #
  Cov <- eval(as.name(str_c("Cov.", strt)))
  
  gridID <- Cov[, "gridIndex"]
  nGrid <- max(gridID)

  PctDead.b <- Cov[, "DeadConif"] # Point-level values
  PctDead.b[which(Cov[, "YSO"] > maxYSOForPD)] <- NA # Drop values from later years when (presumably) %Dead starts reflecting snag fall 
  PctDead.b <- PctDead.b %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Re-scale point values
  PctDead.b[which(is.na(PctDead.b) & is.na(Cov[, "YSO"]))] <- # Insert means for imputing PctDead outside ADS polygons
    mean(PctDead.b[which(!is.na(PctDead.b) & is.na(Cov[, "YSO"]))])
  PctDead.b[which(is.na(PctDead.b) & !is.na(Cov[, "YSO"]))] <- # Insert means for imputing PctDead inside ADS polygons
    mean(PctDead.b[which(!is.na(PctDead.b) & !is.na(Cov[, "YSO"]))])
  PctDead.sd <- sd(PctDead.b[which(is.na(Cov[, "YSO"]))]) %>% rep(length(PctDead.b)) # SDs for imputing missing values outside ADS polygons
  PctDead.sd[which(!is.na(Cov[, "YSO"]))] <- sd(PctDead.b[which(!is.na(Cov[, "YSO"]))]) # SDs for imputing missing values inside ADS polygons
  PctDead.lower <- min(PctDead.b, na.rm = T) # Lower bound for imputing missing values
  PctDead.b.missing <- is.na(PctDead.b) %>% as.integer # Index missing values to be imputed
  
  outbreak_grids <- tapply(Cov[, "YSO"], Cov[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
  YSO.b <- Cov[, "YSO"] %>% (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)) # Point-level values
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
  
  heatload.d <- Cov[, "heatload"]  %>% # Grid-level values only
    tapply(gridID, mean, na.rm = T) %>%
    (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
  
  TWI.d <- Cov[, "TWI"]  %>% # Grid-level values only
    tapply(gridID, mean, na.rm = T) %>%
    (function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
  
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
  if(strt == "LP") Time.b[point.list.LP == "LP-057-06"] <- 
    mean(c(Time.b[point.list.LP == "LP-057-07"], Time.b[point.list.LP == "LP-057-03"]))
  # Imputed missing time value in LP data based on survey times for point likely completed just before and just after.
  # Inferred which points after viewing the other points and their associated times in ArcMap.
  
  for(spp in SPP) {
    spp.ind <- which(SPP == spp)
    
    save.out <- str_c(save.out.loc, "mod_", spp, "_abundance_outbreak")
    
    area.band <- eval(as.name(str_c("area.band.", strt)))[[spp.ind]]
    area.prop <- eval(as.name(str_c("area.prop.", strt)))[[spp.ind]]
    breaks <- eval(as.name(str_c("breaks.", strt)))[[spp.ind]]
    area.circle <- eval(as.name(str_c("area.circle.", strt)))[spp.ind]
    cutoff <- eval(as.name(str_c("cutoff.", strt)))[spp.ind]
    
    Y <- eval(as.name(str_c("Y.", spp, ".", strt, ".dist")))
    dclass <- eval(as.name(str_c("dclass.", spp, ".", strt)))
    mean.cl <- max(1.001, mean(dclass[, "CL_Count"]))
    sd.cl <- max(0.001, sd(dclass[, "CL_Count"]))
    dimnames(dclass) <- NULL
    nInd <- nrow(dclass)
    nPoint <- length(Y)
    
    st.time <- Sys.time()
    out <- jagsUI(data, inits, parameters.to.save = parameters, model.file, n.thin=nt, n.chains=nc,
                  n.burnin=nb, n.iter=ni, parallel=TRUE)
    end.time <- Sys.time()
    run.time <- end.time - st.time
    run.time
    rm(st.time,end.time)
    
    library(R.utils)
    saveObject(out, save.out)
  }
}


# max(out$summary[which(!is.na(out$summary[ ,"Rhat"])) ,"Rhat"])
# sort(out$summary[,"Rhat"], decreasing = T)[1:100]
# 
# min(out$summary[,"n.eff"])
# sort(out$summary[,"n.eff"])[1:100]
# 
# sum(out$sims.list$test) / out$mcmc.info$n.samples

