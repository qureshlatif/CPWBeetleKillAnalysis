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
stratum <- "SF"
mod <- loadObject("mod_SFcommunity_outbreak_reduced2")
maxyso <- 9 # Set to 12 for LP and 9 for SF
theme_set(theme_cowplot())
#___________________________#
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval

#_______ Functions _________#
dat.fn <- function(ind.spp, pdead, mod) {
  x.pd <- seq(0, 100, length.out = 20)
  z.pd <- ((x.pd / 100) - mean(pdead, na.rm = T)) / sd(pdead, na.rm = T)
  dat.pred <- data.frame(x.pd = x.pd, z.pd = z.pd,
                         pred = NA, pred.lo = NA, pred.hi = NA)
  
  D0 <- mod$sims.list$d0[, ind.spp]
  B0 <- mod$sims.list$b0[, ind.spp]
  B.pdead <- mod$sims.list$bb.pdead[, ind.spp]
  for(i in 1:nrow(dat.pred)) {
    psi <- expit(D0)
    theta <- expit(B0 + B.pdead*dat.pred$z.pd[i])
    dat.pred$pred[i] <- median(psi * theta)
    dat.pred$pred.lo[i] <- quantile(psi * theta, prob = 0.05, type = 8)
    dat.pred$pred.hi[i] <- quantile(psi * theta, prob = 0.95, type = 8)
  }
  return(dat.pred)
}
#___________________________#

spp.plot <- c("STJA", "WAVI", "RBNU")

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  pdd <- Cov[, "DeadConif"]
  dat <- dat.fn(ind.spp, pdd, mod)
  p <- ggplot(dat, aes(x = x.pd, y = pred)) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi), alpha = 0.3) +
    geom_line(size = 1) +
    ylim(0, 0.31) +
    xlab(NULL) + ylab(NULL)
  assign(str_c("pp", spp), p)
}

