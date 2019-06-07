library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

#setwd("/home/RMBO.LOCAL/quresh.latif/CFLRP")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled_abundance.RData")

#_______ Functions _________#
dat.fn <- function(yso, pdead, mod) {
  x.yso <- rep(-1:18, 2)
  z.yso <- (x.yso - mean(yso, na.rm = T)) / sd(yso, na.rm = T)
  x.pd <- c(rep(0, 20), rep(100, 20))
  z.pd <- ((x.pd / 100) - mean(pdead, na.rm = T)) / sd(pdead, na.rm = T)
  dat.pred <- data.frame(x.yso = x.yso, z.yso = z.yso,
                         x.pd = x.pd, z.pd = z.pd,
                         pred = NA, pred.lo = NA, pred.hi = NA)
  
  B0 <- mod$sims.list$beta0.mean
  B.pdead <- mod$sims.list$bl.pdead
  B.yso <- mod$sims.list$bl.YSO
  B.yso2 <- mod$sims.list$bl.YSO2
  B.pdXYSO <- mod$sims.list$bl.pdXYSO
  for(i in 1:nrow(dat.pred)) {
    lambda <- exp(B0 + B.pdead*dat.pred$z.pd[i] + B.yso*dat.pred$z.yso[i] +
                     B.yso2*dat.pred$z.yso[i]^2 + B.pdXYSO*dat.pred$z.pd[i]*dat.pred$z.yso[i])
    dat.pred$pred[i] <- median(lambda)
    dat.pred$pred.lo[i] <- quantile(lambda, prob = 0.05, type = 8)
    dat.pred$pred.hi[i] <- quantile(lambda, prob = 0.95, type = 8)
  }
  return(dat.pred)
}
#___________________________#

#________Lodgepole pine plots ___________#
Spp <- c("AMRO", "DEJU", "HETH",
            "PISI", "RCKI", "YRWA")
strt <- "LP"
maxyso <- 12 # Set to 12 for LP and 9 for SF
Cov <- Cov.LP

pdd <- Cov[, "DeadConif"]
pdd[which(Cov[, "YSO"] > maxyso)] <- NA
yso <- Cov[, "YSO"]
yso <- yso[-which(is.na(yso))]

for(i in 1:length(Spp)) {
  spp <- Spp[i]
  mod <- loadObject(str_c("abund_models/mod_", spp, "_abundance_outbreak_", strt))
  dat <- dat.fn(yso, pdd, mod) %>%
    mutate(x.pd = as.factor(x.pd))
  maxy <- max(dat$pred.hi)
  p <- ggplot(dat, aes(x = x.yso, y = pred)) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi, fill = x.pd), alpha = 0.3) +
    geom_line(aes(color = x.pd), size = 1) +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    scale_color_manual(values = c("#009E73", "#D55E00")) +
    scale_fill_manual(values = c("#009E73", "#D55E00")) +
    xlab(NULL) + ylab(NULL) +
    guides(fill = F, color = F) +
    annotate("text", x = 1, y = maxy, label = spp)
  assign(str_c("pp", spp), p)
}

p1 <- ggdraw() +
  draw_plot(ppAMRO, y = 0.667, x = 0., width = 1, height = 0.333) +
  draw_plot(ppDEJU, y = 0.333, x = 0, width = 1, height = 0.333) +
  draw_plot(ppHETH, y = 0, x = 0, width = 1, height = 0.333)

p2 <- ggdraw() +
  draw_plot(ppPISI, y = 0.667, x = 0., width = 1, height = 0.25) +
  draw_plot(ppRCKI, y = 0.333, x = 0, width = 1, height = 0.25) +
  draw_plot(ppYRWA, y = 0, x = 0, width = 1, height = 0.25)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p2, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Mean abundance", "Years since outbreak"),
                  x = c(0, 0.43), y = c(0.47, 0.05),
                  angle = c(90, 0), size = c(20, 15), hjust = c(0, 0))

save_plot(str_c("Plot_spp_abundance_outbreak_", strt, ".tiff"), p, ncol = 2, nrow = 3, dpi = 200)
#________________________________________#

#________Spruce-fir plots ___________#
Spp <- c("GRAJ", "HETH", "MOCH",
         "RBNU", "RCKI")
strt <- "SF"
maxyso <- 9 # Set to 12 for LP and 9 for SF
Cov <- Cov.SF

pdd <- Cov[, "DeadConif"]
pdd[which(Cov[, "YSO"] > maxyso)] <- NA
yso <- Cov[, "YSO"]
yso <- yso[-which(is.na(yso))]

for(i in 1:length(Spp)) {
  spp <- Spp[i]
  mod <- loadObject(str_c("abund_models/mod_", spp, "_abundance_outbreak_", strt))
  dat <- dat.fn(yso, pdd, mod) %>%
    mutate(x.pd = as.factor(x.pd))
  maxy <- max(dat$pred.hi)
  p <- ggplot(dat, aes(x = x.yso, y = pred)) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi, fill = x.pd), alpha = 0.3) +
    geom_line(aes(color = x.pd), size = 1) +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    scale_color_manual(values = c("#009E73", "#D55E00")) +
    scale_fill_manual(values = c("#009E73", "#D55E00")) +
    xlab(NULL) + ylab(NULL) +
    guides(fill = F, color = F) +
    annotate("text", x = 1, y = maxy, label = spp)
  assign(str_c("pp", spp), p)
}

p1 <- ggdraw() +
  draw_plot(ppGRAJ, y = 0.667, x = 0., width = 1, height = 0.333) +
  draw_plot(ppHETH, y = 0.333, x = 0, width = 1, height = 0.333) +
  draw_plot(ppMOCH, y = 0, x = 0, width = 1, height = 0.333)

p2 <- ggdraw() +
  draw_plot(ppRBNU, y = 0.5, x = 0, width = 1, height = 0.5) +
  draw_plot(ppRCKI, y = 0, x = 0, width = 1, height = 0.5)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p2, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Mean abundance", "Years since outbreak"),
                  x = c(0, 0.43), y = c(0.47, 0.05),
                  angle = c(90, 0), size = c(20, 15), hjust = c(0, 0))

save_plot(str_c("Plot_spp_abundance_outbreak_", strt, ".tiff"), p, ncol = 2, nrow = 2, dpi = 200)
#____________________________________#