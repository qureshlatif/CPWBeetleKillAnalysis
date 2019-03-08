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
mod <- loadObject("mod_SFcommunity_outbreak_reduced3")
maxyso <- 9 # Set to 12 for LP and 9 for SF
#___________________________#
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval

#_______ Functions _________#
dat.pdead.fn <- function(ind.spp, pdead, mod) {
  x.pdead <- seq(min(pdead), max(pdead), length.out = 10)
  z.pdead <- (x.pdead - mean(pdead, na.rm = T)) / sd(pdead, na.rm = T)
  x.pdead <- x.pdead * 100
  dat.pred <- data.frame(x.pdead = x.pdead, z.pdead = z.pdead,
                         pred = NA, pred.lo = NA, pred.hi = NA)

  D0 <- mod$sims.list$d0[, ind.spp]
  B0 <- mod$sims.list$b0[, ind.spp]
  B.pdead <- mod$sims.list$bb.pdead[, ind.spp]
  for(i in 1:nrow(dat.pred)) {
    psi <- expit(D0)
    theta <- expit(B0 + B.pdead*dat.pred$z.pdead[i])
    dat.pred$pred[i] <- median(psi * theta)
    dat.pred$pred.lo[i] <- quantile(psi * theta, prob = 0.05, type = 8)
    dat.pred$pred.hi[i] <- quantile(psi * theta, prob = 0.95, type = 8)
  }
  return(dat.pred)
}

dat.yso.fn <- function(ind.spp, yso, mod) {
  yso <- yso[which(!is.na(yso))]
  x.yso <- seq(min(yso), max(yso))
  z.yso <- (x.yso - mean(yso, na.rm = T)) / sd(yso, na.rm = T)
  dat.pred <- data.frame(x.yso = x.yso, z.yso = z.yso,
                         pred = NA, pred.lo = NA, pred.hi = NA)
  
  D0 <- mod$sims.list$d0[, ind.spp]
  B0 <- mod$sims.list$b0[, ind.spp]
  B.yso <- mod$sims.list$bb.YSO[, ind.spp]
  B.yso2 <- mod$sims.list$bb.YSO2[, ind.spp]
  for(i in 1:nrow(dat.pred)) {
    psi <- expit(D0)
    theta <- expit(B0 + B.yso*dat.pred$z.yso[i] + B.yso2*dat.pred$z.yso[i]^2)
    dat.pred$pred[i] <- median(psi * theta)
    dat.pred$pred.lo[i] <- quantile(psi * theta, prob = 0.05, type = 8)
    dat.pred$pred.hi[i] <- quantile(psi * theta, prob = 0.95, type = 8)
  }
  return(dat.pred)
}
#___________________________#

##____ Group 1 - supported relations with pdead ____##
spp.plot <- c("ATTW", "COFL", "WAVI", "STJA", "RBNU", "PISI")

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  pdd <- Cov[, "DeadConif"]
  pdd[which(Cov[, "YSO"] > maxyso)] <- NA
  pdd <- pdd[-which(is.na(pdd))]
  dat.pdead <- dat.pdead.fn(ind.spp, pdd, mod)
  p1 <- ggplot(dat.pdead, aes(x = x.pdead, y = pred)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi), alpha = 0.4) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  yso <- Cov[, "YSO"]
  yso <- yso[-which(is.na(yso))]
  dat.yso <- dat.yso.fn(ind.spp, yso, mod)
  p2 <- ggplot(dat.yso, aes(x = x.yso, y = pred)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi), alpha = 0.4) +
    scale_x_continuous(breaks = seq(min(yso), max(yso), by = 2)) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p1, x = 0, y = 0, width = 0.5, height = 0.95) +
    draw_plot(p2, x = 0.5, y = 0, width = 0.5, height = 0.95) +
    draw_plot_label(c(spp, spp), x = c(0.2, 0.7), y = c(1, 1))
  assign(str_c("pp", i), p)
}

p <- ggdraw() +
  draw_plot(pp1, x = 0, y = 0.6666667, width = 0.5, height = 0.3333) +
  draw_plot(pp2, x = 0.5, y = 0.6666667, width = 0.5, height = 0.3333) +
  draw_plot(pp3, x = 0, y = 0.3333333, width = 0.5, height = 0.3333) +
  draw_plot(pp4, x = 0.5, y = 0.3333333, width = 0.5, height = 0.3333) +
  draw_plot(pp5, x = 0, y = 0, width = 0.5, height = 0.3333) +
  draw_plot(pp6, x = 0.5, y = 0, width = 0.5, height = 0.3333)

p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0.05, width = 0.95, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent dead conifer", "Years since outbreak",
                    "Percent dead conifer", "Years since outbreak"),
                  x = c(0, 0.13, 0.35, 0.6, 0.83), y = c(0.47, 0.05, 0.05, 0.05, 0.05),
                  angle = c(90, 0, 0, 0, 0), size = c(20, 15, 15, 15, 15), hjust = c(0, 0, 0, 0, 0))

save_plot(str_c("Plot_spp_outbreak_", stratum, "1.tiff"), p, ncol = 4, nrow = 3, dpi = 200)

##____ Group 2 - supported relations with YSO but not % dead ____##
spp.plot <- c("BTLH", "NOFL", "GRAJ", "GCKI",
              "RCKI", "HETH", "PIGR", "CAFI",
              "RECR", "CHSP", "LISP", "WCSP",
              "DEJU", "YRWA", "WETA")

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  ind.spp <- which(spp.list == spp)
  
  pdd <- Cov[, "DeadConif"]
  pdd[which(Cov[, "YSO"] > maxyso)] <- NA
  pdd <- pdd[-which(is.na(pdd))]
  dat.pdead <- dat.pdead.fn(ind.spp, pdd, mod)
  p1 <- ggplot(dat.pdead, aes(x = x.pdead, y = pred)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi), alpha = 0.4) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  yso <- Cov[, "YSO"]
  yso <- yso[-which(is.na(yso))]
  dat.yso <- dat.yso.fn(ind.spp, yso, mod)
  p2 <- ggplot(dat.yso, aes(x = x.yso, y = pred)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi), alpha = 0.4) +
    scale_x_continuous(breaks = seq(min(yso), max(yso), by = 2)) +
    ylim(0, 1) +
    xlab(NULL) + ylab(NULL)
  
  p <- ggdraw() +
    draw_plot(p1, x = 0, y = 0, width = 0.5, height = 0.95) +
    draw_plot(p2, x = 0.5, y = 0, width = 0.5, height = 0.95) +
    draw_plot_label(c(spp, spp), x = c(0.2, 0.7), y = c(1, 1))
  assign(str_c("pp", i), p)
}

p <- ggdraw() +
  draw_plot(pp1, x = 0, y = 0.875, width = 0.5, height = 0.125) +
  draw_plot(pp2, x = 0, y = 0.750, width = 0.5, height = 0.125) +
  draw_plot(pp3, x = 0.5, y = 0.750, width = 0.5, height = 0.125) +
  draw_plot(pp4, x = 0, y = 0.625, width = 0.5, height = 0.125) +
  draw_plot(pp5, x = 0.5, y = 0.625, width = 0.5, height = 0.125) +
  draw_plot(pp6, x = 0, y = 0.500, width = 0.5, height = 0.125) +
  draw_plot(pp7, x = 0.5, y = 0.500, width = 0.5, height = 0.125) +
  draw_plot(pp8, x = 0, y = 0.375, width = 0.5, height = 0.125) +
  draw_plot(pp9, x = 0.5, y = 0.375, width = 0.5, height = 0.125) +
  draw_plot(pp10, x = 0, y = 0.250, width = 0.5, height = 0.125) +
  draw_plot(pp11, x = 0.5, y = 0.250, width = 0.5, height = 0.125) +
  draw_plot(pp12, x = 0, y = 0.125, width = 0.5, height = 0.125) +
  draw_plot(pp13, x = 0.5, y = 0.125, width = 0.5, height = 0.125) +
  draw_plot(pp14, x = 0, y = 0, width = 0.5, height = 0.125) +
  draw_plot(pp15, x = 0.5, y = 0, width = 0.5, height = 0.125)
  
p <- ggdraw() +
  draw_plot(p, x = 0.05, y = 0.05, width = 0.95, height = 0.95) +
  draw_plot_label(c("Occupancy", "Percent dead conifer", "Years since outbreak",
                    "Percent dead conifer", "Years since outbreak"),
                  x = c(0, 0.13, 0.35, 0.6, 0.83), y = c(0.47, 0.05, 0.05, 0.05, 0.05),
                  angle = c(90, 0, 0, 0, 0), size = c(20, 15, 15, 15, 15), hjust = c(0, 0, 0, 0, 0))

save_plot(str_c("Plot_spp_outbreak_", stratum, "2.tiff"), p, ncol = 4, nrow = 6, dpi = 200)
