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
stratum <- "LP"
mod <- loadObject("mod_LPcommunity_outbreak_reduced2")
maxyso <- 12 # Set to 12 for LP and 9 for SF
theme_set(theme_cowplot())
#___________________________#
Cov <- str_c("Cov.", stratum) %>% as.name %>% eval

#_______ Functions _________#
dat.fn <- function(ind.spp, yso, pdead, mod) {
  x.yso <- rep(-1:18, 2)
  z.yso <- (x.yso - mean(yso, na.rm = T)) / sd(yso, na.rm = T)
  x.pd <- c(rep(0, 20), rep(100, 20))
  z.pd <- ((x.pd / 100) - mean(pdead, na.rm = T)) / sd(pdead, na.rm = T)
  dat.pred <- data.frame(x.yso = x.yso, z.yso = z.yso,
                         x.pd = x.pd, z.pd = z.pd,
                         pred = NA, pred.lo = NA, pred.hi = NA)

  D0 <- mod$sims.list$d0[, ind.spp]
  B0 <- mod$sims.list$b0[, ind.spp]
  B.pdead <- mod$sims.list$bb.pdead[, ind.spp]
  B.yso <- mod$sims.list$bb.YSO[, ind.spp]
  B.yso2 <- mod$sims.list$bb.YSO2[, ind.spp]
  B.pdXYSO <- mod$sims.list$bb.pdXYSO[, ind.spp]
  for(i in 1:nrow(dat.pred)) {
    psi <- expit(D0)
    theta <- expit(B0 + B.pdead*dat.pred$z.pd[i] + B.yso*dat.pred$z.yso[i] +
                     B.yso2*dat.pred$z.yso[i]^2 + B.pdXYSO*dat.pred$z.pd[i]*dat.pred$z.yso[i])
    dat.pred$pred[i] <- median(psi * theta)
    dat.pred$pred.lo[i] <- quantile(psi * theta, prob = 0.05, type = 8)
    dat.pred$pred.hi[i] <- quantile(psi * theta, prob = 0.95, type = 8)
  }
  return(dat.pred)
}
#___________________________#

spp.plot <- c("MODO", "ATTW", "OSFL", "WEWP", "DUFL", "COFL", "VGSW",
              "HOWR", "RCKI", "TOSO", "HETH", "AMRO", "GTTO", "LISP",
              "WCSP", "DEJU", "BHCO", "MGWA", "YRWA", "WIWA", "WETA")
relat.labs <- c("YSO+", "YSO+", "YSO+", "YSO+", "YSO+", "YSO+", "YSO-",
              "YSO+", "DCon-", "YSOlag", "YSOpeak", "YSO+", "YSOlag", "YSO+",
              "YSO+", "YSO+", "YSO+", "YSO+", "DConXYSO", "YSO+", "YSOlag")

for(i in 1:length(spp.plot)) {
  spp <- spp.plot[i]
  rlab <- relat.labs[i]
  ind.spp <- which(spp.list == spp)
  
  pdd <- Cov[, "DeadConif"]
  pdd[which(Cov[, "YSO"] > maxyso)] <- NA
  yso <- Cov[, "YSO"]
  yso <- yso[-which(is.na(yso))]
  dat <- dat.fn(ind.spp, yso, pdd, mod) %>%
    mutate(x.pd = as.factor(x.pd))
  p <- ggplot(dat, aes(x = x.yso, y = pred)) +
    geom_ribbon(aes(ymin = pred.lo, ymax = pred.hi, fill = x.pd), alpha = 0.3) +
    geom_line(aes(color = x.pd), size = 1) +
    ylim(0, 1) +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    scale_color_manual(values = c("#009E73", "#D55E00")) +
    scale_fill_manual(values = c("#009E73", "#D55E00")) +
    xlab(NULL) + ylab(NULL) +
    guides(fill = F, color = F) +
    annotate("text", x = 7.5, y = 1, label = spp) +
    annotate("text", x = 7.5, y = 0.8, label = rlab)
  assign(str_c("pp", spp), p)
}

p1 <- ggdraw() +
  draw_plot(ppMODO, y = 0.8571429, x = 0., width = 1, height = 0.1428571) +
  draw_plot(ppWEWP, y = 0.7142857, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppVGSW, y = 0.5714286, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppTOSO, y = 0.4285714, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppGTTO, y = 0.2857143, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppDEJU, y = 0.1428571, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppYRWA, y = 0, x = 0, width = 1, height = 0.1428571)

p2 <- ggdraw() +
  draw_plot(ppATTW, y = 0.8571429, x = 0., width = 1, height = 0.1428571) +
  draw_plot(ppDUFL, y = 0.7142857, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppHOWR, y = 0.5714286, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppHETH, y = 0.4285714, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppLISP, y = 0.2857143, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppBHCO, y = 0.1428571, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppWIWA, y = 0, x = 0, width = 1, height = 0.1428571)

p3 <- ggdraw() +
  draw_plot(ppOSFL, y = 0.8571429, x = 0., width = 1, height = 0.1428571) +
  draw_plot(ppCOFL, y = 0.7142857, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppRCKI, y = 0.5714286, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppAMRO, y = 0.4285714, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppWCSP, y = 0.2857143, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppMGWA, y = 0.1428571, x = 0, width = 1, height = 0.1428571) +
  draw_plot(ppWETA, y = 0, x = 0, width = 1, height = 0.1428571)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.3166667, height = 0.95) +
  draw_plot(p2, x = 0.3666667, y = 0.05, width = 0.3166667, height = 0.95) +
  draw_plot(p3, x = 0.6833333, y = 0.05, width = 0.3166667, height = 0.95) +
  draw_plot_label(c("Occupancy", "Years since outbreak"),
                  x = c(0, 0.43), y = c(0.47, 0.05),
                  angle = c(90, 0), size = c(20, 15), hjust = c(0, 0))

save_plot(str_c("Plot_spp_outbreak_", stratum, ".tiff"), p, ncol = 2, nrow = 3, dpi = 200)
