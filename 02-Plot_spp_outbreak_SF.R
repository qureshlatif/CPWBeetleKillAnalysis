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

spp.plot <- c("ATTW", "WAVI", "STJA", "CLNU", "CORA",
              "RBNU", "GCKI", "HETH", "PISI", "DEJU")
relat.labs <- c("DCon+", "DCon-", "DCon-", "DConXYSO-", "DConXYSO-",
                "DCon-", "DConXYSO-", "YSO-", "DCon-,YSOlag", "YSOlag")

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
  draw_plot(ppATTW, y = 0.8, x = 0., width = 1, height = 0.2) +
  draw_plot(ppSTJA, y = 0.6, x = 0, width = 1, height = 0.2) +
  draw_plot(ppCORA, y = 0.4, x = 0, width = 1, height = 0.2) +
  draw_plot(ppGCKI, y = 0.2, x = 0, width = 1, height = 0.2) +
  draw_plot(ppPISI, y = 0, x = 0, width = 1, height = 0.2)

p2 <- ggdraw() +
  draw_plot(ppWAVI, y = 0.8, x = 0., width = 1, height = 0.2) +
  draw_plot(ppCLNU, y = 0.6, x = 0, width = 1, height = 0.2) +
  draw_plot(ppRBNU, y = 0.4, x = 0, width = 1, height = 0.2) +
  draw_plot(ppHETH, y = 0.2, x = 0, width = 1, height = 0.2) +
  draw_plot(ppDEJU, y = 0, x = 0, width = 1, height = 0.2)

p <- ggdraw() +
  draw_plot(p1, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p2, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Occupancy", "Years since outbreak"),
                  x = c(0, 0.43), y = c(0.47, 0.05),
                  angle = c(90, 0), size = c(20, 15), hjust = c(0, 0))

save_plot(str_c("Plot_spp_outbreak_", stratum, ".tiff"), p, ncol = 2, nrow = 3, dpi = 200)
