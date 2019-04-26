library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

mod_LP <- loadObject("mod_LPcommunity_habitat_reduced")
mod_SF <- loadObject("mod_SFcommunity_habitat_reduced")
spp.outbreak <- c("MODO", "BTLH", "ATTW", "NOFL",
                  "OSFL", "WEWP", "COFL", "WAVI",
                  "GRAJ", "STJA", "RBNU", "GCKI",
                  "RCKI", "TOSO", "HETH", "AMRO",
                  "PIGR", "CAFI", "RECR", "PISI",
                  "GTTO", "CHSP", "LISP", "WCSP",
                  "DEJU", "YRWA", "WETA")

# Tabulate parameter estimates
pars <- c("theta", "bb.RCovAS", "bb.RCovPine", "bb.RCovES", "bb.CanCov", "bb.ShCov", "bb.RSC_Con", "bb.GHerb", "bb.Gwoody", "bb.GDD")
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           select(Var2, Var1) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
parest_LP <- parest_SF <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

for(par in pars[-which(pars %in% c("theta", "bb.RCovES"))]) {
  parm <- mod_LP$sims.list[[par]]
  parest_LP[, par] <- apply(parm, 2, median)
  parest_LP[, str_c(par, ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.05, type = 8))
  parest_LP[, str_c(par, ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.95, type = 8))
}
parm <- expit(mod_LP$sims.list[["b0"]])
parest_LP[, "theta"] <- apply(parm, 2, median)
parest_LP[, "theta.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.05, type = 8))
parest_LP[, "theta.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.95, type = 8))

for(par in pars[-which(pars %in% c("theta"))]) {
  parm <- mod_SF$sims.list[[par]]
  parest_SF[, par] <- apply(parm, 2, median)
  parest_SF[, str_c(par, ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.05, type = 8))
  parest_SF[, str_c(par, ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.95, type = 8))
}
parm <- expit(mod_SF$sims.list[["b0"]])
parest_SF[, "theta"] <- apply(parm, 2, median)
parest_SF[, "theta.lo"] <- apply(parm, 2, function(x) quantile(x, prob = 0.05, type = 8))
parest_SF[, "theta.hi"] <- apply(parm, 2, function(x) quantile(x, prob = 0.95, type = 8))

rm(parm)

beta.cols <- dimnames(parest_LP)[[2]][-which(dimnames(parest_LP)[[2]] %in% c("theta", "theta.lo", "theta.hi"))]
ind.spp <- c(which(spp.list %in% spp.outbreak),
             parest_LP[, beta.cols[(beta.cols %>% str_detect(".lo") & !beta.cols %>% str_detect("RCovES")) %>% which]] %>%
               apply(1, function(x) any(x > 0)) %>%
               which,
             parest_LP[, beta.cols[(beta.cols %>% str_detect(".hi") & !beta.cols %>% str_detect("RCovES")) %>% which]] %>%
               apply(1, function(x) any(x < 0)) %>%
               which,
             parest_SF[, beta.cols[(beta.cols %>% str_detect(".lo") & !beta.cols %>% str_detect("RCovPine")) %>% which]] %>%
               apply(1, function(x) any(x > 0)) %>%
               which,
             parest_SF[, beta.cols[(beta.cols %>% str_detect(".hi") & !beta.cols %>% str_detect("RCovPine")) %>% which]] %>%
               apply(1, function(x) any(x < 0)) %>%
               which) %>% unique %>% sort %>% length
spp.plt <- spp.list[ind.spp]

dat.plt.LP <- parest_LP %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  filter(Spp %in% spp.plt) %>%
  mutate(index = row_number()) %>%
  mutate(index = (max(index) - index) + 1) %>%
  select(theta:bb.RCovPine.hi, bb.CanCov:index)
dat.plt.LP$Spp[which(dat.plt.LP$Spp %in% spp.outbreak)] <-
  dat.plt.LP$Spp[which(dat.plt.LP$Spp %in% spp.outbreak)] %>% str_c("*")

cols <- pars[-c(1, 4)] %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.LP), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.LP[, which(str_detect(names(dat.plt.LP), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}

dat.plt.LP <- dat.plt.LP %>%
  bind_cols(dat.supp %>% as.data.frame)

dat.plt.SF <- parest_SF %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  filter(Spp %in% spp.plt) %>%
  mutate(index = row_number()) %>%
  mutate(index = (max(index) - index) + 1)
dat.plt.SF$Spp[which(dat.plt.SF$Spp %in% spp.outbreak)] <-
  dat.plt.SF$Spp[which(dat.plt.SF$Spp %in% spp.outbreak)] %>% str_c("*")

cols <- pars[-1] %>% str_c(".supp")
dat.supp <- matrix("none", nrow = nrow(dat.plt.SF), ncol = length(cols),
                   dimnames = list(NULL, cols))
for(i in 1:length(cols)) {
  col.chck <- str_sub(cols[i], 1, -6)
  chck <- dat.plt.SF[, which(str_detect(names(dat.plt.SF), col.chck))]
  dat.supp[which(chck[, 2] > 0), cols[i]] <- "pos"
  dat.supp[which(chck[, 3] < 0), cols[i]] <- "neg"
}

dat.plt.SF <- dat.plt.SF %>%
  bind_cols(dat.supp %>% as.data.frame)

rm(dat.supp, col.chck, chck)

#### Summarize relationships by covariate ####
cols <- c("LP_pos", "LP_neg", "LP_tot", "SF_pos", "SF_neg", "SF_tot")
sum.tab <- matrix(0, nrow = length(pars[-1]), ncol = length(cols),
                  dimnames = list(pars[-1], cols))
for(i in pars[-1]) {
  if(any(str_detect(dimnames(parest_LP)[[2]], i))) {
    lo <- parest_LP[, which(str_detect(dimnames(parest_LP)[[2]], i) & str_detect(dimnames(parest_LP)[[2]], ".lo"))]
    hi <- parest_LP[, which(str_detect(dimnames(parest_LP)[[2]], i) & str_detect(dimnames(parest_LP)[[2]], ".hi"))]
    sum.tab[i, "LP_pos"] <- sum(lo > 0)
    sum.tab[i, "LP_neg"] <- sum(hi < 0)
  }
  if(any(str_detect(dimnames(parest_SF)[[2]], i))) {
    lo <- parest_SF[, which(str_detect(dimnames(parest_SF)[[2]], i) & str_detect(dimnames(parest_SF)[[2]], ".lo"))]
    hi <- parest_SF[, which(str_detect(dimnames(parest_SF)[[2]], i) & str_detect(dimnames(parest_SF)[[2]], ".hi"))]
    sum.tab[i, "SF_pos"] <- sum(lo > 0)
    sum.tab[i, "SF_neg"] <- sum(hi < 0)
  }
}
sum.tab[, "LP_tot"] <- sum.tab[, c("LP_pos", "LP_neg")] %>% apply(1, sum)
sum.tab[, "SF_tot"] <- sum.tab[, c("SF_pos", "SF_neg")] %>% apply(1, sum)
rm(i, lo, hi)

#### Plots ####
## Lodgepole ##
p.theta <- ggplot(dat = dat.plt.LP, aes(x = index, y = theta)) +
  geom_errorbar(aes(ymin = theta.lo, ymax = theta.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(0, 1)) +
  ylab(expression(hat(theta)["mean"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25))

min.y <- min(c(dat.plt.LP$bb.RCovAS.lo, dat.plt.LP$bb.RCovPine.lo, dat.plt.LP$bb.CanCov.lo, dat.plt.LP$bb.ShCov.lo,
               dat.plt.LP$bb.RSC_Con.lo, dat.plt.LP$bb.GHerb.lo, dat.plt.LP$bb.Gwoody.lo, dat.plt.LP$bb.GDD.lo))
max.y <- max(c(dat.plt.LP$bb.RCovAS.hi, dat.plt.LP$bb.RCovPine.hi, dat.plt.LP$bb.CanCov.hi, dat.plt.LP$bb.ShCov.hi,
               dat.plt.LP$bb.RSC_Con.hi, dat.plt.LP$bb.GHerb.hi, dat.plt.LP$bb.Gwoody.hi, dat.plt.LP$bb.GDD.hi))

p.CanCov <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.CanCov, color = bb.CanCov.supp)) +
  geom_errorbar(aes(ymin = bb.CanCov.lo, ymax = bb.CanCov.hi, color = bb.CanCov.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["CanCov"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RCovAS <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.RCovAS, color = bb.RCovAS.supp)) +
  geom_errorbar(aes(ymin = bb.RCovAS.lo, ymax = bb.RCovAS.hi, color = bb.RCovAS.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Aspen"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RCovPine <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.RCovPine, color = bb.RCovPine.supp)) +
  geom_errorbar(aes(ymin = bb.RCovPine.lo, ymax = bb.RCovPine.hi, color = bb.RCovPine.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Pine"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.ShCov <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.ShCov, color = bb.ShCov.supp)) +
  geom_errorbar(aes(ymin = bb.ShCov.lo, ymax = bb.ShCov.hi, color = bb.ShCov.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["ShrubCov"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RSC_Con <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.RSC_Con, color = bb.RSC_Con.supp)) +
  geom_errorbar(aes(ymin = bb.RSC_Con.lo, ymax = bb.RSC_Con.hi, color = bb.RSC_Con.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["ConShrb"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.GHerb <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.GHerb, color = bb.GHerb.supp)) +
  geom_errorbar(aes(ymin = bb.GHerb.lo, ymax = bb.GHerb.hi, color = bb.GHerb.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Herb"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.Gwoody <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.Gwoody, color = bb.Gwoody.supp)) +
  geom_errorbar(aes(ymin = bb.Gwoody.lo, ymax = bb.Gwoody.hi, color = bb.Gwoody.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Woody"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.GDD <- ggplot(dat = dat.plt.LP, aes(x = index, y = bb.GDD, color = bb.GDD.supp)) +
  geom_errorbar(aes(ymin = bb.GDD.lo, ymax = bb.GDD.hi, color = bb.GDD.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.LP), labels = dat.plt.LP$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["DeadDown"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.theta, x = 0.0000000, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.CanCov, x = 0.1111111, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.RCovAS, x = 0.2222222, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.RCovPine, x = 0.3333333, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.ShCov, x = 0.4444444, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.RSC_Con, x = 0.5555556, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.GHerb, x = 0.6666667, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.Gwoody, x = 0.7777778, y = 0, width = 0.1111111, height = 1) +
  draw_plot(p.GDD, x = 0.8888889, y = 0, width = 0.1111111, height = 1)
save_plot("Plot_habitat_effects_LP.tiff", p, ncol = 6, nrow = 3, dpi = 200)


## Spruce-fir ##
p.theta <- ggplot(dat = dat.plt.SF, aes(x = index, y = theta)) +
  geom_errorbar(aes(ymin = theta.lo, ymax = theta.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(0, 1)) +
  ylab(expression(hat(theta)["mean"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25))

min.y <- min(c(dat.plt.SF$bb.RCovAS.lo, dat.plt.SF$bb.RCovPine.lo, dat.plt.SF$bb.RCovES.lo, dat.plt.SF$bb.CanCov.lo, dat.plt.SF$bb.ShCov.lo,
               dat.plt.SF$bb.RSC_Con.lo, dat.plt.SF$bb.GHerb.lo, dat.plt.SF$bb.Gwoody.lo, dat.plt.SF$bb.GDD.lo))
max.y <- max(c(dat.plt.SF$bb.RCovAS.hi, dat.plt.SF$bb.RCovPine.hi, dat.plt.SF$bb.RCovES.hi, dat.plt.SF$bb.CanCov.hi, dat.plt.SF$bb.ShCov.hi,
               dat.plt.SF$bb.RSC_Con.hi, dat.plt.SF$bb.GHerb.hi, dat.plt.SF$bb.Gwoody.hi, dat.plt.SF$bb.GDD.hi))

p.CanCov <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.CanCov, color = bb.CanCov.supp)) +
  geom_errorbar(aes(ymin = bb.CanCov.lo, ymax = bb.CanCov.hi, color = bb.CanCov.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["CanCov"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RCovAS <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.RCovAS, color = bb.RCovAS.supp)) +
  geom_errorbar(aes(ymin = bb.RCovAS.lo, ymax = bb.RCovAS.hi, color = bb.RCovAS.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Aspen"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RCovPine <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.RCovPine, color = bb.RCovPine.supp)) +
  geom_errorbar(aes(ymin = bb.RCovPine.lo, ymax = bb.RCovPine.hi, color = bb.RCovPine.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Pine"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RCovES <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.RCovES, color = bb.RCovES.supp)) +
  geom_errorbar(aes(ymin = bb.RCovES.lo, ymax = bb.RCovES.hi, color = bb.RCovES.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Spruce"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.ShCov <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.ShCov, color = bb.ShCov.supp)) +
  geom_errorbar(aes(ymin = bb.ShCov.lo, ymax = bb.ShCov.hi, color = bb.ShCov.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["ShrubCov"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.RSC_Con <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.RSC_Con, color = bb.RSC_Con.supp)) +
  geom_errorbar(aes(ymin = bb.RSC_Con.lo, ymax = bb.RSC_Con.hi, color = bb.RSC_Con.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["ConShrb"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.GHerb <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.GHerb, color = bb.GHerb.supp)) +
  geom_errorbar(aes(ymin = bb.GHerb.lo, ymax = bb.GHerb.hi, color = bb.GHerb.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Herb"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.Gwoody <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.Gwoody, color = bb.Gwoody.supp)) +
  geom_errorbar(aes(ymin = bb.Gwoody.lo, ymax = bb.Gwoody.hi, color = bb.Gwoody.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["Woody"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.GDD <- ggplot(dat = dat.plt.SF, aes(x = index, y = bb.GDD, color = bb.GDD.supp)) +
  geom_errorbar(aes(ymin = bb.GDD.lo, ymax = bb.GDD.hi, color = bb.GDD.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt.SF), labels = dat.plt.SF$Spp %>% rev, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("dark gray", "#D55E00")) +
  ylab(expression(hat(beta)["DeadDown"])) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.theta, x = 0, y = 0, width = 0.1, height = 1) +
  draw_plot(p.CanCov, x = 0.1, y = 0, width = 0.1, height = 1) +
  draw_plot(p.RCovAS, x = 0.2, y = 0, width = 0.1, height = 1) +
  draw_plot(p.RCovES, x = 0.3, y = 0, width = 0.1, height = 1) +
  draw_plot(p.RCovPine, x = 0.4, y = 0, width = 0.1, height = 1) +
  draw_plot(p.ShCov, x = 0.5, y = 0, width = 0.1, height = 1) +
  draw_plot(p.RSC_Con, x = 0.6, y = 0, width = 0.1, height = 1) +
  draw_plot(p.GHerb, x = 0.7, y = 0, width = 0.1, height = 1) +
  draw_plot(p.Gwoody, x = 0.8, y = 0, width = 0.1, height = 1) +
  draw_plot(p.GDD, x = 0.9, y = 0, width = 0.1, height = 1)
save_plot("Plot_habitat_effects_SF.tiff", p, ncol = 6, nrow = 3, dpi = 200)
