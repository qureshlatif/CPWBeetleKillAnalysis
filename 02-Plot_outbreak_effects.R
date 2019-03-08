library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

mod_LP <- loadObject("mod_LPcommunity_outbreak_reduced3")
mod_SF <- loadObject("mod_SFcommunity_outbreak_reduced3")

# Tabulate parameter estimates
pars <- c("bb.pdead", "bb.YSO", "bb.YSO2")
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           select(Var2, Var1) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
parest_LP <- parest_SF <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

for(i in 1:length(pars)) {
  parm <- mod_LP$sims.list[[pars[i]]]
  parest_LP[, pars[i]] <- apply(parm, 2, median)
  parest_LP[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.05, type = 8))
  parest_LP[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.95, type = 8))
}
for(i in 1:length(pars)) {
  parm <- mod_SF$sims.list[[pars[i]]]
  parest_SF[, pars[i]] <- apply(parm, 2, median)
  parest_SF[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.05, type = 8))
  parest_SF[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.95, type = 8))
}
rm(parm)

#### Plot all species outbreak effects ####
## Percent dead, lodgepole ##
dat <- parest_LP %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.pdead) %>%
  mutate(index = row_number())

p.pdead.LP <- ggplot(dat = dat, aes(x = index, y = bb.pdead)) +
  geom_errorbar(aes(ymin = bb.pdead.lo, ymax = bb.pdead.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.pdead.lo), max(dat$bb.pdead.hi))) +
  ylab(expression(hat(beta)["pdead"])) + xlab(NULL)

## Years since outbreak, lodgepole ##
dat <- parest_LP %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.YSO) %>%
  mutate(index = row_number())

p.YSO.LP <- ggplot(dat = dat, aes(x = index, y = bb.YSO)) +
  geom_errorbar(aes(ymin = bb.YSO.lo, ymax = bb.YSO.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.YSO.lo), max(dat$bb.YSO.hi))) +
  ylab(expression(hat(beta)["YSO"])) + xlab(NULL)

## Years since outbreak - squared, lodgepole ##
dat <- parest_LP %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.YSO2) %>%
  mutate(index = row_number())

p.YSO2.LP <- ggplot(dat = dat, aes(x = index, y = bb.YSO2)) +
  geom_errorbar(aes(ymin = bb.YSO2.lo, ymax = bb.YSO2.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.YSO2.lo), max(dat$bb.YSO2.hi))) +
  ylab(expression(hat(beta)["YSO"^2])) + xlab(NULL)

## Percent dead, spruce-fir ##
dat <- parest_SF %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.pdead) %>%
  mutate(index = row_number())

p.pdead.SF <- ggplot(dat = dat, aes(x = index, y = bb.pdead)) +
  geom_errorbar(aes(ymin = bb.pdead.lo, ymax = bb.pdead.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.pdead.lo), max(dat$bb.pdead.hi))) +
  ylab(expression(hat(beta)["pdead"])) + xlab(NULL)

## Years since outbreak, spruce-fir ##
dat <- parest_SF %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.YSO) %>%
  mutate(index = row_number())

p.YSO.SF <- ggplot(dat = dat, aes(x = index, y = bb.YSO)) +
  geom_errorbar(aes(ymin = bb.YSO.lo, ymax = bb.YSO.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.YSO.lo), max(dat$bb.YSO.hi))) +
  ylab(expression(hat(beta)["YSO"])) + xlab(NULL)

## Years since outbreak - squared, spruce-fir ##
dat <- parest_SF %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.YSO2) %>%
  mutate(index = row_number())

p.YSO2.SF <- ggplot(dat = dat, aes(x = index, y = bb.YSO2)) +
  geom_errorbar(aes(ymin = bb.YSO2.lo, ymax = bb.YSO2.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.YSO2.lo), max(dat$bb.YSO2.hi))) +
  ylab(expression(hat(beta)["YSO"^2])) + xlab(NULL)

p.LP <- ggdraw() + 
  draw_plot(p.pdead.LP, x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(p.YSO.LP, x = 0.3333, y = 0, width = 0.25, height = 1) +
  draw_plot(p.YSO2.LP, x = 0.6667, y = 0, width = 0.25, height = 1)
p.SF <- ggdraw() + 
  draw_plot(p.pdead.SF, x = 0, y = 0, width = 0.3333, height = 1) +
  draw_plot(p.YSO.SF, x = 0.3333, y = 0, width = 0.3333, height = 1) +
  draw_plot(p.YSO2.SF, x = 0.6667, y = 0, width = 0.3333, height = 1)
p <- ggdraw() +
  draw_plot(p.LP, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(p.SF, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(c("Spruce-fir", "Lodgepole pine"),
                  x = c(0.44, 0.41),
                  y = c(0.48, 0.98))

save_plot("Plot_outbreak_effects_allspp.tiff", p, ncol = 2, nrow = 3, dpi = 200)


#### Plots for spp with supported outbreak effects ####

ind.supp <- which(parest_LP[, "bb.pdead.lo"] > 0 |
                    parest_LP[, "bb.pdead.hi"] < 0 |
                    parest_LP[, "bb.YSO.lo"] > 0 |
                    parest_LP[, "bb.YSO.hi"] < 0 |
                    parest_LP[, "bb.YSO2.lo"] > 0 |
                    parest_LP[, "bb.YSO2.hi"] < 0 |
                    parest_SF[, "bb.pdead.lo"] > 0 |
                    parest_SF[, "bb.pdead.hi"] < 0 |
                    parest_SF[, "bb.YSO.lo"] > 0 |
                    parest_SF[, "bb.YSO.hi"] < 0 |
                    parest_SF[, "bb.YSO2.lo"] > 0 |
                    parest_SF[, "bb.YSO2.hi"] < 0)

dat.supp <- parest_LP %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  slice(ind.supp) %>%
  mutate(stratum = "LP") %>%
  bind_rows(
    parest_SF %>% tbl_df() %>%
      mutate(Spp = spp.list) %>%
      slice(ind.supp) %>%
      mutate(stratum = "SF")
  ) %>%
  mutate(index = c(27:1, 27:1)) %>%
  mutate(bb.pdead.supp = "none") %>%
  mutate(bb.pdead.supp = replace(bb.pdead.supp, which(bb.pdead.lo > 0), "pos")) %>%
  mutate(bb.pdead.supp = replace(bb.pdead.supp, which(bb.pdead.hi < 0), "neg")) %>%
  mutate(bb.YSO.supp = "none") %>%
  mutate(bb.YSO.supp = replace(bb.YSO.supp, which(bb.YSO.lo > 0), "pos")) %>%
  mutate(bb.YSO.supp = replace(bb.YSO.supp, which(bb.YSO.hi < 0), "neg")) %>%
  mutate(bb.YSO2.supp = "none") %>%
  mutate(bb.YSO2.supp = replace(bb.YSO2.supp, which(bb.YSO2.lo > 0), "pos")) %>%
  mutate(bb.YSO2.supp = replace(bb.YSO2.supp, which(bb.YSO2.hi < 0), "neg"))

min.y <- min(c(dat.supp$bb.pdead.lo, dat.supp$bb.YSO.lo, dat.supp$bb.YSO2.lo))
max.y <- max(c(dat.supp$bb.pdead.hi, dat.supp$bb.YSO.hi, dat.supp$bb.YSO2.hi))

## Lodgepole stratum plots##
dat.plt <- dat.supp %>% filter(stratum == "LP") %>% arrange(index)

# Percent dead#
p.pdead.LP <- ggplot(dat = dat.plt, aes(x = index, y = bb.pdead, color = bb.pdead.supp)) +
  geom_errorbar(aes(ymin = bb.pdead.lo, ymax = bb.pdead.hi, color = bb.pdead.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["pdead"])) + xlab(NULL) +
  guides(color = F)

# Years since outbreak #
p.YSO.LP <- ggplot(dat = dat.plt, aes(x = index, y = bb.YSO, color = bb.YSO.supp)) +
  geom_errorbar(aes(ymin = bb.YSO.lo, ymax = bb.YSO.hi, color = bb.YSO.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["YSO"])) + xlab(NULL) +
  guides(color = F)

# Years since outbreak - squared #
p.YSO2.LP <- ggplot(dat = dat.plt, aes(x = index, y = bb.YSO2, color = bb.YSO2.supp)) +
  geom_errorbar(aes(ymin = bb.YSO2.lo, ymax = bb.YSO2.hi, color = bb.YSO2.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(beta)["YSO"^2])) + xlab(NULL) +
  guides(color = F)

## Percent dead, spruce-fir ##
dat.plt <- dat.supp %>% filter(stratum == "SF") %>% arrange(index)

# Percent dead#
p.pdead.SF <- ggplot(dat = dat.plt, aes(x = index, y = bb.pdead, color = bb.pdead.supp)) +
  geom_errorbar(aes(ymin = bb.pdead.lo, ymax = bb.pdead.hi, color = bb.pdead.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["pdead"])) + xlab(NULL) +
  guides(color = F)

# Years since outbreak #
p.YSO.SF <- ggplot(dat = dat.plt, aes(x = index, y = bb.YSO, color = bb.YSO.supp)) +
  geom_errorbar(aes(ymin = bb.YSO.lo, ymax = bb.YSO.hi, color = bb.YSO.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
  ylab(expression(hat(beta)["YSO"])) + xlab(NULL) +
  guides(color = F)

# Years since outbreak - squared #
p.YSO2.SF <- ggplot(dat = dat.plt, aes(x = index, y = bb.YSO2, color = bb.YSO2.supp)) +
  geom_errorbar(aes(ymin = bb.YSO2.lo, ymax = bb.YSO2.hi, color = bb.YSO2.supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = 1:nrow(dat.plt), labels = dat.plt$Spp, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  ylab(expression(hat(beta)["YSO"^2])) + xlab(NULL) +
  guides(color = F)

p.LP <- ggdraw() + 
  draw_plot(p.pdead.LP, x = 0, y = 0, width = 0.3333, height = 1) +
  draw_plot(p.YSO.LP, x = 0.3333, y = 0, width = 0.3333, height = 1) +
  draw_plot(p.YSO2.LP, x = 0.6667, y = 0, width = 0.3333, height = 1)
p.SF <- ggdraw() + 
  draw_plot(p.pdead.SF, x = 0, y = 0, width = 0.3333, height = 1) +
  draw_plot(p.YSO.SF, x = 0.3333, y = 0, width = 0.3333, height = 1) +
  draw_plot(p.YSO2.SF, x = 0.6667, y = 0, width = 0.3333, height = 1)
p <- ggdraw() +
  draw_plot(p.LP, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(p.SF, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(c("Spruce-fir", "Lodgepole pine"),
                  x = c(0.44, 0.41),
                  y = c(0.48, 0.98))

save_plot("Plot_outbreak_effects_SuppSpp.tiff", p, ncol = 3, nrow = 4.5, dpi = 200)
