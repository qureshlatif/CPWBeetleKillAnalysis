library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

mod_LP <- loadObject("mod_LPcommunity_outbreak_reduced2")
mod_SF <- loadObject("mod_SFcommunity_outbreak_reduced2")

# Tabulate parameter estimates
pars <- c("bb.pdead", "bb.YSO", "bb.YSO2", "bb.pdXYSO")
cols <- (c("", ".lo", ".hi") %>%
           expand.grid(pars, stringsAsFactors = F) %>%
           select(Var2, Var1) %>%
           mutate(Var3 = str_c(Var2, Var1, sep = "")))$Var3
parest_LP <- parest_SF <- matrix(NA, nrow = length(spp.list), ncol = length(cols), dimnames = list(NULL, cols))

for(i in 1:length(pars)) {
  parm <- mod_LP$sims.list[[pars[i]]]
  parest_LP[, pars[i]] <- apply(parm, 2, median)
  parest_LP[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  parest_LP[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
}
for(i in 1:length(pars)) {
  parm <- mod_SF$sims.list[[pars[i]]]
  parest_SF[, pars[i]] <- apply(parm, 2, median)
  parest_SF[, str_c(pars[i], ".lo")] <- apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8))
  parest_SF[, str_c(pars[i], ".hi")] <- apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8))
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

## Percent dead X Years since outbreak, lodgepole ##
dat <- parest_LP %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.pdXYSO) %>%
  mutate(index = row_number())

p.pdXYSO.LP <- ggplot(dat = dat, aes(x = index, y = bb.pdXYSO)) +
  geom_errorbar(aes(ymin = bb.pdXYSO.lo, ymax = bb.pdXYSO.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.pdXYSO.lo), max(dat$bb.pdXYSO.hi))) +
  ylab(expression(hat(beta)["pdead X YSO"])) + xlab(NULL)

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

## Percent dead X Years since outbreak, spruce-fir ##
dat <- parest_SF %>% tbl_df() %>%
  mutate(Spp = spp.list) %>%
  arrange(bb.pdXYSO) %>%
  mutate(index = row_number())

p.pdXYSO.SF <- ggplot(dat = dat, aes(x = index, y = bb.pdXYSO)) +
  geom_errorbar(aes(ymin = bb.pdXYSO.lo, ymax = bb.pdXYSO.hi), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min(dat$bb.pdXYSO.lo), max(dat$bb.pdXYSO.hi))) +
  ylab(expression(hat(beta)["pdead X YSO"])) + xlab(NULL)

p.LP <- ggdraw() + 
  draw_plot(p.pdead.LP, x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(p.YSO.LP, x = 0.25, y = 0, width = 0.25, height = 1) +
  draw_plot(p.YSO2.LP, x = 0.5, y = 0, width = 0.25, height = 1) +
  draw_plot(p.pdXYSO.LP, x = 0.75, y = 0, width = 0.25, height = 1)
p.SF <- ggdraw() + 
  draw_plot(p.pdead.SF, x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(p.YSO.SF, x = 0.25, y = 0, width = 0.25, height = 1) +
  draw_plot(p.YSO2.SF, x = 0.5, y = 0, width = 0.25, height = 1) +
  draw_plot(p.pdXYSO.SF, x = 0.75, y = 0, width = 0.25, height = 1)
p <- ggdraw() +
  draw_plot(p.LP, x = 0, y = 0.5, width = 1, height = 0.45) +
  draw_plot(p.SF, x = 0, y = 0, width = 1, height = 0.45) +
  draw_plot_label(c("Spruce-fir", "Lodgepole pine"),
                  x = c(0.44, 0.41),
                  y = c(0.48, 0.98))

save_plot("Plot_outbreak_effects_allspp.tiff", p, ncol = 3, nrow = 3, dpi = 200)


#### Plots for spp with supported outbreak effects ####

# ind.supp <- which(parest_LP[, "bb.pdead.lo"] > 0 |
#                     parest_LP[, "bb.pdead.hi"] < 0 |
#                     parest_LP[, "bb.YSO.lo"] > 0 |
#                     parest_LP[, "bb.YSO.hi"] < 0 |
#                     parest_LP[, "bb.YSO2.lo"] > 0 |
#                     parest_LP[, "bb.YSO2.hi"] < 0 |
#                     parest_LP[, "bb.pdXYSO.lo"] > 0 |
#                     parest_LP[, "bb.pdXYSO.hi"] < 0 |
#                     parest_SF[, "bb.pdead.lo"] > 0 |
#                     parest_SF[, "bb.pdead.hi"] < 0 |
#                     parest_SF[, "bb.YSO.lo"] > 0 |
#                     parest_SF[, "bb.YSO.hi"] < 0 |
#                     parest_SF[, "bb.YSO2.lo"] > 0 |
#                     parest_SF[, "bb.YSO2.hi"] < 0 |
#                     parest_SF[, "bb.pdXYSO.lo"] > 0 |
#                     parest_SF[, "bb.pdXYSO.hi"] < 0 )
# 
# dat.supp <- parest_LP %>% tbl_df() %>%
#   mutate(Spp = spp.list) %>%
#   slice(ind.supp) %>%
#   mutate(stratum = "LP") %>%
#   bind_rows(
#     parest_SF %>% tbl_df() %>%
#       mutate(Spp = spp.list) %>%
#       slice(ind.supp) %>%
#       mutate(stratum = "SF")
#   ) %>%
#   mutate(index = c(16:1, 16:1)) %>%
#   mutate(bb.pdead.supp = "none") %>%
#   mutate(bb.pdead.supp = replace(bb.pdead.supp, which(bb.pdead.lo > 0), "pos")) %>%
#   mutate(bb.pdead.supp = replace(bb.pdead.supp, which(bb.pdead.hi < 0), "neg")) %>%
#   mutate(bb.YSO.supp = "none") %>%
#   mutate(bb.YSO.supp = replace(bb.YSO.supp, which(bb.YSO.lo > 0), "pos")) %>%
#   mutate(bb.YSO.supp = replace(bb.YSO.supp, which(bb.YSO.hi < 0), "neg")) %>%
#   mutate(bb.YSO2.supp = "none") %>%
#   mutate(bb.YSO2.supp = replace(bb.YSO2.supp, which(bb.YSO2.lo > 0), "pos")) %>%
#   mutate(bb.YSO2.supp = replace(bb.YSO2.supp, which(bb.YSO2.hi < 0), "neg")) %>%
#   mutate(bb.pdXYSO.supp = "none") %>%
#   mutate(bb.pdXYSO.supp = replace(bb.pdXYSO.supp, which(bb.pdXYSO.lo > 0), "pos")) %>%
#   mutate(bb.pdXYSO.supp = replace(bb.pdXYSO.supp, which(bb.pdXYSO.hi < 0), "neg"))
#   
# p.pdead <- ggplot(dat = dat.supp, aes(x = index, y = bd.ptrt)) +
#   geom_errorbar(aes(ymin = bd.ptrt.lo, ymax = bd.ptrt.hi), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bd.ptrt.lo), max(dat.plt$bd.ptrt.hi))) +
#   ylab(NULL) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15))
# 
# pd.ptrt.supp <- ggplot(dat = dat.supp, aes(x = index, y = bd.ptrt)) +
#   geom_errorbar(aes(ymin = bd.ptrt.lo, ymax = bd.ptrt.hi, color = bd.ptrt.supp), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bd.ptrt.lo), max(dat.plt$bd.ptrt.hi))) +
#   scale_color_manual(values = c("#000000", "#D55E00")) +
#   ylab(expression(hat(beta)["percent treated"])) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15)) +
#   guides(color = F)
# 
# pd.YST <- ggplot(dat = dat.plt, aes(x = index, y = bd.YST)) +
#   geom_errorbar(aes(ymin = bd.YST.lo, ymax = bd.YST.hi), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bd.YST.lo), max(dat.plt$bd.YST.hi))) +
#   ylab(NULL) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15))
# 
# pd.YST.supp <- ggplot(dat = dat.supp, aes(x = index, y = bd.YST)) +
#   geom_errorbar(aes(ymin = bd.YST.lo, ymax = bd.YST.hi), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bd.YST.lo), max(dat.plt$bd.YST.hi))) +
#   ylab(expression(hat(beta)["mean years since treatment"])) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15)) +
#   guides(color = F)
# 
# pb.trt <- ggplot(dat = dat.plt, aes(x = index, y = bb.trt)) +
#   geom_errorbar(aes(ymin = bb.trt.lo, ymax = bb.trt.hi), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bb.trt.lo), max(dat.plt$bb.trt.hi))) +
#   ylab(NULL) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15))
# 
# pb.trt.supp <- ggplot(dat = dat.supp, aes(x = index, y = bb.trt)) +
#   geom_errorbar(aes(ymin = bb.trt.lo, ymax = bb.trt.hi, color = bb.trt.supp), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bb.trt.lo), max(dat.plt$bb.trt.hi))) +
#   scale_color_manual(values = c("#0072B2", "#000000", "#D55E00")) +
#   ylab(expression(hat(beta)["treatment status"])) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15)) +
#   guides(color = F)
# 
# pb.YST <- ggplot(dat = dat.plt, aes(x = index, y = bb.YST)) +
#   geom_errorbar(aes(ymin = bb.YST.lo, ymax = bb.YST.hi), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bb.YST.lo), max(dat.plt$bb.YST.hi))) +
#   ylab(NULL) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15))
# 
# pb.YST.supp <- ggplot(dat = dat.supp, aes(x = index, y = bb.YST)) +
#   geom_errorbar(aes(ymin = bb.YST.lo, ymax = bb.YST.hi, color = bb.YST.supp), size=1, width=0) +
#   geom_point(size = 2.5) + 
#   geom_hline(yintercept = 0) +
#   coord_flip() +
#   scale_x_continuous(breaks = 1:nrow(dat.supp), labels = dat.supp$Spp, expand=c(0, 1)) +
#   scale_y_continuous(lim = c(min(dat.plt$bb.YST.lo), max(dat.plt$bb.YST.hi))) +
#   scale_color_manual(values = c("#0072B2", "#000000")) +
#   ylab(expression(hat(beta)["years since treatment"])) + xlab(NULL) +
#   theme(axis.title.y=element_text(size=30)) +
#   theme(axis.title.x=element_text(size=30)) +
#   theme(axis.text.x=element_text(size=15)) +
#   theme(axis.text.y=element_text(size=15)) +
#   guides(color = F)
# 
# p <- ggdraw() + 
#   draw_plot(pd.ptrt, x = 0.05, y = 0.6333333, width = 0.2375, height = 0.3166667) +
#   draw_plot(pd.ptrt.supp, x = 0.05, y = 0, width = 0.2375, height = 0.6333333) +
#   draw_plot(pd.YST, x = 0.2875, y = 0.6333333, width = 0.2375, height = 0.3166667) +
#   draw_plot(pd.YST.supp, x = 0.2875, y = 0, width = 0.2375, height = 0.6333333) +
#   draw_plot(pb.trt, x = 0.525, y = 0.6333333, width = 0.2375, height = 0.3166667) +
#   draw_plot(pb.trt.supp, x = 0.525, y = 0, width = 0.2375, height = 0.6333333) +
#   draw_plot(pb.YST, x = 0.7625, y = 0.6333333, width = 0.2375, height = 0.3166667) +
#   draw_plot(pb.YST.supp, x = 0.7625, y = 0, width = 0.2375, height = 0.6333333) +
#   draw_plot_label(c("Species", "Grid scale", "Point scale"),
#                   x = c(0, 0.25, 0.7),
#                   y = c(0.5, 0.98, 0.98),
#                   size = c(40, 30, 30),
#                   angle = c(90, 0, 0),
#                   hjust = c(0, 0, 0))
# 
# save_plot("Plot_trt_effects.tiff", p, ncol = 4, nrow = 4.5, dpi = 200)
