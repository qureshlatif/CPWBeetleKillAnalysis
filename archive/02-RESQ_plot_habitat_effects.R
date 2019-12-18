library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)
library(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/CPW/")
load("Data_compiled_RESQ.RData")

mod_LP <- loadObject("mod_RESQ_habitat_HZdist_LP")
mod_SF <- loadObject("mod_RESQ_habitat_HZdist_SF")
samp.ha <- sum(area.band) * 0.0001

# Tabulate slope parameter estimates
pars <- c("bl.ccov", "bl.RCovAS", "bl.RCovPine", "bl.RCovES", "bl.shcov", "bl.RSC_Con", "bl.GHerb", "bl.Gwoody", "bl.GDD")
cols <- c("md", "lo", "hi")
parest_LP <- parest_SF <- matrix(NA, nrow = length(pars), ncol = length(cols), dimnames = list(pars, cols))

for(par in pars[-4]) {
  parm <- mod_LP$sims.list[[par]]
  parest_LP[par, "md"] <- median(parm)
  parest_LP[par, "lo"] <- quantile(parm, prob = 0.05, type = 8)
  parest_LP[par, "hi"] <- quantile(parm, prob = 0.95, type = 8)
}
dat.pars.LP <- parest_LP %>% tbl_df() %>%
  mutate(vars = c("CanCov", "Aspen", "Pine", "Spruce", "ShrubCov",
                  "ConShrb", "Herb", "Woody", "DeadDown"),
         x = length(pars):1,
         supp = "none") %>%
  mutate(supp = replace(supp, which(lo > 0), "pos")) %>%
  mutate(supp = replace(supp, which(hi < 0), "neg")) %>%
  mutate(supp = as.factor(supp))
  
for(par in pars) {
  parm <- mod_SF$sims.list[[par]]
  parest_SF[par, "md"] <- median(parm)
  parest_SF[par, "lo"] <- quantile(parm, prob = 0.05, type = 8)
  parest_SF[par, "hi"] <- quantile(parm, prob = 0.95, type = 8)
}
dat.pars.SF <- parest_SF %>% tbl_df() %>%
  mutate(vars = c("CanCov", "Aspen", "Pine", "Spruce", "ShrubCov",
                  "ConShrb", "Herb", "Woody", "DeadDown"),
         x = length(pars):1,
         supp = "none") %>%
  mutate(supp = replace(supp, which(lo > 0), "pos")) %>%
  mutate(supp = replace(supp, which(hi < 0), "neg")) %>%
  mutate(supp = as.factor(supp))
rm(parm)

min.y <- min(c(dat.pars.LP$lo, dat.pars.SF$lo))
max.y <- max(c(dat.pars.LP$hi, dat.pars.SF$hi))

#### Plots ####
p.LP <- ggplot(dat = dat.pars.LP, aes(x = x, y = md, color = supp)) +
  geom_errorbar(aes(ymin = lo, ymax = hi, color = supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = dat.pars.LP$x, labels = dat.pars.LP$vars, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta))) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p.SF <- ggplot(dat = dat.pars.SF, aes(x = x, y = md, color = supp)) +
  geom_errorbar(aes(ymin = lo, ymax = hi, color = supp), size=1, width=0) +
  geom_point(size = 2.5) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_x_continuous(breaks = dat.pars.SF$x, labels = dat.pars.SF$vars, expand=c(0, 1)) +
  scale_y_continuous(lim = c(min.y, max.y)) +
  scale_color_manual(values = c("#0072B2", "dark gray", "#D55E00")) +
  ylab(expression(hat(beta))) + xlab(NULL) +
  theme(axis.title.x=element_text(size=25)) +
  guides(color = F)

p <- ggdraw() + 
  draw_plot(p.LP, x = 0, y = 0, width = 0.5, height = 0.95) +
  draw_plot(p.SF, x = 0.5, y = 0, width = 0.5, height = 0.95) +
  draw_plot_label(c("Lodgepole", "Spruce-fir"),
                  x = c(0.25, 0.7), y = c(0.98, 0.98))

save_plot("Plot_habitat_effects_RESQ.tiff", p, ncol = 2, nrow = 2, dpi = 200)

