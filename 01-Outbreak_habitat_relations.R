require(dplyr)
require(stringr)
require(ggplot2)
require(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

dat.LP <- Cov.LP %>% tbl_df %>%
  mutate(Point = row.names(Cov.LP)) %>%
  select(Point, gridIndex:Rd_dens1km) %>%
  mutate(DeadConif = DeadConif * 100)

dat.SF <- Cov.SF %>% tbl_df %>%
  mutate(Point = row.names(Cov.SF)) %>%
  select(Point, gridIndex:Rd_dens1km) %>%
  mutate(DeadConif = DeadConif * 100)

#### Scatterplots - DeadConif VS YSO ####
p.LP <- ggplot(data = dat.LP, aes(x = YSO, y = DeadConif)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 105, label = "Lodgpole pine stratum", size = 4)

p.SF <- ggplot(data = dat.SF, aes(x = YSO, y = DeadConif)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 105, label = "Spruce-fir stratum", size = 4)

p <- ggdraw() + 
  draw_plot(p.LP, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p.SF, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Percent conifer dead", "Years since outbreak"),
                  x = c(0.03, 0.4),
                  y = c(0.25, 0.03),
                  size = c(15, 15),
                  angle = c(90, 0),
                  hjust = c(0, 0),
                  vjust = c(0, 0))

save_plot("figure_DeadConif_VS_YSO.tiff", p, ncol = 2, nrow = 1, dpi = 200)

## Compare DeadConif between outbreak vs non-outbreak grids ##
#dat.LP %>%
#  group_by(gridIndex) %>%
#  summarise(DConif = mean(DeadConif, na.rm = T), DConif_max = max(DeadConif, na.rm = T), YSI = mean(YSO, na.rm = T)) %>%
#  ungroup %>%
#  mutate(Outbreak = !is.na(YSI)) %>%
#  group_by(Outbreak) %>%
#  summarise(DConif_mn = mean(DConif, na.rm = T), DConif_SD = sd(DConif, na.rm = T), DConif_max = max(DConif, na.rm = T)) %>%
#  View

#dat.SF %>%
#  group_by(gridIndex) %>%
#  summarise(DConif = mean(DeadConif, na.rm = T), DConif_max = max(DeadConif, na.rm = T), YSI = mean(YSO, na.rm = T)) %>%
#  ungroup %>%
#  mutate(Outbreak = !is.na(YSI)) %>%
#  group_by(Outbreak) %>%
#  summarise(DConif_mn = mean(DConif, na.rm = T), DConif_SD = sd(DConif, na.rm = T), DConif_max = max(DConif, na.rm = T)) %>%
#  View

#### Regression of vegetation (mechanistic) factors VS outbreak metrics ####
vars <- c("CanCov", "RCOV_AS", "RCOV_ES", "RCOV_Pine", "shrub_cover", "RCShrb_UC", "HerbCov", "WoodyCov", "DDCov")
cols <- c("Int.LP", "YSO.LP", "YSO2.LP", "DCon.LP", "YSOXDCon.LP", "Int.SF", "YSO.SF", "YSO2.SF", "DCon.SF", "YSOXDCon.SF")
out <- matrix(NA, nrow = length(vars), ncol = length(cols),
              dimnames = list(vars, cols))

for(i in 1:length(vars)) {
  mod <- lm(as.formula(str_c(vars[i], "~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO)")), data = dat.LP)
  sum.tab <- summary(mod)$coefficients
  out[vars[i], c("Int.LP", "YSO.LP", "YSO2.LP", "DCon.LP", "YSOXDCon.LP")] <-
    str_c(round(sum.tab[, "Estimate"], digits = 3), " (",
          round(sum.tab[, "Std. Error"], digits = 3), ")")
  if(any(sum.tab[-1, "Pr(>|t|)"] < 0.05)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] < 0.05)
    out[vars[i], c("YSO.LP", "YSO2.LP", "DCon.LP", "YSOXDCon.LP")][sig.ind] <-
      str_c(out[vars[i], c("YSO.LP", "YSO2.LP", "DCon.LP", "YSOXDCon.LP")][sig.ind], "*")
  }
  if(any(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)
    out[vars[i], c("YSO.LP", "YSO2.LP", "DCon.LP", "YSOXDCon.LP")][sig.ind] <-
      str_c(out[vars[i], c("YSO.LP", "YSO2.LP", "DCon.LP", "YSOXDCon.LP")][sig.ind], ".")
  }
  mod <- lm(as.formula(str_c(vars[i], "~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO)")), data = dat.SF)
  sum.tab <- summary(mod)$coefficients
  out[vars[i], c("Int.SF", "YSO.SF", "YSO2.SF", "DCon.SF", "YSOXDCon.SF")] <-
    str_c(round(sum.tab[, "Estimate"], digits = 3), " (",
          round(sum.tab[, "Std. Error"], digits = 3), ")")
  if(any(sum.tab[-1, "Pr(>|t|)"] < 0.05)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] < 0.05)
    out[vars[i], c("YSO.SF", "YSO2.SF", "DCon.SF", "YSOXDCon.SF")][sig.ind] <-
      str_c(out[vars[i], c("YSO.SF", "YSO2.SF", "DCon.SF", "YSOXDCon.SF")][sig.ind], "*")
  }
  if(any(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)
    out[vars[i], c("YSO.SF", "YSO2.SF", "DCon.SF", "YSOXDCon.SF")][sig.ind] <-
      str_c(out[vars[i], c("YSO.SF", "YSO2.SF", "DCon.SF", "YSOXDCon.SF")][sig.ind], ".")
  }
}

write.csv(out, "Outbreak_habitat_relations.csv", row.names = T)

#### Plot vegetation (mehanistic factors) VS outbreak metrics ####
## Canopy relationships ##
dat.obs <- dat.LP
mod <- lm(CanCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC.LP <- ggplot(dat.obs, aes(x = YSO, y = CanCov)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(0, 60) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon.LP <- ggplot(dat.obs, aes(x = DeadConif, y = CanCov)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 60) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.obs <- dat.SF
mod <- lm(CanCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC.SF <- ggplot(dat.obs, aes(x = YSO, y = CanCov)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(0, 60) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon.SF <- ggplot(dat.obs, aes(x = DeadConif, y = CanCov)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 60) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

p.CCov <- ggdraw() + 
  draw_plot(p.CC.LP, x = 0, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon.LP, x = .25, y = 0, width = .25, height = 1) +
  draw_plot(p.CC.SF, x = 0.5, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon.SF, x = .75, y = 0, width = .25, height = 1)

## Aspen relationships ##
dat.obs <- dat.LP
mod <- lm(RCOV_AS ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC <- ggplot(dat.obs, aes(x = YSO, y = RCOV_AS)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(c(0, 100)) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon <- ggplot(dat.obs, aes(x = DeadConif, y = RCOV_AS)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                       max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in unique(dat.pred$YSO)) {
  maxdc <- max(dat.obs$DeadConif[which(dat.obs$YSO == i)], na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))

p.heat <- ggplot(dat.pred, aes(YSO, DeadConif)) +
  geom_tile(aes(fill = Y), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Dead conifer (%)") +
  xlab("Years since outbreak") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(fill = "% Aspen")

p.AS.LP <- ggdraw() +
  draw_plot(p.CC, x = 0, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon, x = .25, y = 0, width = .25, height = 1) +
  draw_plot(p.heat, x = 0.5, y = 0, width = .5, height = 1)

dat.obs <- dat.SF
mod <- lm(RCOV_AS ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC <- ggplot(dat.obs, aes(x = YSO, y = RCOV_AS)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(c(0, 100)) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon <- ggplot(dat.obs, aes(x = DeadConif, y = RCOV_AS)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #heme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in unique(dat.pred$YSO)) {
  maxdc <- max(dat.obs$DeadConif[which(dat.obs$YSO == i)], na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))

p.heat <- ggplot(dat.pred, aes(YSO, DeadConif)) +
  geom_tile(aes(fill = Y), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Dead conifer (%)") +
  xlab("Years since outbreak") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(fill = "% Aspen")

p.AS.SF <- ggdraw() +
  draw_plot(p.CC, x = 0, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon, x = .25, y = 0, width = .25, height = 1) +
  draw_plot(p.heat, x = 0.5, y = 0, width = .5, height = 1)

## Spruce dominance relationships (SF stratum only) ##
dat.obs <- dat.SF
mod <- lm(RCOV_ES ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC <- ggplot(dat.obs, aes(x = YSO, y = RCOV_ES)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(c(0, 100)) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon <- ggplot(dat.obs, aes(x = DeadConif, y = RCOV_ES)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in unique(dat.pred$YSO)) {
  maxdc <- max(dat.obs$DeadConif[which(dat.obs$YSO == i)], na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))

p.heat <- ggplot(dat.pred, aes(YSO, DeadConif)) +
  geom_tile(aes(fill = Y), color = "white") +
  scale_fill_gradient(low = "green", high = "red") +
  ylab("Dead conifer (%)") +
  xlab("Years since outbreak") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(fill = "% Spruce")

p.ES.SF <- ggdraw() +
  draw_plot(p.CC, x = 0, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon, x = .25, y = 0, width = .25, height = 1) +
  draw_plot(p.heat, x = 0.5, y = 0, width = .5, height = 1)

## Pine dominance relationships ##
dat.obs <- dat.LP
mod <- lm(RCOV_Pine ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC.LP <- ggplot(dat.obs, aes(x = YSO, y = RCOV_Pine)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon.LP <- ggplot(dat.obs, aes(x = DeadConif, y = RCOV_Pine)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.obs <- dat.SF
mod <- lm(RCOV_Pine ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.CC.SF <- ggplot(dat.obs, aes(x = YSO, y = RCOV_Pine)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = YSO, y = fit)) +
  geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Years since outbreak") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))

p.DCon.SF <- ggplot(dat.obs, aes(x = DeadConif, y = RCOV_Pine)) +
  geom_point(alpha = 0.1) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = fit)) +
  geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed") +
  geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed") +
  ylim(0, 100) +
  xlab("Dead conifer (%)") + ylab(NULL)# +
  #theme(axis.title.x = element_text(size = 20)) +
  #theme(axis.title.y = element_text(size = 20)) +
  #theme(axis.text.x = element_text(size = 15)) +
  #theme(axis.text.y = element_text(size = 15))

p.Pine <- ggdraw() + 
  draw_plot(p.CC.LP, x = 0, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon.LP, x = .25, y = 0, width = .25, height = 1) +
  draw_plot(p.CC.SF, x = 0.5, y = 0, width = .25, height = 1) +
  draw_plot(p.DCon.SF, x = .75, y = 0, width = .25, height = 1)

p.Canopy <- ggdraw() + 
  draw_plot(p.CCov, x = 0.05, y = 0.78, width = .95, height = 0.17) +
  draw_plot(p.Pine, x = .05, y = 0.61, width = .95, height = 0.17) +
  draw_plot(p.AS.LP, x = .05, y = 0.39, width = .95, height = 0.17) +
  draw_plot(p.AS.SF, x = 0.05, y = 0.17, width = .95, height = 0.17) +
  draw_plot(p.ES.SF, x = 0.05, y = 0, width = .95, height = 0.17) +
  draw_plot_label(c("% Spruce", "% Aspen", "% Pine", "Canopy cover",
                    "Lodgepole pine", "Spruce-fir", "Lodgepole pine", "Spruce-fir"),
                  x = c(0.03, 0.03, 0.03, 0.03, 0.25, 0.7, 0.48, 0.5),
                  y = c(0.065, 0.35, 0.67, 0.8, 0.97, 0.97, 0.58, 0.36),
                  size = c(15, 15, 15, 15, 17, 17, 17, 17),
                  angle = c(90, 90, 90, 90, 0, 0, 0, 0),
                  hjust = c(0, 0, 0, 0, 0, 0, 0, 0),
                  vjust = c(0, 0, 0, 0, 0, 0, 0, 0))
  
save_plot("figure_Canopy_VS_outbreak.tiff", p.Canopy, ncol = 2, nrow = 2.5, dpi = 300)
