require(dplyr)
require(stringr)
require(ggplot2)
require(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

dat.LP <- Cov.LP %>% tbl_df %>%
  mutate(Point = row.names(Cov.LP)) %>%
  select(Point, gridIndex:Rd_dens1km) %>%
  mutate(DeadConif = DeadConif * 100) %>%
  left_join(Cov.LP %>% tbl_df %>%
              group_by(gridIndex) %>%
              summarise(Outbreak = any(!is.na(YSO)) %>% as.integer,
                        YSO.med = median(YSO, na.rm = T) %>% round), by = "gridIndex") %>%
  mutate(YSO = replace(YSO, which(is.na(YSO) & Outbreak == 0), 0))
ind <- which(is.na(dat.LP$YSO) & dat.LP$Outbreak == 1)
dat.LP$YSO[ind] <- dat.LP$YSO.med[ind]
dat.LP <- dat.LP %>% select(-YSO.med)

dat.SF <- Cov.SF %>% tbl_df %>%
  mutate(Point = row.names(Cov.SF)) %>%
  select(Point, gridIndex:Rd_dens1km) %>%
  mutate(DeadConif = DeadConif * 100) %>%
  left_join(Cov.SF %>% tbl_df %>%
              group_by(gridIndex) %>%
              summarise(Outbreak = any(!is.na(YSO)) %>% as.integer,
                        YSO.med = median(YSO, na.rm = T) %>% round), by = "gridIndex") %>%
  mutate(YSO = replace(YSO, which(is.na(YSO) & Outbreak == 0), 0))
ind <- which(is.na(dat.SF$YSO) & dat.SF$Outbreak == 1)
dat.SF$YSO[ind] <- dat.SF$YSO.med[ind]
dat.SF <- dat.SF %>% select(-YSO.med)

#### Scatterplots - DeadConif VS YSO ####
p.LP <- ggplot(data = dat.LP %>%
                 mutate(YSO = replace(YSO, which(Outbreak == 0), -2)), aes(x = YSO, y = DeadConif)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 105, label = "Lodgpole pine stratum", size = 4) +
  geom_vline(xintercept = -1, linetype = "dashed") +
  annotate("text", x = -2, y = 90, label = "No outbreak", size = 3, angle = 90)

p.SF <- ggplot(data = dat.SF %>%
                 mutate(YSO = replace(YSO, which(Outbreak == 0), -2)), aes(x = YSO, y = DeadConif)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 105, label = "Spruce-fir stratum", size = 4) +
  geom_vline(xintercept = -1, linetype = "dashed") +
  annotate("text", x = -2, y = 90, label = "No outbreak", size = 3, angle = 90)

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

# Adjust data following review of this figure
dat.LP <- dat.LP %>%
  mutate(YSO = replace(YSO, which(Outbreak == 0), -1))

dat.SF <- dat.SF %>%
  mutate(YSO = replace(YSO, which(Outbreak == 0), -1))

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

# Rearrange #
out <- rbind(out[, 1:5], out[, 6:10]) %>%
  cbind(Stratum = c(rep("LP", 9), rep("SF", 9)), .) %>%
  cbind(HabFeature = dimnames(out)[[1]], .) %>%
  tbl_df %>%
  slice((1:18) %>% matrix(ncol = 2) %>% t %>% as.integer)

write.csv(out, "Outbreak_habitat_relations.csv", row.names = T)

#### Plot vegetation (mehanistic factors) VS outbreak metrics ####
#___________ Plotting functions ____________#
p.YSO <- function(dat.obs, dat.pred, ymax) {
  p <- ggplot(dat.obs, aes(x = YSO, y = habvar)) +
    geom_point(alpha = 0.1) +
    geom_line(data = dat.pred, aes(x = YSO, y = fit), color = "red") +
    geom_line(data = dat.pred, aes(x = YSO, y = lwr), linetype = "dashed", color = "red") +
    geom_line(data = dat.pred, aes(x = YSO, y = upr), linetype = "dashed", color = "red") +
    ylim(0, ymax) +
    xlab("Years since outbreak") + ylab(NULL)
  return(p)
}

p.DCON <- function(dat.obs, dat.pred, ymax) {
  p <- ggplot(dat.obs, aes(x = DeadConif, y = habvar)) +
    geom_point(alpha = 0.1) +
    geom_line(data = dat.pred, aes(x = DeadConif, y = fit), color = "red") +
    geom_line(data = dat.pred, aes(x = DeadConif, y = lwr), linetype = "dashed", color = "red") +
    geom_line(data = dat.pred, aes(x = DeadConif, y = upr), linetype = "dashed", color = "red") +
    ylim(0, ymax) +
    xlab("Dead conifer (%)") + ylab(NULL)
  return(p)
}

p.HEAT <- function(dat.pred, heat.legend.title) {
  p <- ggplot(dat.pred, aes(YSO, DeadConif)) +
    geom_tile(aes(fill = Y), color = "white") +
    scale_fill_gradient(low = "green", high = "red") +
    ylab("Dead conifer (%)") +
    xlab("Years since outbreak") +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 12)) +
    labs(fill = heat.legend.title)
  return(p)
}

p.stitch.row <- function(p.CC, p.DCon, p.heat) {
  p <- ggdraw() +
    draw_plot(p.CC, x = 0, y = 0, width = .25, height = 1) +
    draw_plot(p.DCon, x = .25, y = 0, width = .25, height = 1) +
    draw_plot(p.heat, x = 0.5, y = 0, width = .5, height = 1)
  return(p)
}
#__________________________________#

## Canopy layer relationships ##
# Canopy cover #
  # LP #
dat.obs <- dat.LP
mod <- lm(CanCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = CanCov), dat.pred, ymax = 60)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = CanCov), dat.pred, ymax = 60)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "Canopy cover")

p.CC.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

  # SF #
dat.obs <- dat.SF
mod <- lm(CanCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = CanCov), dat.pred, ymax = 60)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = CanCov), dat.pred, ymax = 60)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "Canopy cover")

p.CC.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

# Aspen relationships #
# LP #
dat.obs <- dat.LP
mod <- lm(RCOV_AS ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCOV_AS), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCOV_AS), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Aspen")

p.AS.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(RCOV_AS ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCOV_AS), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCOV_AS), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Aspen")

p.AS.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

# Spruce dominance relationships (SF stratum only) #
# SF #
dat.obs <- dat.SF
mod <- lm(RCOV_ES ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCOV_ES), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCOV_ES), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Spruce")

p.ES.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

# Pine dominance relationships #
# LP #
dat.obs <- dat.LP
mod <- lm(RCOV_Pine ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCOV_Pine), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCOV_Pine), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Pine")

p.Pine.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(RCOV_Pine ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCOV_Pine), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCOV_Pine), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Pine")

p.Pine.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

p.Canopy <- ggdraw() + 
  draw_plot(p.CC.LP, x = .05, y = 0.7125, width = .475, height = 0.2375) +
  draw_plot(p.CC.SF, x = .525, y = 0.7125, width = .475, height = 0.2375) +
  draw_plot(p.AS.LP, x = .05, y = 0.475, width = .475, height = 0.2375) +
  draw_plot(p.AS.SF, x = .525, y = 0.475, width = .475, height = 0.2375) +
  draw_plot(p.ES.SF, x = .525, y = 0.2375, width = .475, height = 0.2375) +
  draw_plot(p.Pine.LP, x = .05, y = 0, width = .475, height = 0.2375) +
  draw_plot(p.Pine.SF, x = .525, y = 0, width = .475, height = 0.2375) +
  draw_plot_label(c("% Pine", "% Spruce", "% Aspen", "Canopy cover",
                    "Lodgepole pine", "Spruce-fir"),
                  x = c(.03, .52, .03, .03, .25, .7),
                  y = c(.1, .33, .58, .8, .97, .97),
                  size = c(15, 15, 15, 15, 17, 17),
                  angle = c(90, 90, 90, 90, 0, 0),
                  hjust = c(0, 0, 0, 0, 0, 0),
                  vjust = c(0, 0, 0, 0, 0, 0))


save_plot("figure_Canopy_VS_outbreak.tiff", p.Canopy, ncol = 4, nrow = 2.5, dpi = 300)

## Shrub layer relationships ##
# Shrub cover #
# LP #
dat.obs <- dat.LP
mod <- lm(shrub_cover ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = shrub_cover), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = shrub_cover), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "Shrub cover")

p.Shrub.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(shrub_cover ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = shrub_cover), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = shrub_cover), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "Shrub cover")

p.Shrub.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

# Conifer shrubs #
# LP #
dat.obs <- dat.LP
mod <- lm(RCShrb_UC ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCShrb_UC), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCShrb_UC), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Conifer")

p.CShrb.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(RCShrb_UC ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = RCShrb_UC), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = RCShrb_UC), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Conifer")

p.CShrb.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

p.Shrub <- ggdraw() + 
  draw_plot(p.Shrub.LP, x = .05, y = 0.475, width = .475, height = 0.475) +
  draw_plot(p.Shrub.SF, x = .525, y = 0.475, width = .475, height = 0.475) +
  draw_plot(p.CShrb.LP, x = .05, y = 0, width = .475, height = 0.475) +
  draw_plot(p.CShrb.SF, x = .525, y = 0, width = .475, height = 0.475) +
  draw_plot_label(c("% Conifer shrubs", "Shrub Cover",
                    "Lodgepole pine", "Spruce-fir"),
                  x = c(.03, .03, .25, .7),
                  y = c(.25, .75, .97, .97),
                  size = c(15, 15, 17, 17),
                  angle = c(90, 90, 0, 0),
                  hjust = c(0, 0, 0, 0),
                  vjust = c(0, 0, 0, 0))

save_plot("figure_Shrub_VS_outbreak.tiff", p.Shrub, ncol = 4, nrow = 2, dpi = 300)

## Ground cover relationships ##
# Woody stems #
# LP #
dat.obs <- dat.LP
mod <- lm(WoodyCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = WoodyCov), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = WoodyCov), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Woody")

p.Woody.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(WoodyCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = WoodyCov), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = WoodyCov), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Woody")

p.Woody.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

# Dead and down #
# LP #
dat.obs <- dat.LP
mod <- lm(DDCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = DDCov), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = DDCov), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Dead and down")

p.DDCov.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(DDCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = DDCov), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = DDCov), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Dead and down")

p.DDCov.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

# Herbaceous #
# LP #
dat.obs <- dat.LP
mod <- lm(HerbCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:6

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = HerbCov), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = HerbCov), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Herbaceous")

p.Herb.LP <- p.stitch.row(p.CC, p.DCon, p.heat)

# SF #
dat.obs <- dat.SF
mod <- lm(HerbCov ~ YSO + I(YSO^2) + DeadConif + I(DeadConif*YSO), data = dat.obs)
YSO.limited <- -1:3

dat.pred <- data.frame(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                       DeadConif = mean(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.CC <- p.YSO(dat.obs %>% rename(habvar = HerbCov), dat.pred, ymax = 100)

dat.pred <- data.frame(YSO = mean(dat.obs$YSO, na.rm = T),
                       DeadConif = min(dat.obs$DeadConif, na.rm = T):
                         max(dat.obs$DeadConif, na.rm = T))
dat.pred <- cbind(dat.pred, predict.lm(mod, dat.pred, interval = "confidence"))
p.DCon <- p.DCON(dat.obs %>% rename(habvar = HerbCov), dat.pred, ymax = 100)

dat.pred <- expand.grid(YSO = min(dat.obs$YSO, na.rm = T):max(dat.obs$YSO, na.rm = T),
                        DeadConif = seq(min(dat.obs$DeadConif, na.rm = T),
                                        max(dat.obs$DeadConif, na.rm = T), length.out = 10))
for(i in YSO.limited) {
  maxdc <- quantile(dat.obs$DeadConif[which(dat.obs$YSO == i)], prob = 0.95, type = 8, na.rm = T)
  dat.pred <- dat.pred %>% filter(!(YSO == i & DeadConif > maxdc))
}
dat.pred <- cbind(dat.pred, Y = predict.lm(mod, dat.pred))
p.heat <- p.HEAT(dat.pred, "% Herbaceous")

p.Herb.SF <- p.stitch.row(p.CC, p.DCon, p.heat)

p.Ground <- ggdraw() + 
  draw_plot(p.Woody.LP, x = .05, y = 0.6333334, width = .475, height = 0.3166667) +
  draw_plot(p.Woody.SF, x = .525, y = 0.6333334, width = .475, height = 0.3166667) +
  draw_plot(p.DDCov.LP, x = .05, y = 0.3166667, width = .475, height = 0.3166667) +
  draw_plot(p.DDCov.SF, x = .525, y = 0.3166667, width = .475, height = 0.3166667) +
  draw_plot(p.Herb.LP, x = .05, y = 0, width = .475, height = 0.3166667) +
  draw_plot(p.Herb.SF, x = .525, y = 0, width = .475, height = 0.3166667) +
  draw_plot_label(c("% Herbaceous", "% Dead and down", "% Woody",
                    "Lodgepole pine", "Spruce-fir"),
                  x = c(.03, .03, .03, .25, .7),
                  y = c(.12, .43, .77, .97, .97),
                  size = c(15, 15, 15, 17, 17),
                  angle = c(90, 90, 90, 0, 0),
                  hjust = c(0, 0, 0, 0, 0),
                  vjust = c(0, 0, 0, 0, 0))

save_plot("figure_Ground_VS_outbreak.tiff", p.Ground, ncol = 4, nrow = 2.5, dpi = 300)
