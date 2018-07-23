require(dplyr)
require(stringr)
require(ggplot2)
require(cowplot)
require(QSLpersonal)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

dat.LP <- Cov.LP %>% tbl_df %>%
  mutate(Point = row.names(Cov.LP)) %>%
  select(Point, gridIndex:Rd_dens1km) %>%
  mutate(DeadConif = DeadConif * 100) %>%
  left_join(Cov.LP %>% tbl_df %>%
              group_by(gridIndex) %>%
              summarise(Outbreak = any(!is.na(YSO)) %>% as.integer,
                        YSO.med = median(YSO, na.rm = T)), by = "gridIndex") %>%
  mutate(YSO = replace(YSO, which(is.na(YSO) & Outbreak == 0), 0))
ind <- which(is.na(dat.LP$YSO) & dat.LP$Outbreak == 1)
dat.LP$YSO[ind] <- dat.LP$YSO.med[ind]
dat.LP <- dat.LP %>% select(-YSO.med) %>%
  mutate(YSO = replace(YSO, which(Outbreak == 0), -1)) %>%
  # consolidate to grid level
  group_by(gridIndex) %>%
  summarise(DeadConif = DeadConif %>% mean(na.rm = T),
            YSO = YSO %>% mean(),
            CanCov = CanCov %>% mean(na.rm = T),
            RCOV_AS = RCOV_AS %>% mean(na.rm = T),
            RCOV_ES = RCOV_ES %>% mean(na.rm = T),
            RCOV_Pine = RCOV_Pine %>% mean(na.rm = T),
            shrub_cover = shrub_cover %>% mean(na.rm = T),
            RCShrb_UC = RCShrb_UC %>% mean(na.rm = T),
            HerbCov = HerbCov %>% mean(na.rm = T),
            WoodyCov = WoodyCov %>% mean(na.rm = T),
            DDCov = DDCov %>% mean(na.rm = T),
            WILD = WILD %>% mean(na.rm = T) %>% round,
            Rd_dens1km = Rd_dens1km %>% mean(na.rm = T),
            Outbreak = Outbreak %>% mean(na.rm = T))

dat.SF <- Cov.SF %>% tbl_df %>%
  mutate(Point = row.names(Cov.SF)) %>%
  select(Point, gridIndex:Rd_dens1km) %>%
  mutate(DeadConif = DeadConif * 100) %>%
  left_join(Cov.SF %>% tbl_df %>%
              group_by(gridIndex) %>%
              summarise(Outbreak = any(!is.na(YSO)) %>% as.integer,
                        YSO.med = median(YSO, na.rm = T)), by = "gridIndex") %>%
  mutate(YSO = replace(YSO, which(is.na(YSO) & Outbreak == 0), 0))
ind <- which(is.na(dat.SF$YSO) & dat.SF$Outbreak == 1)
dat.SF$YSO[ind] <- dat.SF$YSO.med[ind]
dat.SF <- dat.SF %>% select(-YSO.med) %>%
  mutate(YSO = replace(YSO, which(Outbreak == 0), -1)) %>%
  # consolidate to grid level
  group_by(gridIndex) %>%
  summarise(DeadConif = DeadConif %>% mean(na.rm = T),
            YSO = YSO %>% mean(),
            CanCov = CanCov %>% mean(na.rm = T),
            RCOV_AS = RCOV_AS %>% mean(na.rm = T),
            RCOV_ES = RCOV_ES %>% mean(na.rm = T),
            RCOV_Pine = RCOV_Pine %>% mean(na.rm = T),
            shrub_cover = shrub_cover %>% mean(na.rm = T),
            RCShrb_UC = RCShrb_UC %>% mean(na.rm = T),
            HerbCov = HerbCov %>% mean(na.rm = T),
            WoodyCov = WoodyCov %>% mean(na.rm = T),
            DDCov = DDCov %>% mean(na.rm = T),
            WILD = WILD %>% mean(na.rm = T) %>% round,
            Rd_dens1km = Rd_dens1km %>% mean(na.rm = T),
            Outbreak = Outbreak %>% mean(na.rm = T))

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

save_plot("figure_DeadConif_VS_YSO_grid.tiff", p, ncol = 2, nrow = 1, dpi = 200)

#### Variance inflation factors for habitat attributes ####
vif.LP <- VIF(dat.LP, vars = names(dat.LP)[-c(1:3, 6)])
vif.SF <- VIF(dat.SF, vars = names(dat.SF)[-c(1:3)])

#### Regression of vegetation (mechanistic) factors VS outbreak metrics ####
vars <- c("CanCov", "RCOV_AS", "RCOV_ES", "RCOV_Pine", "shrub_cover", "RCShrb_UC", "HerbCov", "WoodyCov", "DDCov")
cols <- c("Int.LP", "YSO.LP", "DCon.LP", "Int.SF", "YSO.SF", "DCon.SF")
out <- matrix(NA, nrow = length(vars), ncol = length(cols),
              dimnames = list(vars, cols))

for(i in 1:length(vars)) {
  mod <- lm(as.formula(str_c(vars[i], "~ YSO + DeadConif")), data = dat.LP)
  sum.tab <- summary(mod)$coefficients
  out[vars[i], c("Int.LP", "YSO.LP", "DCon.LP")] <-
    str_c(round(sum.tab[, "Estimate"], digits = 3), " (",
          round(sum.tab[, "Std. Error"], digits = 3), ")")
  if(any(sum.tab[-1, "Pr(>|t|)"] < 0.05)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] < 0.05)
    out[vars[i], c("YSO.LP", "DCon.LP")][sig.ind] <-
      str_c(out[vars[i], c("YSO.LP", "DCon.LP")][sig.ind], "*")
  }
  if(any(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)
    out[vars[i], c("YSO.LP", "DCon.LP")][sig.ind] <-
      str_c(out[vars[i], c("YSO.LP", "DCon.LP")][sig.ind], ".")
  }
  mod <- lm(as.formula(str_c(vars[i], "~ YSO + DeadConif")), data = dat.SF)
  sum.tab <- summary(mod)$coefficients
  out[vars[i], c("Int.SF", "YSO.SF", "DCon.SF")] <-
    str_c(round(sum.tab[, "Estimate"], digits = 3), " (",
          round(sum.tab[, "Std. Error"], digits = 3), ")")
  if(any(sum.tab[-1, "Pr(>|t|)"] < 0.05)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] < 0.05)
    out[vars[i], c("YSO.SF", "DCon.SF")][sig.ind] <-
      str_c(out[vars[i], c("YSO.SF", "DCon.SF")][sig.ind], "*")
  }
  if(any(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)) {
    sig.ind <- which(sum.tab[-1, "Pr(>|t|)"] >= 0.05 & sum.tab[-1, "Pr(>|t|)"] < 0.1)
    out[vars[i], c("YSO.SF", "DCon.SF")][sig.ind] <-
      str_c(out[vars[i], c("YSO.SF", "DCon.SF")][sig.ind], ".")
  }
}

# Rearrange #
out <- rbind(out[, 1:3], out[, 4:6]) %>%
  cbind(Stratum = c(rep("LP", 9), rep("SF", 9)), .) %>%
  cbind(HabFeature = dimnames(out)[[1]], .) %>%
  tbl_df %>%
  slice((1:18) %>% matrix(ncol = 2) %>% t %>% as.integer)

write.csv(out, "Outbreak_habitat_relations_grid.csv", row.names = T)

plot(dat.LP$YSO, dat.LP$RCOV_AS)
plot(dat.SF$YSO, dat.SF$RCOV_AS)
plot(dat.SF$DeadConif, dat.SF$RCOV_ES)
