require(dplyr)
require(stringr)
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

dat <- cov_tab_import %>%
  rename(YSI = YSO) %>%
  rename(RoadDens = Rd_dens1km) %>%
  mutate(Stratum = str_sub(Point, 1, 9)) %>%
  filter(!is.na(TWIP)) # The two grids dropped here were in burned areas.

dat.LP <- dat %>% filter(Stratum == "CO-CPW-LP") %>%
  mutate(Transect = Point %>% str_sub(8, 13)) %>%
  select(Point, Transect, DeadConif:RoadDens, TWIP)

# Exploration of relative canopy covers #
#rcov.tot <- (dat.LP %>% mutate(TOT = RCOV_ES + RCOV_SU + RCOV_Pine + RCOV_AS))$TOT
#hist(rcov.tot)
#sum(rcov.tot == 100, na.rm = T) / sum(!is.na(rcov.tot)) # 86% add up to 100
#sum(rcov.tot >= 90, na.rm = T) / sum(!is.na(rcov.tot)) # 94% are at least 90

# Exploration of relative ground covers #
#rcov.tot <- (dat.LP %>% mutate(TOT = gc_herb + gc_woody + gc_deadAndDown + gc_bare_litter))$TOT
#hist(rcov.tot)
#sum(rcov.tot == 100, na.rm = T) / sum(!is.na(rcov.tot)) # 74% add up to 100
#sum(rcov.tot >= 90, na.rm = T) / sum(!is.na(rcov.tot)) # 95% are at least 90

# Explore relationships of wilderness with other covariates #
#dat.LP$DeadConif_RCov[which(dat.LP$WILD == 1)] %>% mean(na.rm = T) # 26%
#dat.LP$DeadConif_RCov[which(dat.LP$WILD == 1)] %>% sd(na.rm = T)
#dat.LP$DeadConif_RCov[which(dat.LP$WILD == 0)] %>% mean(na.rm = T) # 18%
#dat.LP$DeadConif_RCov[which(dat.LP$WILD == 0)] %>% sd(na.rm = T)

dat.SF <- dat %>% filter(Stratum == "CO-CPW-SF") %>%
  mutate(Transect = Point %>% str_sub(8, 13)) %>%
  select(Point, Transect, DeadConif:RoadDens, TWIP)

# Exploration  of relative canopy covers #
#rcov.tot <- (dat.SF %>% mutate(TOT = RCOV_ES + RCOV_SU + RCOV_Pine + RCOV_AS))$TOT
#hist(rcov.tot)
#sum(rcov.tot == 100, na.rm = T) / sum(!is.na(rcov.tot)) # 90% add up to 100
#sum(rcov.tot >= 90, na.rm = T) / sum(!is.na(rcov.tot)) # 96% are at least 90

# Exploration of relative ground covers #
#rcov.tot <- (dat.SF %>% mutate(TOT = gc_herb + gc_woody + gc_deadAndDown + gc_bare_litter))$TOT
#hist(rcov.tot)
#sum(rcov.tot == 100, na.rm = T) / sum(!is.na(rcov.tot)) # 62% add up to 100
#sum(rcov.tot >= 90, na.rm = T) / sum(!is.na(rcov.tot)) # 91% are at least 90

# Explore relationships of wilderness with other covariates #
#dat.SF$DeadConif_RCov[which(dat.SF$WILD == 1)] %>% mean(na.rm = T) # 30%
#dat.SF$DeadConif_RCov[which(dat.SF$WILD == 1)] %>% sd(na.rm = T)
#dat.SF$DeadConif_RCov[which(dat.SF$WILD == 0)] %>% mean(na.rm = T) # 14%
#dat.SF$DeadConif_RCov[which(dat.SF$WILD == 0)] %>% sd(na.rm = T)

rows <- c("DeadConif", "DeadConif_SD", "YSI", "CanCov", "CanCov_SD", "RCOV_AS", "RCOV_ES", "RCOV_Pine", "shrub_cover", 
          "RCShrb_UC", "HerbCov", "WoodyCov", "DDCov", "RoadDens", "TWIP", "WILD")
cols <- c("LP", "ntr.LP", "npt.LP", "SF", "ntr.SF", "npt.SF")
out <- matrix("", nrow = length(rows), ncol  = length(cols),
              dimnames = list(rows, cols))

for(r in which(!rows %in% c("DeadConif_SD", "CanCov_SD", "WILD"))) {
  vec <- dat.LP[, rows[r]] %>% as.matrix %>% as.numeric
  out[rows[r], "LP"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
          vec %>% min(na.rm = T) %>% round(digits = 2), "-",
          vec %>% max(na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "npt.LP"] <- sum(!is.na(vec))
  out[rows[r], "ntr.LP"] <- vec %>% tapply(dat.LP$Transect, function(x)
    any(!is.na(x)) %>% as.integer) %>% sum
  vec <- dat.SF[, rows[r]] %>% as.matrix %>% as.numeric
  out[rows[r], "SF"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
          vec %>% min(na.rm = T) %>% round(digits = 2), "-",
          vec %>% max(na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "npt.SF"] <- sum(!is.na(vec))
  out[rows[r], "ntr.SF"] <- vec %>% tapply(dat.SF$Transect, function(x)
    any(!is.na(x)) %>% as.integer) %>% sum
}

# Grid-level Heterogeneity #
vec <- (dat.LP %>%
    group_by(Transect) %>%
    summarise(DeadConif = sd(DeadConif, na.rm = T)))$DeadConif
out["DeadConif_SD", "LP"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["DeadConif_SD", "ntr.LP"] <- sum(!is.na(vec))
out["DeadConif_SD", "npt.LP"] <- NA

vec <- (dat.SF %>%
          group_by(Transect) %>%
          summarise(DeadConif = sd(DeadConif, na.rm = T)))$DeadConif
out["DeadConif_SD", "SF"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["DeadConif_SD", "ntr.SF"] <- sum(!is.na(vec))
out["DeadConif_SD", "npt.SF"] <- NA

vec <- (dat.LP %>%
          group_by(Transect) %>%
          summarise(CanCov = sd(CanCov, na.rm = T)))$CanCov
out["CanCov_SD", "LP"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["CanCov_SD", "ntr.LP"] <- sum(!is.na(vec))
out["CanCov_SD", "npt.LP"] <- NA

vec <- (dat.SF %>%
          group_by(Transect) %>%
          summarise(CanCov = sd(CanCov, na.rm = T)))$CanCov
out["CanCov_SD", "SF"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["CanCov_SD", "ntr.SF"] <- sum(!is.na(vec))
out["CanCov_SD", "npt.SF"] <- NA

# Wilderness #
vec <- dat.LP[, "WILD"] %>% as.matrix %>% as.numeric
out["WILD", "LP"] <- vec %>%
  (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
  str_c("%")
out["WILD", "npt.LP"] <- sum(!is.na(vec))
out["WILD", "ntr.LP"] <- vec %>% tapply(dat.LP$Transect, function(x)
  any(!is.na(x)) %>% as.integer) %>% sum
vec <- dat.SF[, "WILD"] %>% as.matrix %>% as.numeric
out["WILD", "SF"] <- vec %>%
  (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
  str_c("%")
out["WILD", "npt.SF"] <- sum(!is.na(vec))
out["WILD", "ntr.SF"] <- vec %>% tapply(dat.SF$Transect, function(x)
  any(!is.na(x)) %>% as.integer) %>% sum

write.csv(out, "Covariate_summary.csv", row.names = T)

# Scatterplots - DeadConif VS YSO #
library(cowplot)

p.LP <- ggplot(data = dat.LP, aes(x = YSI, y = DeadConif_RCov)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 1.05, label = "Lodgpole pine stratum", size = 4)

p.SF <- ggplot(data = (dat.SF %>% filter(DeadConif_RCov <= 1)), aes(x = YSI, y = DeadConif_RCov)) +
  geom_point(alpha = 0.2) +
  xlab(NULL) + ylab(NULL) +
  annotate("text", x = 5, y = 1.05, label = "Spruce-fir stratum", size = 4)

p <- ggdraw() + 
  draw_plot(p.LP, x = 0.05, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot(p.SF, x = 0.525, y = 0.05, width = 0.475, height = 0.95) +
  draw_plot_label(c("Proportion trees dead", "Years since outbreak"),
                  x = c(0.03, 0.4),
                  y = c(0.2, 0.03),
                  size = c(15, 15),
                  angle = c(90, 0),
                  hjust = c(0, 0),
                  vjust = c(0, 0))

save_plot("figure_DeadConif_VS_YSO.tiff", p, ncol = 2, nrow = 1, dpi = 200)

### Summarize grid-level heterogeneity for select variables ###
## LP ##
# Dead conifer
tapply(dat.LP$DeadConif_RCov, dat.LP$Transect, sd, na.rm = T) %>% hist

# Canopy cover
tapply(dat.LP$CanCov, dat.LP$Transect, sd, na.rm = T) %>% hist

# Aspen
tapply(dat.LP$RCOV_AS, dat.LP$Transect, sd, na.rm = T) %>% hist

# Shrub cover
tapply(dat.LP$shrub_cover, dat.LP$Transect, sd, na.rm = T) %>% hist

## SF ##
# Dead conifer
tapply(dat.SF$DeadConif_RCov, dat.SF$Transect, sd, na.rm = T) %>% hist

# Canopy cover
tapply(dat.SF$CanCov, dat.SF$Transect, sd, na.rm = T) %>% hist

# Aspen
tapply(dat.SF$RCOV_AS, dat.SF$Transect, sd, na.rm = T) %>% hist

# Shrub cover
tapply(dat.SF$shrub_cover, dat.SF$Transect, sd, na.rm = T) %>% hist
