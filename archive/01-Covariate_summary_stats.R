require(dplyr)
require(stringr)
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

dat <- cov_tab_import %>%
  rename(YSI = YSO) %>%
  rename(RoadDens = Rd_dens1km) %>%
  mutate(Stratum = str_sub(Point, 1, 9)) %>%
  mutate(DeadConif = replace(DeadConif, which(DeadConif == 1.25), 1)) %>% # fix incorrect DeadCon value
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
cols <- c("LP.pt", "npt.LP", "LP.tr", "ntr.LP", "SF.pt", "npt.SF", "SF.tr", "ntr.SF")
out <- matrix("", nrow = length(rows), ncol  = length(cols),
              dimnames = list(rows, cols))

for(r in which(!rows %in% c("DeadConif_SD", "CanCov_SD", "WILD"))) {
  vec <- dat.LP[, rows[r]] %>% as.matrix %>% as.numeric
  out[rows[r], "LP.pt"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
          vec %>% min(na.rm = T) %>% round(digits = 2), "-",
          vec %>% max(na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "npt.LP"] <- sum(!is.na(vec))
  vec.tr <- vec %>%
    tapply(dat.LP$Transect, mean, na.rm = T)
  out[rows[r], "LP.tr"] <- vec.tr %>%
    mean(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec.tr %>% sd(na.rm = T) %>% round(digits = 2), ", ",
          vec.tr %>% min(na.rm = T) %>% round(digits = 2), "-",
          vec.tr %>% max(na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "ntr.LP"] <- vec %>% tapply(dat.LP$Transect, function(x)
    any(!is.na(x)) %>% as.integer) %>% sum
  vec <- dat.SF[, rows[r]] %>% as.matrix %>% as.numeric
  out[rows[r], "SF.pt"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
          vec %>% min(na.rm = T) %>% round(digits = 2), "-",
          vec %>% max(na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "npt.SF"] <- sum(!is.na(vec))
  vec.tr <- vec %>%
    tapply(dat.SF$Transect, mean, na.rm = T)
  out[rows[r], "SF.tr"] <- vec.tr %>%
    mean(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec.tr %>% sd(na.rm = T) %>% round(digits = 2), ", ",
          vec.tr %>% min(na.rm = T) %>% round(digits = 2), "-",
          vec.tr %>% max(na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "ntr.SF"] <- vec %>% tapply(dat.SF$Transect, function(x)
    any(!is.na(x)) %>% as.integer) %>% sum
}

# Grid-level heterogeneity metrics #
vec <- (dat.LP %>%
    group_by(Transect) %>%
    summarise(DeadConif = sd(DeadConif, na.rm = T)))$DeadConif
out["DeadConif_SD", "LP.tr"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["DeadConif_SD", "ntr.LP"] <- sum(!is.na(vec))
out["DeadConif_SD", "LP.pt"] <-
  out["DeadConif_SD", "npt.LP"] <- NA

vec <- (dat.SF %>%
          group_by(Transect) %>%
          summarise(DeadConif = sd(DeadConif, na.rm = T)))$DeadConif
out["DeadConif_SD", "SF.tr"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["DeadConif_SD", "ntr.SF"] <- sum(!is.na(vec))
out["DeadConif_SD", "SF.pt"] <-
  out["DeadConif_SD", "npt.SF"] <- NA

vec <- (dat.LP %>%
          group_by(Transect) %>%
          summarise(CanCov = sd(CanCov, na.rm = T)))$CanCov
out["CanCov_SD", "LP.tr"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["CanCov_SD", "ntr.LP"] <- sum(!is.na(vec))
out["CanCov_SD", "LP.pt"] <-
  out["CanCov_SD", "npt.LP"] <- NA

vec <- (dat.SF %>%
          group_by(Transect) %>%
          summarise(CanCov = sd(CanCov, na.rm = T)))$CanCov
out["CanCov_SD", "SF.tr"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
  str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
        vec %>% min(na.rm = T) %>% round(digits = 2), "-",
        vec %>% max(na.rm = T) %>% round(digits = 2), ")")
out["CanCov_SD", "ntr.SF"] <- sum(!is.na(vec))
out["CanCov_SD", "SF.pt"] <-
  out["CanCov_SD", "npt.SF"] <- NA

# Wilderness #
vec <- dat.LP[, "WILD"] %>% as.matrix %>% as.numeric
out["WILD", "LP.pt"] <- vec %>%
  (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
  str_c("%")
out["WILD", "npt.LP"] <- sum(!is.na(vec))
out["WILD", "LP.tr"] <- vec %>%
  tapply(dat.LP$Transect, mean) %>% round %>%
  (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
  str_c("%")
out["WILD", "ntr.LP"] <- vec %>% tapply(dat.LP$Transect, function(x)
  any(!is.na(x)) %>% as.integer) %>% sum
vec <- dat.SF[, "WILD"] %>% as.matrix %>% as.numeric
out["WILD", "SF.pt"] <- vec %>%
  (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
  str_c("%")
out["WILD", "npt.SF"] <- sum(!is.na(vec))
out["WILD", "SF.tr"] <- vec %>%
  tapply(dat.SF$Transect, mean) %>% round %>%
  (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
  str_c("%")
out["WILD", "ntr.SF"] <- vec %>% tapply(dat.SF$Transect, function(x)
  any(!is.na(x)) %>% as.integer) %>% sum

write.csv(out, "Covariate_summary.csv", row.names = T)
