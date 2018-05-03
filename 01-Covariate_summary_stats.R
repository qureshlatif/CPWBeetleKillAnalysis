require(dplyr)
require(stringr)
setwd("F:/research stuff/BCR/projects/CPW")

dat <- read.csv("Covariates.csv", header = T, stringsAsFactors = F) %>%
  tbl_df() %>%
  mutate(Stratum = Point %>% str_sub(1, 9)) %>%
  mutate(YSI = YSI %>% replace(which(YSI == -1), NA)) %>%
  filter(!is.na(TWIP)) # The two grids dropped here were in burned areas.

dat.LP <- dat %>% filter(Stratum == "CO-CPW-LP") %>%
  mutate(Transect = Point %>% str_sub(8, 13)) %>%
  select(Point, Transect, TWIP, CanCov, DeadConif_RCov,
         YSI, RCOV_ES:HerbCov, primaryHabitat)

dat.SF <- dat %>% filter(Stratum == "CO-CPW-SF") %>%
  mutate(Transect = Point %>% str_sub(8, 13)) %>%
  select(Point, Transect, TWIP, CanCov, DeadConif_RCov,
         YSI, RCOV_ES:HerbCov, primaryHabitat)

rows <- c("DeadConif_RCov", "YSI", "CanCov", "RCOV_AS", "RCOV_ES", "RCOV_Pine", "shrub_cover", 
          "RCShrb_UC", "RCShrb_UD", "HerbCov", "TWIP")
cols <- c("LP", "ntr.LP", "npt.LP", "SF", "ntr.SF", "npt.SF")
out <- matrix("", nrow = length(rows), ncol  = length(cols),
              dimnames = list(rows, cols))

for(r in 1:length(rows)) {
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

write.csv(out, "Covariate_summary.csv", row.names = T)
