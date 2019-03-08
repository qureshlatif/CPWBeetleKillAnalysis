require(dplyr)
require(stringr)
setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

maxyso4PD <- 12
gridID <- Cov.LP[, "gridIndex"]
outbreak_grids <- tapply(Cov.LP[, "YSO"], Cov.LP[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
Cov.LP[which(!outbreak_grids[gridID]), "YSO"] <- -1
YSO.check <- Cov.LP[, "YSO"] %>% # This is for screening out late-outbreak PctDead values
  replace(., which(is.na(.)), tapply(., gridID, median, na.rm = T)[gridID][which(is.na(.))] %>% round)
dat.LP <- Cov.LP %>%
  as.data.frame %>%
  mutate(DeadConif = replace(DeadConif, which(YSO.check > maxyso4PD), NA)) %>%
  select(gridIndex, DayOfYear, Time, DeadConif:DDCov)

maxyso4PD <- 9
gridID <- Cov.SF[, "gridIndex"]
outbreak_grids <- tapply(Cov.SF[, "YSO"], Cov.SF[, "gridIndex"], function(x) any(!is.na(x))) # index grids intersecting ADS outbreaks
Cov.SF[which(!outbreak_grids[gridID]), "YSO"] <- -1
YSO.check <- Cov.SF[, "YSO"] %>% # This is for screening out late-outbreak PctDead values
  replace(., which(is.na(.)), tapply(., gridID, median, na.rm = T)[gridID][which(is.na(.))] %>% round)
dat.SF <- Cov.SF %>%
  as.data.frame %>%
  mutate(DeadConif = replace(DeadConif, which(YSO.check > maxyso4PD), NA)) %>%
  select(gridIndex, DayOfYear, Time, DeadConif:DDCov)

rm(gridID, outbreak_grids, maxyso4PD)

rows <- c("DeadConif", "YSO", "CanCov", "RCOV_AS", "RCOV_ES", "RCOV_Pine", "shrub_cover", 
          "RCShrb_UC", "HerbCov", "WoodyCov", "DDCov")
cols <- c("LP.pt", "npt.LP", "SF.pt", "npt.SF")
out <- matrix("", nrow = length(rows), ncol  = length(cols),
              dimnames = list(rows, cols))

for(r in 1:length(rows)) {
  vec <- dat.LP[, rows[r]] %>% as.matrix %>% as.numeric
  out[rows[r], "LP.pt"] <- vec %>% median(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec %>% quantile(prob = 0.05, type = 8, na.rm = T) %>% round(digits = 2), ", ",
          vec %>% quantile(prob = 0.25, type = 8, na.rm = T) %>% round(digits = 2), ", ",
          vec %>% quantile(prob = 0.75, type = 8, na.rm = T) %>% round(digits = 2), ", ",
          vec %>% quantile(prob = 0.95, type = 8, na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "npt.LP"] <- sum(!is.na(vec))
  vec <- dat.SF[, rows[r]] %>% as.matrix %>% as.numeric
  out[rows[r], "SF.pt"] <- vec %>% median(na.rm = T) %>% round(digits = 2) %>%
    str_c(" (", vec %>% quantile(prob = 0.05, type = 8, na.rm = T) %>% round(digits = 2), ", ",
          vec %>% quantile(prob = 0.25, type = 8, na.rm = T) %>% round(digits = 2), ", ",
          vec %>% quantile(prob = 0.75, type = 8, na.rm = T) %>% round(digits = 2), ", ",
          vec %>% quantile(prob = 0.95, type = 8, na.rm = T) %>% round(digits = 2), ")")
  out[rows[r], "npt.SF"] <- sum(!is.na(vec))
}

write.csv(out, "Covariate_summary.csv", row.names = T)

