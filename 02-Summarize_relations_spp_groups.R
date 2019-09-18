library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")
spp.list[which(spp.list == "GRAJ")] <- "CAJA"

# Get spp table with predictions and life histories #
spp_pred <- read.csv("Spp_apriori.csv", header = T, stringsAsFactors = F) %>% tbl_df

# Get outbreak models
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

# Categorize estimated relationships #
beta.sum <- function(x) 
  if(quantile(x, type = 8, prob = 0.05) > 0 | quantile(x, type = 8, prob = 0.95) < 0) {
    str_c(
    median(x) %>% round(digits = 2),
    "(",
    quantile(x, type = 8, prob = 0.05) %>% round(digits = 2),
    ",",
    quantile(x, type = 8, prob = 0.95) %>% round(digits = 2),
    ")*")
  } else {
    str_c(
      median(x) %>% round(digits = 2),
      "(",
      quantile(x, type = 8, prob = 0.05) %>% round(digits = 2),
      ",",
      quantile(x, type = 8, prob = 0.95) %>% round(digits = 2),
      ")")
  }

relations_cat_LP <- data.frame(Spp = spp.list,
                               Outbreak_pos_LP = as.numeric(parest_LP[, "bb.pdead.lo"] > 0 |
                                                              parest_LP[, "bb.YSO.lo"] > 0 |
                                                              parest_LP[, "bb.YSO2.lo"] > 0), 
                               Outbreak_neg_LP = as.numeric(parest_LP[, "bb.pdead.hi"] < 0 |
                                                              parest_LP[, "bb.YSO.hi"] < 0 |
                                                              parest_LP[, "bb.YSO2.hi"] < 0 |
                                                              parest_LP[, "bb.pdXYSO.hi"] < 0),
                               YSO_lag_pos_LP = as.numeric(parest_LP[, "bb.YSO2.lo"] > 0),
                               YSO_peak_pos_LP = as.numeric(parest_LP[, "bb.YSO2.hi"] < 0 &
                                                           parest_LP[, "bb.YSO"] > 0),
                               DConXYSO_neg_LP = as.numeric(parest_LP[, "bb.pdXYSO.hi"] < 0),
                               any_relation_LP = as.numeric(parest_LP[, "bb.pdead.lo"] > 0 |
                                                           parest_LP[, "bb.pdead.hi"] < 0 |
                                                           parest_LP[, "bb.YSO.lo"] > 0 |
                                                           parest_LP[, "bb.YSO.hi"] < 0 |
                                                           parest_LP[, "bb.YSO2.lo"] > 0 |
                                                           parest_LP[, "bb.YSO2.hi"] < 0 |
                                                           parest_LP[, "bb.pdXYSO.lo"] > 0 |
                                                           parest_LP[, "bb.pdXYSO.hi"] < 0),
                               stringsAsFactors = F) %>% tbl_df

relations_cat_SF <- data.frame(Spp = spp.list,
                               Outbreak_pos_SF = as.numeric(parest_SF[, "bb.pdead.lo"] > 0 |
                                                              parest_SF[, "bb.YSO.lo"] > 0 |
                                                              parest_SF[, "bb.YSO2.lo"] > 0), 
                               Outbreak_neg_SF = as.numeric(parest_SF[, "bb.pdead.hi"] < 0 |
                                                              parest_SF[, "bb.YSO.hi"] < 0 |
                                                              parest_SF[, "bb.YSO2.hi"] < 0 |
                                                              parest_SF[, "bb.pdXYSO.hi"] < 0),
                               YSO_lag_pos_SF = as.numeric(parest_SF[, "bb.YSO2.lo"] > 0),
                               YSO_peak_pos_SF = as.numeric(parest_SF[, "bb.YSO2.hi"] < 0 &
                                                              parest_SF[, "bb.YSO"] > 0),
                               DConXYSO_neg_SF = as.numeric(parest_SF[, "bb.pdXYSO.hi"] < 0),
                               any_relation_SF = as.numeric(parest_SF[, "bb.pdead.lo"] > 0 |
                                                              parest_SF[, "bb.pdead.hi"] < 0 |
                                                              parest_SF[, "bb.YSO.lo"] > 0 |
                                                              parest_SF[, "bb.YSO.hi"] < 0 |
                                                              parest_SF[, "bb.YSO2.lo"] > 0 |
                                                              parest_SF[, "bb.YSO2.hi"] < 0 |
                                                              parest_SF[, "bb.pdXYSO.lo"] > 0 |
                                                              parest_SF[, "bb.pdXYSO.hi"] < 0),
                               stringsAsFactors = F) %>% tbl_df

spp_pred <- spp_pred %>%
  left_join(relations_cat_LP, by = "Spp") %>%
  left_join(relations_cat_SF, by = "Spp")
write.csv(spp_pred, "Species_pred&obs_outbreak_summary.csv", row.names = F)

rm(relations_cat_LP, relations_cat_SF)

## Summary tables ##
min(spp_pred$LP_det[which(spp_pred$any_relation_LP == 1)])
min(spp_pred$SF_det[which(spp_pred$any_relation_SF == 1)])

# Outbreak relationships #
rows <- c("Pos_outbreak", "Neg_outbreak", "Unk_outbreak", "LagPos_YSO", "PeakPos_YSO",
          "Cav_Nest", "CP_Nest", "CS_Nest", "OC_Nest", "OU_Nest",
          "AI_Forage", "U_Forage", "FB_Forage", "S_Forage", "CS_Forage", "Predator")
cols <- c("n", "n10", "Outbreak_pos", "Outbreak_neg",
          "YSO_lag_pos", "YSO_peak_pos", "DConXYSO_neg",
          "bb.pdead.mn", "bb.YSO.mn", "bb.YSO2.mn", "bb.pdXYSO.mn")
sum.LP <- sum.SF <- matrix(NA, nrow = length(rows), ncol = length(cols),
                           dimnames = list(rows, cols))

sp.lsts <- list(Pos_outbreak = spp_pred %>% filter(Dcon == "+" | YSO %in% c("+", "peaked +", "lagged +")) %>% pull(Spp),
                Neg_outbreak = spp_pred %>% filter(Dcon == "-" | YSO == "-") %>% pull(Spp),
                Unk_outbreak = spp_pred %>% filter(Dcon == "U" | YSO == "U") %>% pull(Spp),
                LagPos_YSO = spp_pred %>% filter(YSO == "lagged +") %>% pull(Spp),
                PeakPos_YSO = spp_pred %>% filter(YSO == "peaked +") %>% pull(Spp),
                Cav_Nest = spp_pred %>% filter(LH_Nest %in% c("CP", "CS")) %>% pull(Spp),
                CP_Nest = spp_pred %>% filter(LH_Nest == "CP") %>% pull(Spp),
                CS_Nest = spp_pred %>% filter(LH_Nest == "CS") %>% pull(Spp),
                OC_Nest = spp_pred %>% filter(LH_Nest %in% c("OC", "OC,U")) %>% pull(Spp),
                OU_Nest = spp_pred %>% filter(LH_Nest %in% c("OU", "OC,U")) %>% pull(Spp),
                AI_Forage = spp_pred %>% filter(LH_Forage == "AI") %>% pull(Spp),
                U_Forage = spp_pred %>% filter(LH_Forage == "U") %>% pull(Spp),
                FB_Forage = spp_pred %>% filter(LH_Forage %in% c("FB", "FB, CSA")) %>% pull(Spp),
                S_Forage = spp_pred %>% filter(LH_Forage == "S") %>% pull(Spp),
                CS_Forage = spp_pred %>% filter(LH_Forage %in% c("CS", "FB, CSA")) %>% pull(Spp),
                Predator = c("CAJA", "STJA", "CLNU", "BBMA", "AMCR", "CORA", "HOWR", "BHCO"))

for(i in 1:length(rows)) {
  sp <- sp.lsts[[i]]
  sp.ind <- which(spp_pred$Spp %in% sp)
  sum.LP[i, "n"] <- length(sp)
  sum.LP[i, "n10"] <- sum(spp_pred$LP_det[sp.ind] >= 10)
  sum.LP[i, "Outbreak_pos"] <- sum(spp_pred$Outbreak_pos_LP[sp.ind] == 1)
  sum.LP[i, "Outbreak_neg"] <- sum(spp_pred$Outbreak_neg_LP[sp.ind] == 1)
  sum.LP[i, "YSO_lag_pos"] <- sum(spp_pred$YSO_lag_pos_LP[sp.ind] == 1)
  sum.LP[i, "YSO_peak_pos"] <- sum(spp_pred$YSO_peak_pos_LP[sp.ind] == 1)
  sum.LP[i, "DConXYSO_neg"] <- sum(spp_pred$DConXYSO_neg_LP[sp.ind] == 1)
  sum.SF[i, "n"] <- length(sp)
  sum.SF[i, "n10"] <- sum(spp_pred$SF_det[sp.ind] >= 10)
  sum.SF[i, "Outbreak_pos"] <- sum(spp_pred$Outbreak_pos_SF[sp.ind] == 1)
  sum.SF[i, "Outbreak_neg"] <- sum(spp_pred$Outbreak_neg_SF[sp.ind] == 1)
  sum.SF[i, "YSO_lag_pos"] <- sum(spp_pred$YSO_lag_pos_SF[sp.ind] == 1)
  sum.SF[i, "YSO_peak_pos"] <- sum(spp_pred$YSO_peak_pos_SF[sp.ind] == 1)
  sum.SF[i, "DConXYSO_neg"] <- sum(spp_pred$DConXYSO_neg_SF[sp.ind] == 1)
  sp.ind <- which(spp.list %in% sp)

  bb.mn <- mod_LP$sims.list$bb.pdead[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "bb.pdead.mn"] <- beta.sum(bb.mn)

  bb.mn <- mod_LP$sims.list$bb.YSO[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "bb.YSO.mn"] <- beta.sum(bb.mn)

  bb.mn <- mod_LP$sims.list$bb.YSO2[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "bb.YSO2.mn"] <- beta.sum(bb.mn)

  bb.mn <- mod_LP$sims.list$bb.pdXYSO[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "bb.pdXYSO.mn"] <- beta.sum(bb.mn)

  bb.mn <- mod_SF$sims.list$bb.pdead[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "bb.pdead.mn"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.YSO[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "bb.YSO.mn"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.YSO2[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "bb.YSO2.mn"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.pdXYSO[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "bb.pdXYSO.mn"] <- beta.sum(bb.mn)
}

write.csv(sum.LP, "SpGroup_outbreak_summary_LP.csv", row.names = T)
write.csv(sum.SF, "SpGroup_outbreak_summary_SF.csv", row.names = T)

# Habitat relationships #
  # Get models
mod_LP <- loadObject("mod_LPcommunity_habitat_reduced")
mod_SF <- loadObject("mod_SFcommunity_habitat_reduced")

cols <- c("n", "n10", "CanCov", "Aspen", "Spruce", "Pine", "ShrubCov", "ConShrb",
          "Herb", "Woody", "DeadDown")
sum.LP <- sum.SF <- matrix(NA, nrow = length(rows), ncol = length(cols),
                           dimnames = list(rows, cols))

for(i in 1:length(rows)) {
  sp <- sp.lsts[[i]]
  sp.ind <- which(spp_pred$Spp %in% sp)
  sum.LP[i, "n"] <- length(sp)
  sum.LP[i, "n10"] <- sum(spp_pred$LP_det[sp.ind] >= 10)
  sum.SF[i, "n"] <- length(sp)
  sum.SF[i, "n10"] <- sum(spp_pred$SF_det[sp.ind] >= 10)
  sp.ind <- which(spp.list %in% sp)
  
  bb.mn <- mod_LP$sims.list$bb.CanCov[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "CanCov"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_LP$sims.list$bb.RCovAS[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "Aspen"] <- beta.sum(bb.mn)

  bb.mn <- mod_LP$sims.list$bb.RCovPine[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "Pine"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_LP$sims.list$bb.ShCov[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "ShrubCov"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_LP$sims.list$bb.RSC_Con[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "ConShrb"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_LP$sims.list$bb.GHerb[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "Herb"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_LP$sims.list$bb.Gwoody[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "Woody"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_LP$sims.list$bb.GDD[, sp.ind] %>%
    apply(1, mean)
  sum.LP[i, "DeadDown"] <- beta.sum(bb.mn)

  bb.mn <- mod_SF$sims.list$bb.CanCov[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "CanCov"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.RCovAS[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "Aspen"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.RCovES[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "Spruce"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.RCovPine[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "Pine"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.ShCov[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "ShrubCov"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.RSC_Con[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "ConShrb"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.GHerb[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "Herb"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.Gwoody[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "Woody"] <- beta.sum(bb.mn)
  
  bb.mn <- mod_SF$sims.list$bb.GDD[, sp.ind] %>%
    apply(1, mean)
  sum.SF[i, "DeadDown"] <- beta.sum(bb.mn)
}

write.csv(sum.LP, "SpGroup_habitat_summary_LP.csv", row.names = T)
write.csv(sum.SF, "SpGroup_habitat_summary_SF.csv", row.names = T)
