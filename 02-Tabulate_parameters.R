library(jagsUI)
library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")
load("Data_compiled.RData")

#__________ Script inputs _____________#
stratum <- "LP" # Select LP or SF
mod <- loadObject("mod_LPcommunity_outbreak_QuadAtPnt (marg converge)")
params <- c("bd.pdead",
            #"bd.pdead2",
            "bd.outbrk",
            "bd.YSO",
            #"bd.YSO2",
            #"bd.pdXYSO",
            "bd.TWIP",
            "bb.pdead",
            "bb.pdead2",
            "bb.outbrk",
            "bb.YSO",
            "bb.YSO2",
            #"bb.pdXYSO",
            "bb.TWIP",
            "ba.Time",
            "ba.Time2",
            "ba.DOY",
            "ba.DOY2",
            "ba.ccov",
            "ba.shcov")
out.vals <- c("est", "f")
#______________________________________#

cols <- (expand.grid(out.vals, params, stringsAsFactors = F) %>%
  select(Var2, Var1) %>%
  mutate(Var3 = str_c(Var2, Var1, sep = ".")))$Var3
out <- matrix(NA, nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols))

for(i in 1:length(params)) {
  parm <- mod$sims.list[[params[i]]]
  out[, (i*2 - 1)] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
    ")")
  out[, (i*2)] <- apply(parm, 2, function(x) max(c(sum(x > 0), sum(x < 0))) / length(x)) %>%
    round(digits = 2)
}

write.csv(out, "Parameter_est.csv", row.names = T)
