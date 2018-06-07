require(dplyr)
require(stringr)
require(ggplot2)
require(cowplot)

setwd("C:/Users/Quresh.Latif/files/projects/CPW")

dat <- read.csv("Covariates.csv", header = T, stringsAsFactors = F) %>%
  tbl_df() %>%
  mutate(Stratum = Point %>% str_sub(1, 9)) %>%
  mutate(YSI = YSI %>% replace(which(YSI == -1), NA)) %>%
  filter(!is.na(TWIP)) # The two grids dropped here were in burned areas.

dat.LP <- dat %>% filter(Stratum == "CO-CPW-LP") %>%
  mutate(Transect = Point %>% str_sub(8, 13)) %>%
  select(Point, Transect, TWIP, CanCov, DeadConif_RCov,
         YSI, WILD, RCOV_ES:gc_bare_litter, primaryHabitat) %>%
  mutate(YSI = replace(YSI, which(is.na(YSI)), -2)) %>%
  filter(!is.na(DeadConif_RCov))

dat.SF <- dat %>% filter(Stratum == "CO-CPW-SF") %>%
  mutate(Transect = Point %>% str_sub(8, 13)) %>%
  select(Point, Transect, TWIP, CanCov, DeadConif_RCov,
         YSI, WILD, RCOV_ES:gc_bare_litter, primaryHabitat) %>%
  mutate(YSI = replace(YSI, which(is.na(YSI)), -2)) %>%
  filter(!is.na(DeadConif_RCov))

p <- ggplot(dat.LP, aes(x = YSI, y = DeadConif_RCov, color = CanCov)) +
  geom_point()
