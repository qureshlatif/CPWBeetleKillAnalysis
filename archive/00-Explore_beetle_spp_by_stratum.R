require(dplyr)
require(stringr)
setwd("C:/Users/Quresh.Latif/files/projects/CPW")

### Compile beetle covariates from Jake Ivan's files ###
dat <- read.csv("JIvan_files/BCR_Dead_Conifer_Qry_Headers.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>%
  mutate(Point = str_c(RMBOName, "-", str_pad(PointNumber, width = 2, side = "left", pad = "0"))) %>%
  select(Point, PctAbundance:SpeciesCode)

newdat <- dat %>%
  group_by(Point) %>%
  summarise(DeadConif = sum(PctAbundance*(PctRed + PctSilver)*0.01^2)) %>%
  ungroup

add <- dat %>% filter(SpeciesCode == "LP") %>%
  group_by(Point) %>%
  summarise(DeadLP = sum(PctAbundance*(PctRed + PctSilver)*0.01^2)) %>%
  ungroup
newdat <- newdat %>% left_join(add, by = "Point")

add <- dat %>% filter(SpeciesCode %in% c("ES", "SU")) %>%
  group_by(Point) %>%
  summarise(DeadSF = sum(PctAbundance*(PctRed + PctSilver)*0.01^2)) %>%
  ungroup
newdat <- newdat %>% left_join(add, by = "Point")

# Explore relation between stratum and beetle type #
newdat <- newdat %>% mutate(Stratum = str_sub(Point, 8, 9))
max(newdat$DeadSF[which(newdat$Stratum == "LP")], na.rm = T)
sort(unique(newdat$DeadLP[which(newdat$Stratum == "SF")]))
