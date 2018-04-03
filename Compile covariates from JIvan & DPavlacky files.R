require(dplyr)
require(stringr)
setwd("F:/research stuff/BCR/projects/CPW")

### Compile covariates from Jake Ivan's file ###
dat <- read.csv("JIvan_files/BCR_Dead_Conifer_Qry_Headers.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>%
  mutate(Point = str_c(RMBOName, "-", str_pad(PointNumber, width = 2, side = "left", pad = "0"))) %>%
  select(Point, PctAbundance:SpeciesCode)

dat.Ivan <- dat %>%
  group_by(Point) %>%
  summarise(DeadConif = sum(PctAbundance*(PctRed + PctSilver)*0.01^2)) %>%
  ungroup

add <- dat %>% filter(SpeciesCode == "LP") %>%
  group_by(Point) %>%
  summarise(CovLP = sum(PctAbundance)) %>%
  ungroup
dat.Ivan <- dat.Ivan %>% left_join(add, by = "Point")

add <- dat %>% filter(SpeciesCode %in% "ES") %>%
  group_by(Point) %>%
  summarise(CovES = sum(PctAbundance)) %>%
  ungroup
dat.Ivan <- dat.Ivan %>% left_join(add, by = "Point")

add <- dat %>% filter(SpeciesCode %in% "SU") %>%
  group_by(Point) %>%
  summarise(CovSU = sum(PctAbundance)) %>%
  ungroup
dat.Ivan <- dat.Ivan %>% left_join(add, by = "Point")

dat.Ivan <- dat.Ivan %>%
  mutate(CovLP = replace(CovLP, is.na(CovLP), 0)) %>%
  mutate(CovES = replace(CovES, is.na(CovES), 0)) %>%
  mutate(CovSU = replace(CovSU, is.na(CovSU), 0))

# Error check -- look for duplicate records by point #
#checkDups <- dat %>% group_by(Point) %>%
#  summarise(noLP = sum(SpeciesCode == "LP"),
#            noES = sum(SpeciesCode == "ES"),
#            noSU = sum(SpeciesCode == "SU")) %>%
#  filter(noLP > 1 | noES > 1 | noSU > 1)
#dat.check <- dat %>% filter(Point %in% checkDups$Point)
#write.csv(dat.check, "JIvan_files/Points_w_dupSppCovers.csv", row.names = F)

# Explore cover by stratum #
dat.Ivan <- dat.Ivan %>% mutate(Stratum = str_sub(Point, 8, 9))
hist(dat.Ivan$CovLP[which(dat.Ivan$Stratum == "LP")])
hist(dat.Ivan$CovLP[which(dat.Ivan$Stratum == "SF")])
dat.Ivan %>% filter(CovLP > 100 & Stratum == "SF") %>% View

hist(dat.Ivan$CovES[which(dat.Ivan$Stratum == "SF")])
hist(dat.Ivan$CovSU[which(dat.Ivan$Stratum == "SF")])
hist(dat.Ivan$CovES[which(dat.Ivan$Stratum == "LP")])
hist(dat.Ivan$CovSU[which(dat.Ivan$Stratum == "LP")])

### Compile covariates from D. Pavlacky ###
dat <- read.csv("CPW_veg_query.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>% 
  mutate(Point = str_c(TransectNum, "-", str_pad(Point, width = 2, side = "left", pad = "0"))) %>%
  select(TransectNum, Point, o_canopy_percent, OVSPP:primaryHabitat)

# Check for duplicates #
dups <- dat[dat %>% select(Point, OVSPP, SHSPP) %>% duplicated %>% which, ] %>%
  mutate(ID = str_c(Point, "_", OVSPP, "_", SHSPP))
dat %>% mutate(ID = str_c(Point, "_", OVSPP, "_", SHSPP)) %>%
  filter(ID %in% dups$ID) %>%
  arrange(ID) %>%
  select(-ID) %>%
  write.csv("CPW_Veg_query_dups.csv", row.names = F)
