require(dplyr)
require(stringr)
setwd("F:/research stuff/BCR/projects/CPW")

### Get covariates from Jake Ivan's file ###
dat.Ivan <- read.csv("JIvan_files/PCPointsWithConifData.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>%
  mutate(Point = str_c(Transect, "-", str_pad(PointNumber, width = 2, side = "left", pad = "0"))) %>%
  select(Transect, Point, PointCountDate, CanopyPct:SilverConifCover)
dat.missing <- read.csv("JIvan_files/MissingHabitatData.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>%
  mutate(Point = str_c(Transect, "-", str_pad(Point, width = 2, side = "left", pad = "0")))
dat.Ivan <- dat.Ivan %>% filter(!Point %in% dat.missing$Point)
#dat.Ivan %>% select(CanopyPct:SilverConifCover) %>% as.data.frame %>% # Check number of NAs
#  as.matrix %>% apply(2, function(x) sum(is.na(x)))

check <- dat.Ivan %>%
  mutate(DeadConifPct2 = RedConifPct + SilverConifPct) %>%
  mutate(Diff = (DeadConifPct2 - DeadConifPct) %>% round(digits = 5))
blanksToFill <- (dat.Ivan %>% filter(check$Diff %>% is.na))$Point
#dat.Ivan %>% filter(check$Diff %>% is.na)) %>%
#  write.csv("JIvan_files/blanks.csv", row.names = F)
#blanksToFill <- read.csv("JIvan_files/blanks.csv", header=T, stringsAsFactors = F)$Point

# Corrections #
dat.Ivan <- dat.Ivan %>%
  mutate(CanopyPct = replace(CanopyPct, which(Point %in% blanksToFill), 0))

blanks_more <- (dat.Ivan %>% filter(Point %in%
                      (dat.Ivan %>%
                         filter(is.na(dat.Ivan$CanopyPct)) %>%
                         filter(!Point %in% blanksToFill))$Point))$Point
dat.Ivan <- dat.Ivan %>%
  mutate(DeadConifCover = replace(DeadConifCover, which(Point %in% blanks_more), 0))

dat.DeadConif <- dat.Ivan %>% select(Point, DeadConifPct, DeadConifCover) %>%
  rename(DeadConif_RCov = DeadConifPct, DeadConif_Cov = DeadConifCover)

### Compile covariates from D. Pavlacky ###
# Plant codes #
# AL = Alder
# AS = Aspen
# BC = burnt conifer
# BI = Birch
# BR = Bristlecone pine
# BS = Blue Spruce
# DF = Douglas fir
# ES = Engelmann spruce
# JU = Juniper
# LM = Limber pine
# MA = ?? (4 records; OVABUND = -1 for CO-CPW-LP-053-14)
# PC = Plains Cottonwood
# PP = Ponderosa pine
# PY = Pinyon Pine
# SU = Subalpine fir
# UC = uknown conifer
# UD = unknown deciduous
# WF = White fir
# WI = willow
#_____________#
dat.overstory <- read.csv("CPW_veg_query_2_overstory.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>% 
  mutate(Point = str_c(TransectNum, "-", str_pad(Point, width = 2, side = "left", pad = "0"))) %>%
  select(TransectNum, Point, o_canopy_percent, OVSPP:PointVisitID) %>%
  mutate(OVABUND = OVABUND %>% as.integer) %>%
  mutate(OVABUND = replace(OVABUND, which(OVABUND == -1), 1)) %>% #***Need to ask about this***
  unique
dat.overstory <- dat.overstory %>%
  mutate(o_canopy_percent = replace(o_canopy_percent, which(Point %in% blanks_more), NA))
#sum(duplicated(dat.overstory %>% select(Point, OVSPP))) # no duplicates found
#check <- dat.overstory %>%
#  group_by(Point) %>%
#  summarise(SumOVABUND = sum(OVABUND)) %>%
#  filter(SumOVABUND != 100)
#write.csv(check, "OVABUND_TotNot100.csv", row.names = F)
dat.overstory2 <- dat.overstory %>%
  filter(!OVSPP == "NULL") %>%
  select(Point, Stratum, o_canopy_percent, primaryHabitat) %>%
  rename(CanCov = o_canopy_percent) %>%
  unique %>%
  # Engelmann spruce cover #
  left_join(dat.overstory %>%
              filter(OVSPP == "ES") %>%
              select(Point, OVABUND) %>%
              rename(RCOV_ES = OVABUND), by = "Point") %>%
  # Subalpine fir cover #
  left_join(dat.overstory %>%
              filter(OVSPP == "SU") %>%
              select(Point, OVABUND) %>%
              rename(RCOV_SU = OVABUND), by = "Point") %>%
  # Pine cover #
  left_join(dat.overstory %>%
              filter(OVSPP %in% c("LP", "LM", "PP", "PY", "BR")) %>%
              select(Point, OVABUND) %>%
              rename(RCOV_Pine = OVABUND), by = "Point") %>%
  # Aspen cover #
  left_join(dat.overstory %>%
              filter(OVSPP == "AS") %>%
              select(Point, OVABUND) %>%
              rename(RCOV_AS = OVABUND), by = "Point") %>%
  # Aspen cover #
  left_join(dat.overstory %>%
              filter(OVSPP == "DF") %>%
              select(Point, OVABUND) %>%
              rename(RCOV_DF = OVABUND), by = "Point") %>%
  # Everything else in other #
  left_join(dat.overstory %>%
              filter(!OVSPP %in% c("AS", "LP", "SU", "ES", "NULL", "LM", "PP", "PY", "BR", "DF")) %>%
              group_by(Point) %>%
              summarise(RCOV_OT = OVABUND %>% sum(na.rm = T)) %>%
              ungroup, by = "Point") %>%
  mutate_at(vars(RCOV_ES:RCOV_OT), funs(replace(.,is.na(.),0))) %>%
  mutate(RCOV_TOT = RCOV_ES + RCOV_SU + RCOV_Pine + RCOV_AS + RCOV_DF + RCOV_OT) %>%
  mutate(RCOV_ES = RCOV_ES/RCOV_TOT, RCOV_SU = RCOV_SU/RCOV_TOT, RCOV_Pine = RCOV_Pine/RCOV_TOT,
         RCOV_AS = RCOV_AS/RCOV_TOT, RCOV_DF = RCOV_DF/RCOV_TOT, RCOV_OT = RCOV_OT/RCOV_TOT) %>%
  select(-RCOV_TOT)

dat.shrub <- read.csv("CPW_veg_query_2_shrub.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>% 
  mutate(Point = str_c(TransectNum, "-", str_pad(Point, width = 2, side = "left", pad = "0"))) %>%
  unique
#sum(duplicated(dat.shrub %>% select(Point, SHSPP))) # 6 found
#pntsWDups <- dat.shrub$Point[which(duplicated(dat.shrub %>% select(Point, SHSPP)))]
#dat.shrub.dups <- dat.shrub %>% filter(Point %in% pntsWDups)
#write.csv(dat.shrub.dups, "Shrub_dups.csv", row.names = F)
# Correct errors (temporary) #
dat.shrub$SHSPP[which(dat.shrub$Point == "CO-CPW-LP-016-09" &
                  dat.shrub$SHSPP == "UD" &
                  dat.shrub$SHABUND == 40)] <- "UC"
dat.shrub$SHSPP[which(dat.shrub$Point == "CO-CPW-LP-031-11" &
                        dat.shrub$SHSPP == "UC" &
                        dat.shrub$SHABUND == 70)] <- "UD"
dat.shrub$SHSPP[which(dat.shrub$Point == "CO-CPW-LP-041-03" &
                        dat.shrub$SHSPP == "UC" &
                        dat.shrub$SHABUND == 99)] <- "UD"
dat.shrub$SHSPP[which(dat.shrub$Point == "CO-CPW-LP-082-11" &
                        dat.shrub$SHSPP == "UC" &
                        dat.shrub$SHABUND == 99)] <- "UD"
dat.shrub$SHSPP[which(dat.shrub$Point == "CO-CPW-SF-021-13" &
                        dat.shrub$SHSPP == "UC" &
                        dat.shrub$SHABUND == 80)] <- "UD"
dat.shrub$SHSPP[which(dat.shrub$Point == "CO-CPW-SF-031-12" &
                        dat.shrub$SHSPP == "UD" &
                        dat.shrub$SHABUND == 30)] <- "UC"
#____________________________#
#check <- dat.shrub %>%
#  filter(SHSPP == "UD") %>%
#  group_by(Point) %>%
#  summarise(UD = sum(SHABUND %>% as.numeric)) %>%
#  full_join(dat.shrub %>%
#              filter(SHSPP == "UC") %>%
#              group_by(Point) %>%
#              summarise(UC = sum(SHABUND %>% as.numeric)), by = "Point")
#check <- check %>%
#  mutate(UC = replace(UC, is.na(UC), 0)) %>%
#  mutate(UD = replace(UD, is.na(UD), 0)) %>%
#  mutate(Tot = UD + UC)
#check <- check %>% filter(Tot != 100)
#write.csv(check, "ShrubAbund_TotNot100.csv", row.names = F)



dat.gcov <- read.csv("CPW_veg_query_2_grndcov.csv", header=T, stringsAsFactors = F) %>%
  tbl_df %>% 
  mutate(Point = str_c(TransectNum, "-", str_pad(Point, width = 2, side = "left", pad = "0"))) %>%
  unique %>%
  select(Point, Year, DATE, PointVisitID, gc_woody:deepAndDown, gc_snow:primaryHabitat) %>%
  mutate(TOT = gc_woody + gc_herb + gc_grass + gc_live_grass +
           gc_bare_litter + deepAndDown + gc_snow + gc_water)
dat.gcov %>% filter(TOT != 100 | is.na(TOT)) %>%
  write.csv("GCov_TotNot100.csv", row.names = F)
