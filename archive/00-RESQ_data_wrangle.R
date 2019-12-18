library(BCRDataAPI)
library(timeDate)
library(stringr)
library(dplyr)

#setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")
#setwd("/home/RMBO.LOCAL/rstudio03")
setwd("C:/Users/Quresh.Latif/files/projects/CPW")

#### Variables ####
trunc.pct <- 0.95
strata <- c("CO-CPW-SF", "CO-CPW-LP")
nG <- 10 # number of distance categories
###################

# Data grab #
BCRDataAPI::set_api_server('analysis.api.birdconservancy.org')

BCRDataAPI::reset_api()
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Stratum|str',
                          'Year|int',
                          'Date|str',
                          'PointVisitStartTime|str',
                          'easting|int',
                          'northing|int',
                          'zone|int',
                          'BirdCode|str',
                          'radialDistance|int',
                          'TimePeriod|int',
                          'CL_ID|str',
                          'CL_Count|int'))
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on(str_c('Stratum in ', str_c(strata, collapse = ",")))
BCRDataAPI::filter_on('BirdCode in RESQ')
BCRDataAPI::filter_on('TransectVisitExcludeAnalysis = FALSE')
BCRDataAPI::filter_on('TransectExcludeAnalysis = FALSE')
BCRDataAPI::filter_on('ninetynine = 0')
BCRDataAPI::filter_on('eightyeight = 0')
BCRDataAPI::filter_on('How <> F')
BCRDataAPI::filter_on('Sex <> J')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('radialDistance > -1')
BCRDataAPI::filter_on('TimePeriod > -1')
grab <- BCRDataAPI::get_data(interpolate_effort=TRUE)

# Derive parameters for distance sampling #
cutoff <- quantile(grab$radialDistance, trunc.pct, na.rm=TRUE) # truncation distance
area.circle <- as.numeric(pi * (cutoff / 1000) ^ 2) # area of point count circle in km^2
breaks <- seq(0, cutoff, length.out = nG + 1) # breaks for distance categories
area.band <- (pi * breaks[-1]^2) - (pi * breaks[-(nG+1)]^2) # area of each distance category
area.prop <- area.band / sum(area.band)

#***Note: No clusters broken up on separate lines for RESQ, so no processing of clusters needed***
grab.proc <- grab %>%
  mutate(CL_ID = replace(CL_ID, which(radialDistance < 0 | radialDistance >= cutoff), NA)) %>%
  mutate(CL_Count = replace(CL_Count, which(radialDistance < 0 | radialDistance >= cutoff), NA)) %>%
  mutate(radialDistance = replace(radialDistance, which(radialDistance < 0 | radialDistance >= cutoff), NA)) %>%
  mutate(radialDistance = replace(radialDistance, which(radialDistance == 0), 0.01)) %>%
  mutate(dclass = ceiling(radialDistance / breaks[2])) %>%

  ## Trim dates, compile day of year & start time in minutes ##
  mutate(Day = Date %>% str_sub(6, 7)) %>%
  mutate(Month = Date %>% str_sub(9, 11)) %>%
  mutate(MM = "05") %>%
  mutate(MM = replace(MM, which(Month == "Jun"), "06")) %>%
  mutate(MM = replace(MM, which(Month == "Jul"), "07")) %>%
  mutate(Date = str_c(MM, Day, sep = "-")) %>%
  mutate(DOY = Year %>% str_c("-", Date) %>% timeDate() %>% dayOfYear()) %>%
  mutate(PointVisitStartTime = PointVisitStartTime %>%
           replace(which(PointVisitStartTime == "0"), NA)) %>%
  mutate(HR = PointVisitStartTime %>% str_sub(1, -3) %>% as.integer) %>%
  mutate(MIN = PointVisitStartTime %>% str_sub(-2, -1) %>% as.integer) %>%
  mutate(Time_min = HR*60 + MIN) %>%
  select(TransectNum:Date, DOY, Time_min, radialDistance, dclass, TimePeriod, PointVisitStartTime:zone, CL_Count)

## Habitat covariate data ##
cov_tab_import <- read.csv("Covariates.csv", header = T, stringsAsFactors = F) %>%
  tbl_df() %>%
  mutate(YSI = replace(YSI, which(YSI == -1), NA)) %>%
  rename(DeadConif = DeadConif_RCov) %>%
  rename(YSO = YSI) %>%
  rename(HerbCov = gc_herb) %>%
  rename(WoodyCov = gc_woody) %>%
  rename(DDCov = gc_deadAndDown) %>%
  select(Point, TWIP, DeadConif, YSO, CanCov, RCOV_AS, RCOV_ES, RCOV_Pine, shrub_cover,
         RCShrb_UD, HerbCov, WoodyCov, DDCov, WILD) %>%
  filter(!is.na(TWIP))

cov_tab_import <- cov_tab_import %>%
  left_join(foreign::read.dbf("C:/Users/Quresh.Latif/files/GIS/CPW/Point_coords.dbf", as.is = T) %>%
              mutate(Point = str_c(TransectNu, "-", str_pad(Point, width = 2, side = "left", pad = "0"))) %>%
              select(Point, Rd_dens1km, heatload, TWI), by = "Point")

## Compile detection data ##
cov.names <- c("gridIndex", "DayOfYear", "Time", names(cov_tab_import)[-1])
grab.proc <- grab.proc %>%
  mutate(Point = TransectNum %>% str_sub(8, -1) %>%
           str_c("-", (Point %>% str_pad(width = 2, pad = 0, side = "left"))))

# Lodgepole pine stratum #
point.list.LP <- grab.proc$Point[which(str_sub(grab.proc$Point, 1, 2) == "LP")] %>%
  unique %>% sort
point.list.LP <- point.list.LP[which(point.list.LP %in%
                                       (cov_tab_import$Point %>%
                                          str_sub(8, -1)))]

ind.LP <- which(grab.proc$Stratum == "CO-CPW-LP" &
                  grab.proc$Point %in% point.list.LP)
Y.LP.dist <- tapply(grab.proc$CL_Count[ind.LP], grab.proc$Point[ind.LP], sum, na.rm = T)
rm(ind.LP)
dclass.LP <- grab.proc %>%
  filter(!is.na(grab.proc$CL_Count) & grab.proc$Stratum == "CO-CPW-LP") %>%
  filter(Point %in% point.list.LP) %>%
  mutate(PntInd = match(Point, names(Y.LP.dist))) %>%
  select(PntInd, CL_Count, dclass) %>%
  as.matrix()

Y.LP.trem <- matrix(0, nrow = length(point.list.LP), max(grab.proc$TimePeriod, na.rm = T),
                    dimnames = list(point.list.LP, NULL)) # Time removal N-mixture observation matrix
for(k in 1:max(grab.proc$TimePeriod, na.rm = T)) {
  obs <- grab.proc %>% filter(TimePeriod == k & str_sub(Stratum, -2, -1) == "LP") %>%
    filter(Point %in% point.list.LP)
  y <- tapply(obs$CL_Count, obs$Point, sum, na.rm = T)
  Y.LP.trem[names(y), k] <- y
}

Cov.LP <- matrix(NA, nrow = length(point.list.LP), ncol = length(cov.names),
                 dimnames = list(point.list.LP, cov.names))
Cov.LP[, "gridIndex"] <- point.list.LP %>% str_sub(1, -4) %>% as.factor %>% as.integer
Cov.LP[, "DayOfYear"] <- (grab.proc %>% filter(Stratum == "CO-CPW-LP") %>%
                            filter(Point %in% point.list.LP) %>%
                            select(Point, DOY) %>% unique %>% arrange(Point))$DOY
Cov.LP[, "Time"] <- (grab.proc %>% filter(Stratum == "CO-CPW-LP") %>%
                       filter(Point %in% point.list.LP) %>%
                       select(Point, Time_min) %>% unique %>% arrange(Point))$Time_min
ind.vals <- which(point.list.LP %in% (cov_tab_import$Point %>% str_sub(8, -1)))
Cov.LP[ind.vals, -c(1:3)] <- cov_tab_import %>%
  mutate(Point = Point %>% str_sub(8, -1)) %>%
  filter(Point %in% point.list.LP) %>%
  arrange(Point) %>%
  select(-Point) %>%
  as.matrix

# Spruce-fir stratum #
point.list.SF <- grab.proc$Point[which(str_sub(grab.proc$Point, 1, 2) == "SF")] %>%
  unique %>% sort
point.list.SF <- point.list.SF[which(point.list.SF %in%
                                       (cov_tab_import$Point %>%
                                          str_sub(8, -1)))]

ind.SF <- which(grab.proc$Stratum == "CO-CPW-SF" &
                  grab.proc$Point %in% point.list.SF)
Y.SF.dist <- tapply(grab.proc$CL_Count[ind.SF], grab.proc$Point[ind.SF], sum, na.rm = T)
rm(ind.SF)
dclass.SF <- grab.proc %>%
  filter(!is.na(grab.proc$CL_Count) & grab.proc$Stratum == "CO-CPW-SF") %>%
  filter(Point %in% point.list.SF) %>%
  mutate(PntInd = match(Point, names(Y.SF.dist))) %>%
  select(PntInd, CL_Count, dclass) %>%
  as.matrix()

Y.SF.trem <- matrix(0, nrow = length(point.list.SF), max(grab.proc$TimePeriod, na.rm = T),
                    dimnames = list(point.list.SF, NULL)) # Time removal N-mixture observation matrix
for(k in 1:max(grab.proc$TimePeriod, na.rm = T)) {
  obs <- grab.proc %>% filter(TimePeriod == k & str_sub(Stratum, -2, -1) == "SF") %>%
    filter(Point %in% point.list.SF)
  y <- tapply(obs$CL_Count, obs$Point, sum, na.rm = T)
  Y.SF.trem[names(y), k] <- y
}

Cov.SF <- matrix(NA, nrow = length(point.list.SF), ncol = length(cov.names),
                 dimnames = list(point.list.SF, cov.names))
Cov.SF[, "gridIndex"] <- point.list.SF %>% str_sub(1, -4) %>% as.factor %>% as.integer
Cov.SF[, "DayOfYear"] <- (grab.proc %>% filter(Stratum == "CO-CPW-SF") %>%
                            filter(Point %in% point.list.SF) %>%
                            select(Point, DOY) %>% unique %>% arrange(Point))$DOY
Cov.SF[, "Time"] <- (grab.proc %>% filter(Stratum == "CO-CPW-SF") %>%
                       filter(Point %in% point.list.SF) %>%
                       select(Point, Time_min) %>% unique %>% arrange(Point))$Time_min
ind.vals <- which(point.list.SF %in% (cov_tab_import$Point %>% str_sub(8, -1)))
Cov.SF[ind.vals, -c(1:3)] <- cov_tab_import %>%
  mutate(Point = Point %>% str_sub(8, -1)) %>%
  filter(Point %in% point.list.SF) %>%
  arrange(Point) %>%
  select(-Point) %>%
  as.matrix
Cov.SF[which(Cov.SF[, "DeadConif"] == 1.25), "DeadConif"] <- 1 # correct DeadConif > 1

rm(ind.vals, obs, k, y)
save.image("Data_compiled_RESQ.RData")
