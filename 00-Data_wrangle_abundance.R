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
Spp_LP <- c("RBNU", "MOCH", "HOWR", "WEWP", "CLNU", "RECR", "PISI", "WAVI",
                  "RCKI", "YRWA", "WETA", "GRAJ", "AMRO", "HETH", "LISP", "DEJU")
Spp_SF <- c("ATTW", "RBNU", "MOCH", "WEWP", "CLNU", "RECR", "PISI", "WAVI",
                  "RCKI", "YRWA", "GRAJ", "STJA", "CORA", "AMRO", "CAFI", "HETH",
                  "CHSP", "LISP", "DEJU")
###################

# Data grab #
BCRDataAPI::set_api_server('192.168.137.180')

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
BCRDataAPI::filter_on(str_c('BirdCode in ', str_c(unique(c(Spp_LP, Spp_SF)), collapse = ",")))
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
any(!c(Spp_LP, Spp_SF) %in% grab$BirdCode) # Should be false

# Derive parameters for distance sampling #
cutoff.LP <- area.circle.LP <- numeric(length = length(Spp_LP))
breaks.LP <- area.band.LP <- area.prop.LP <- rep(list(NULL), length(Spp_LP))
for(i in 1:length(Spp_LP)) {
  dat <- grab %>% filter(BirdCode == Spp_LP[i] & Stratum == "CO-CPW-LP")
  cutoff.LP[i] <- quantile(dat$radialDistance, trunc.pct, na.rm=TRUE) # truncation distance
  area.circle.LP[i] <- as.numeric(pi * (cutoff.LP[i] / 1000) ^ 2) # area of point count circle in km^2
  breaks.LP[[i]] <- seq(0, cutoff.LP[i], length.out = nG + 1) # breaks for distance categories
  area.band.LP[[i]] <- (pi * breaks.LP[[i]][-1]^2) - (pi * breaks.LP[[i]][-(nG+1)]^2) # area of each distance category
  area.prop.LP[[i]] <- area.band.LP[[i]] / sum(area.band.LP[[i]])
}

cutoff.SF <- area.circle.SF <- numeric(length = length(Spp_SF))
breaks.SF <- area.band.SF <- area.prop.SF <- rep(list(NULL), length(Spp_SF))
for(i in 1:length(Spp_SF)) {
  dat <- grab %>% filter(BirdCode == Spp_SF[i] & Stratum == "CO-CPW-SF")
  cutoff.SF[i] <- quantile(dat$radialDistance, trunc.pct, na.rm=TRUE) # truncation distance
  area.circle.SF[i] <- as.numeric(pi * (cutoff.SF[i] / 1000) ^ 2) # area of point count circle in km^2
  breaks.SF[[i]] <- seq(0, cutoff.SF[i], length.out = nG + 1) # breaks for distance categories
  area.band.SF[[i]] <- (pi * breaks.SF[[i]][-1]^2) - (pi * breaks.SF[[i]][-(nG+1)]^2) # area of each distance category
  area.prop.SF[[i]] <- area.band.SF[[i]] / sum(area.band.SF[[i]])
}
rm(dat, i)

## Consolidate clusters ##
#grab %>% filter(!is.na(CL_ID)) %>% # Check within-cluster range of distances
#  dplyr::group_by(TransectNum, Point,
#                  Year, TimePeriod, CL_ID) %>%
#  summarize(DistDiff = max(radialDistance) - min(radialDistance)) %>%
#  View

grab.proc <- grab %>% filter(is.na(CL_ID)) %>%
  select(TransectNum, Point, Year, TimePeriod, BirdCode, Stratum, Date,
         PointVisitStartTime, radialDistance, CL_Count) %>%
  bind_rows(
    grab %>% filter(!is.na(CL_ID)) %>%
      dplyr::group_by(TransectNum, Point, Year,
                    TimePeriod, BirdCode, CL_ID) %>%
      summarize(Stratum = unique(Stratum),
                Date = unique(Date),
                PointVisitStartTime = unique(PointVisitStartTime),
                radialDistance = mean(radialDistance),
                CL_Count = sum(CL_Count)) %>%
      select(-CL_ID)
  )

## Additional processing ##
grab.proc <- grab.proc %>%
  mutate(radialDistance = replace(radialDistance, which(radialDistance == 0), 0.01)) %>%
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
  select(TransectNum:Date, DOY, Time_min, PointVisitStartTime, radialDistance, CL_Count) %>%
  mutate(Point = TransectNum %>% str_sub(8, -1) %>%
           str_c("-", (Point %>% str_pad(width = 2, pad = 0, side = "left"))))

dat.LP <- grab.proc %>% filter(Stratum == "CO-CPW-LP")
dat.SF <- grab.proc %>% filter(Stratum == "CO-CPW-SF")

for(i in 1:length(Spp_LP)) {
  dat.spp <- dat.LP %>%
    filter(BirdCode == Spp_LP[i] & !radialDistance >= cutoff.LP[i]) %>%
    mutate(dclass = ceiling(radialDistance / breaks.LP[[i]][2])) %>%
    select(TransectNum:radialDistance, dclass, CL_Count)
  assign(str_c("detects.", Spp_LP[i], ".LP"), dat.spp)
}

for(i in 1:length(Spp_SF)) {
  dat.spp <- dat.SF %>%
    filter(BirdCode == Spp_SF[i] & !radialDistance >= cutoff.SF[i]) %>%
    mutate(dclass = ceiling(radialDistance / breaks.SF[[i]][2])) %>%
    select(TransectNum:radialDistance, dclass, CL_Count)
  assign(str_c("detects.", Spp_SF[i], ".SF"), dat.spp)
}

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
              select(Point, Rd_dens1km), by = "Point")

## Compile detection data ##
cov.names <- c("gridIndex", "DayOfYear", "Time", names(cov_tab_import)[-1])

# Lodgepole pine stratum #
point.list.LP <- dat.LP$Point %>%
  unique %>% sort
point.list.LP <- point.list.LP[which(point.list.LP %in%
                                       (cov_tab_import$Point %>%
                                          str_sub(8, -1)))]

for(sp in 1:length(Spp_LP)) {
  dat <- eval(as.name(str_c("detects.", Spp_LP[sp], ".LP")))
  dat <- dat %>%
    bind_rows(
      dat.LP %>%
        mutate(BirdCode = NA,
               TimePeriod = NA,
               radialDistance = NA,
               CL_Count = NA) %>%
        unique %>%
        filter(!Point %in% dat$Point)
    ) %>%
    filter(Point %in% point.list.LP)
  Y <- tapply(dat$CL_Count, dat$Point, sum, na.rm = T)
  dclass <- dat %>%
    filter(!is.na(CL_Count)) %>%
    filter(Point %in% point.list.LP) %>%
    select(CL_Count, dclass, Point) %>%
    mutate(Point = match(Point, names(Y)))
  for(i in 2:max(dclass$CL_Count)) { # Replicate rows representing > 1 individual
    dclass <- dclass %>%
      bind_rows(dclass %>%
                  filter(CL_Count == i))
  }
  dclass <- dclass %>%
    select(-CL_Count) %>%
    as.matrix()
  
  assign(str_c("Y.", Spp_LP[sp], ".LP.dist"), Y)
  assign(str_c("dclass.", Spp_LP[sp], ".LP"), dclass)
}

Cov.LP <- matrix(NA, nrow = length(point.list.LP), ncol = length(cov.names),
                 dimnames = list(point.list.LP, cov.names))
Cov.LP[, "gridIndex"] <- point.list.LP %>% str_sub(1, -4) %>% as.factor %>% as.integer
Cov.LP[, "DayOfYear"] <- (dat.LP %>% select(Point, DOY) %>%
                            filter(Point %in% point.list.LP) %>% unique %>% arrange(Point))$DOY
Cov.LP[, "Time"] <- (dat.LP %>% select(Point, Time_min) %>%
                       filter(Point %in% point.list.LP) %>% unique %>% arrange(Point))$Time_min
ind.vals <- which(point.list.LP %in% (cov_tab_import$Point %>% str_sub(8, -1)))
Cov.LP[ind.vals, -c(1:3)] <- cov_tab_import %>%
  mutate(Point = Point %>% str_sub(8, -1)) %>%
  filter(Point %in% point.list.LP) %>%
  arrange(Point) %>%
  select(-Point) %>%
  as.matrix

# Spruce-fir stratum #
point.list.SF <- dat.SF$Point %>%
  unique %>% sort
point.list.SF <- point.list.SF[which(point.list.SF %in%
                                       (cov_tab_import$Point %>%
                                          str_sub(8, -1)))]

for(sp in 1:length(Spp_SF)) {
  dat <- eval(as.name(str_c("detects.", Spp_SF[sp], ".SF")))
  dat <- dat %>%
    bind_rows(
      dat.SF %>%
        mutate(BirdCode = NA,
               TimePeriod = NA,
               radialDistance = NA,
               CL_Count = NA) %>%
        unique %>%
        filter(!Point %in% dat$Point)
    ) %>%
    filter(Point %in% point.list.SF) 
  Y <- tapply(dat$CL_Count, dat$Point, sum, na.rm = T)
  dclass <- dat %>%
    filter(!is.na(CL_Count)) %>%
    filter(Point %in% point.list.SF) %>%
    select(CL_Count, dclass, Point) %>%
    mutate(Point = match(Point, names(Y)))
  for(i in 2:max(dclass$CL_Count)) { # Replicate rows representing > 1 individual
    dclass <- dclass %>%
      bind_rows(dclass %>%
                  filter(CL_Count == i))
  }
  dclass <- dclass %>%
    select(-CL_Count) %>%
    as.matrix()
  
  assign(str_c("Y.", Spp_SF[sp], ".SF.dist"), Y)
  assign(str_c("dclass.", Spp_SF[sp], ".SF"), dclass)
}

Cov.SF <- matrix(NA, nrow = length(point.list.SF), ncol = length(cov.names),
                 dimnames = list(point.list.SF, cov.names))
Cov.SF[, "gridIndex"] <- point.list.SF %>% str_sub(1, -4) %>% as.factor %>% as.integer
Cov.SF[, "DayOfYear"] <- (dat.SF %>% select(Point, DOY) %>%
                            filter(Point %in% point.list.SF) %>% unique %>% arrange(Point))$DOY
Cov.SF[, "Time"] <- (dat.SF %>% select(Point, Time_min) %>%
                       filter(Point %in% point.list.SF) %>% unique %>% arrange(Point))$Time_min
ind.vals <- which(point.list.SF %in% (cov_tab_import$Point %>% str_sub(8, -1)))
Cov.SF[ind.vals, -c(1:3)] <- cov_tab_import %>%
  mutate(Point = Point %>% str_sub(8, -1)) %>%
  filter(Point %in% point.list.SF) %>%
  arrange(Point) %>%
  select(-Point) %>%
  as.matrix
Cov.SF[which(Cov.SF[, "DeadConif"] == 1.25), "DeadConif"] <- 1 # correct DeadConif > 1

rm(ind.vals, i, dat, dat.spp, dclass, Y, grab, grab.proc, sp)
save.image("Data_compiled_abundance.RData")
