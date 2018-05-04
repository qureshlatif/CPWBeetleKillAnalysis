library(stringr)
library(BCRDataAPI)
library(dplyr)

setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")

################## Inputs ####################
# Get latest AOU checklist with tax names and order #
aou.checklist <- read.csv("/home/RMBO.LOCAL/quresh.latif/NACC_list_bird_species_downloaded_20180319.csv",
                          header = T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(tax_ord = row_number())

spp.exclude <- c("Squirrel, Red", "Ruffed Grouse", "Turkey Vulture", "Wild Turkey",
                 "Sandhill Crane", "Bald Eagle", "American Kestrel", "Red-tailed Hawk",
                 "Great Blue Heron", "Swainson's Hawk", "Canada Goose", "Squirrel, Abert's",
                 "Northern Pygmy-Owl", "Northern Goshawk", "Sharp-shinned Hawk", "Green-winged Teal",
                 "Cooper's Hawk", "Great Horned Owl", "Pika", "Gambel's Quail", "Osprey",
                 "Common Merganser", "White-tailed Ptarmigan", "Peregrine Falcon",
                 "Boreal Owl", "Spotted Owl", "Black-crowned Night-Heron", "Ring-necked Duck",
                 "California Gull", "Northern Saw-whet Owl", "Long-eared Owl", "Flammulated Owl",
                 "Prairie Falcon", "Northern Harrier", "American White Pelican", "Western Screech-Owl",
                 "Double-crested Cormorant", "Bufflehead", "Thicket Tinamou", "Dusky Grouse", "Mallard",
                 "Golden Eagle")
strata <- c("CO-CPW-SF", "CO-CPW-LP")
BCR <- 16
SampDesign <- c("IMBCR", "GRTS")
##############################################

######## Functions ##########
# Collapse species to sub-species #
SubSppToSpp <- function(dat) {
  dat <- dat %>%
    mutate(BirdCode = replace(BirdCode, which(BirdCode %in% c("GHJU","PSJU","RBJU","ORJU")), "DEJU")) %>%
    mutate(Species = replace(Species, which(BirdCode == "DEJU"), "Dark-eyed Junco")) %>%
    mutate(BirdCode = replace(BirdCode, which(BirdCode == c("WESJ")), "WOSJ")) %>%
    mutate(Species = replace(Species, which(BirdCode == "WOSJ"), "Woodhouse's Scrub-Jay")) %>%
    mutate(BirdCode = replace(BirdCode, which(BirdCode == c("AUWA")), "YRWA")) %>%
    mutate(Species = replace(Species, which(BirdCode == "YRWA"), "Yellow-rumped Warbler")) %>%
    mutate(BirdCode = replace(BirdCode, which(BirdCode == c("RSFL")), "NOFL")) %>%
    mutate(Species = replace(Species, which(BirdCode == "NOFL"), "Northern Flicker")) %>%
    mutate(BirdCode = replace(BirdCode, which(BirdCode == c("MWCS")), "WCSP")) %>%
    mutate(Species = replace(Species, which(BirdCode == "WCSP"), "White-crowned Sparrow"))
  return(dat)
}
#############################

#### Compile species list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('Year|int',
                          'primaryHabitat|str',
                          'Stratum|str',
                          'BirdCode|str',
                          'Species|str')
)
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('BirdCode <> NOBI')
BCRDataAPI::filter_on('BCR = 16')
grab <- BCRDataAPI::get_data()

spp.list <- grab %>%
  filter(Year %in% 2008:2017) %>%
  filter(primaryHabitat %in% c("LP", "SF", "II")) %>% #II?
  select(primaryHabitat, BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(!Species %in% spp.exclude)
#spp.all$BirdCode[which(str_detect(spp.all$Species, "Dark-eyed Junco"))]

# Collapsing sub-species and renamed species #
spp.list <- spp.list %>% SubSppToSpp

spp.out <- spp.list %>%
  select(BirdCode, Species) %>%
  unique

#sum(!spp.out$Species %in% aou.checklist$common_name) # check - should be zero
spp.out <- spp.out %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

spp.excluded <- grab %>%
  filter(Year %in% 2008:2017) %>%
  filter(primaryHabitat %in% c("LP", "SF", "II")) %>% #II?
  select(primaryHabitat, BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(Species %in% spp.exclude) %>%
  select(BirdCode, Species) %>%
  unique %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

#### Detection data ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'Date|str',
                          'PointVisitStartTime|str',
                          'centerPointEasting|int',
                          'centerPointNorthing|int',
                          'centerPointZone|int',
                          'Stratum|str',
                          'radialDistance|int',
                          'CL_count|int',
                          'BirdCode|str',
                          'Species|str',
                          'How|str',
                          'Sex|str',
                          'TimePeriod|int'
))

BCRDataAPI::filter_on(str_c('BCR = ', BCR))
BCRDataAPI::filter_on(str_c('Stratum in ', str_c(strata, collapse = ",")))
BCRDataAPI::filter_on(str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")))
BCRDataAPI::filter_on('ninetynine = 0')
BCRDataAPI::filter_on('eightyeight = 0')
BCRDataAPI::filter_on('How <> F')
BCRDataAPI::filter_on('Sex <> J')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('TimePeriod > -1')
BCRDataAPI::filter_on('radialDistance < 125')
grab <- BCRDataAPI::get_data()

point.coords <- grab %>%
  select(TransectNum, Point, centerPointEasting, centerPointNorthing, centerPointZone) %>%
  unique
point.list <- unique(str_c(point.coords$TransectNum,
                           str_pad(point.coords$Point, width = 2, side = "left", pad = "0"), sep = "-")) %>%
  sort
grid.list <- unique(point.coords$TransectNum) %>% sort

## Compare grid IDs to those in JIvan's files to identify missing ##
#JIvan_data <- read.csv("Plot_Level_Covariates.csv", header = T, stringsAsFactors = F)
#JIvan_grids <- JIvan_data$Unit
#JIvan_grids <- str_c("CO-CPW-", str_sub(JIvan_grids, 1, 2), "-", str_sub(JIvan_grids, 3, 5))
#grid.list[which(!grid.list %in% JIvan_grids)]
#JIvan_grids[which(!JIvan_grids %in% grid.list)]

## Point X years surveyed; length is same as for point.list, so not needed ##
#pointXyears.list <- unique(str_c(grab$TransectNum,
#                                 str_pad(grab$Point, width = 2,
#                                         side = "left", pad = "0"),
#                                 grab$Year, sep = "-")) %>%
#  sort

## Add number of detections and count summaries to spp.out by stratum ##
smry <- grab %>% filter(Stratum == "CO-CPW-LP") %>%
  select(BirdCode, TransectNum, Point) %>% unique %>%
  group_by(BirdCode) %>% count() %>%
  rename(Detections_LP = n)
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-SF") %>%
  select(BirdCode, TransectNum, Point) %>% unique %>%
  group_by(BirdCode) %>% count() %>%
  rename(Detections_SF = n)
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-LP") %>%
  select(BirdCode, TransectNum, Point, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount_LP = sum(CL_count))
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-SF") %>%
  select(BirdCode, TransectNum, Point, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount_SF = sum(CL_count))
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

spp.out <- spp.out %>% # replace NAs with zeros
  mutate_at(vars(Detections_LP, Detections_SF, sumCount_LP, sumCount_SF),
            (function(x) replace(x, is.na(x), 0)))

spp.out <- spp.out %>% # compile ratio of count totals to number of detections for spp with > 30 detections #
  mutate(RatioCountToDet_LP = sumCount_LP / Detections_LP) %>%
  mutate(RatioCountToDet_LP = replace(RatioCountToDet_LP, which(Detections_LP < 30), NA)) %>%
  mutate(RatioCountToDet_SF = sumCount_SF / Detections_SF) %>%
  mutate(RatioCountToDet_SF = replace(RatioCountToDet_SF, which(Detections_SF < 30), NA))

maxDetPossible <- str_sub(point.list, 8, 9) %>% tapply(., ., length) # max possible by stratum
names(spp.out)[which(names(spp.out) == "Detections_LP")] <-
  str_c("Detections_LP (max = ",maxDetPossible["LP"],")")
names(spp.out)[which(names(spp.out) == "Detections_SF")] <-
  str_c("Detections_SF (max = ",maxDetPossible["SF"],")")

write.csv(spp.out, "Spp_list.csv", row.names = F)
rm(smry)

## Add number of detections and count summaries to excluded species by stratum ##
smry <- grab %>% filter(Stratum == "CO-CPW-LP") %>%
  select(BirdCode, TransectNum, Point) %>% unique %>%
  group_by(BirdCode) %>% count() %>%
  rename(Detections_LP = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-SF") %>%
  select(BirdCode, TransectNum, Point) %>% unique %>%
  group_by(BirdCode) %>% count() %>%
  rename(Detections_SF = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-LP") %>%
  select(BirdCode, TransectNum, Point, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount_LP = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% filter(Stratum == "CO-CPW-SF") %>%
  select(BirdCode, TransectNum, Point, CL_count) %>%
  group_by(BirdCode) %>%
  summarise(sumCount_SF = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

spp.excluded <- spp.excluded %>% # replace NAs with zeros
  mutate_at(vars(Detections_LP, Detections_SF, sumCount_LP, sumCount_SF),
            (function(x) replace(x, is.na(x), 0)))

write.csv(spp.excluded, "Spp_excluded.csv", row.names = F)
rm(smry)

## Trim and summarize dates ##
grab <- grab %>%
  mutate(Day = Date %>% str_sub(6, 7)) %>%
  mutate(Month = Date %>% str_sub(9, 11)) %>%
  mutate(MM = "05") %>%
  mutate(MM = replace(MM, which(Month == "Jun"), "06")) %>%
  mutate(MM = replace(MM, which(Month == "Jul"), "07")) %>%
  mutate(Date = str_c(MM, Day, sep = "-")) %>%
  select(TransectNum:TimePeriod)
  
#### Habitat covariate data ####
## Database connection ## Run by David. Need ODBC connection on laptop to run this. Can't do SQL on the server.
#library(RODBC)
#con <- odbcConnect("mcb", "RMBO", "f1n4b3b3r")

## SQL query
#qryTxt.1 <-  
#  
#  "
#SELECT     Transect.TransectNum, Point.Point, YEAR(TransectVisit.DATE) AS Year, NewVeg.o_canopy_percent, NewVeg.o_mean_height, NewVeg.NumSnags,
#                     Overstory.Species AS OVSPP, Overstory.SpeciesAbundance AS OVABUND, NewVeg.shrub_cover, ShrubLayer.Species AS SHSPP, 
#ShrubLayer.SpeciesAbundance AS SHABUND, NewVeg.shrub_mean_height, NewVeg.gc_woody, NewVeg.gc_herb, NewVeg.gc_grass, NewVeg.gc_live_grass, 
#NewVeg.gc_bare_litter, NewVeg.deepAndDown, NewVeg.gc_grass_height, NewVeg.gc_live_grass_height, Transect.Stratum, NewVeg.gc_bare, NewVeg.gc_snow, 
#NewVeg.gc_water, NewVeg.primaryHabitat
#FROM         Transect INNER JOIN
#TransectVisit ON Transect.TransectID = TransectVisit.TransectID INNER JOIN
#PointVisit ON TransectVisit.TransectVisitID = PointVisit.TransectVisitID INNER JOIN
#NewVeg ON PointVisit.PointVisitID = NewVeg.PointVisitID INNER JOIN
#Point ON PointVisit.PointID = Point.PointID LEFT OUTER JOIN
#Overstory ON NewVeg.VegID = Overstory.NewVegID LEFT OUTER JOIN
#ShrubLayer ON NewVeg.VegID = ShrubLayer.NewVegID
#WHERE     (Transect.Stratum = N'CO-CPW-SF') OR
#(Transect.Stratum = N'CO-CPW-LP')
#ORDER BY Transect.TransectNum, Year, Point.Point

#"
#composition <- sqlQuery( con, qryTxt.1 )
#odbcClose(con)



BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('TransectNum',
                          'Point',
                          'centerPointEasting',
                          'centerPointNorthing',
                          'centerPointZone',
                          'Stratum',
                          'o_canopy_percent',
                          'NumSnags'
))

BCRDataAPI::filter_on(str_c('Stratum in ', str_c(strata, collapse = ",")))
grab <- BCRDataAPI::get_data(interpolate_effort=T) %>%
  formatInteger(integer.fields = c('Point', 'Year', 'centerPointEasting','centerPointNorthing',
                                   'centerPointZone','radialDistance','TimePeriod'))

## Dave's code ##

# Query master list of grids and points

#BCRDataAPI::reset_api()

#BCRDataAPI::set_api_server('192.168.137.180')

#BCRDataAPI::add_columns(c('TransectNum',
#                          'BCR',
#                          'Stratum',
#                          'SelectionMethod',
#                          'TransectVisitExcludeAnalysis',
#                          'Year',
#                          'Point'
#))

#BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
#BCRDataAPI::filter_on('Year in 2015,2016,2017')
#BCRDataAPI::filter_on('BCR in 18,19,35')
#BCRDataAPI::filter_on('TransectExcludeAnalysis = FALSE')

#occ.sample  <- c()
#occ.sample <- BCRDataAPI::get_data(interpolate_effort=T)
#occ.sample <- occ.sample[,c(1:4,6,7)]
#occ.sample <- unique(occ.sample)

#occ.grid.list <- rbind(overlay,imbcr.grids)
#occ.grid.list <- occ.grid.list[,c(1,2)]

#occ.sample.merge <- merge(occ.grid.list,occ.sample,by=c("TransectNum","Year"),all.x=T)


# Replicate TimePeriod

occ.sample.time.period <- c()

for (i in 1:6) {
  
  prelim <- occ.sample.merge
  
  prelim$TimePeriod <- paste(i,sep="")
  
  occ.sample.time.period <- rbind(occ.sample.time.period,prelim)
  
} 

# Join bird data to Grid, Point and TimePeriod  

occ.sample.spp <- c()

species <- list(as.character(guild.2$BirdCode))

species <- species[[1]][1:74]

for (bird in species) {
  
  BCRDataAPI::reset_api()
  
  BCRDataAPI::set_api_server('192.168.137.180')
  
  BCRDataAPI::add_columns(c('TransectNum',
                            'BCR',
                            'Stratum',
                            'SelectionMethod',
                            'Year',
                            'Point',
                            'radialDistance',
                            'BirdCode',
                            'ninetynine',
                            'eightyeight',
                            'How',
                            'Sex',
                            'Migrant',
                            'TimePeriod'
  ))
  
  BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
  BCRDataAPI::filter_on('Year in 2015,2016,2017')
  BCRDataAPI::filter_on('BCR in 18,19,35')
  BCRDataAPI::filter_on(paste('BirdCode = ',bird,sep=''))
  BCRDataAPI::filter_on('ninetynine = 0')
  BCRDataAPI::filter_on('eightyeight = 0')
  BCRDataAPI::filter_on('BirdCode <> NOBI')
  BCRDataAPI::filter_on('How <> F')
  BCRDataAPI::filter_on('Sex <> J')
  BCRDataAPI::filter_on('Migrant = 0')
  BCRDataAPI::filter_on('radialDistance > -1')
  BCRDataAPI::filter_on('TimePeriod > -1')
  BCRDataAPI::filter_on('radialDistance < 125')
  
  spp.prelim <- c()
  spp.prelim <- BCRDataAPI::get_data()  
  
  occ.prelim <- occ.sample.time.period
  
  occ.prelim$Species <- paste(bird,sep="")
  
  occ.prelim <- merge(occ.prelim,spp.prelim,by=c("TransectNum","Year","BCR","Stratum","SelectionMethod","Point","TimePeriod"),all.x=T)
  
  occ.sample.spp <- rbind(occ.sample.spp,occ.prelim)  
  
}


occ.sample.spp.2 <- unique(occ.sample.spp[,c(1:8,10)])
occ.sample.spp.2$Detect <- ifelse(occ.sample.spp.2$BirdCode=="<NA>",occ.sample.spp.2$Detect==0,1)
occ.sample.spp.2$Detect[is.na(occ.sample.spp.2$Detect)] <- 0
occ.sample.spp.2 <- occ.sample.spp.2[,c(1:8,10)]


# Calculate J minute intervals and y detections by species, point, grid

occ.data.wide <- reshape(occ.sample.spp.2,timevar="TimePeriod",idvar=c("TransectNum","Year","BCR","Stratum","SelectionMethod","Point","Species"),direction="wide",sep="")  
occ.data.wide$t1 <- ifelse(occ.data.wide$Detect1 + occ.data.wide$Detect2 >= 1,1,0)
occ.data.wide$t2 <- ifelse(occ.data.wide$Detect3 + occ.data.wide$Detect4 >= 1,1,0)
occ.data.wide$t3 <- ifelse(occ.data.wide$Detect5 + occ.data.wide$Detect6 >= 1,1,0)
occ.data.wide$J <- ifelse(occ.data.wide$t3>0,3,"")
occ.data.wide$J <- ifelse(occ.data.wide$t2>0,2,occ.data.wide$J)
occ.data.wide$J <- ifelse(occ.data.wide$t1>0,1,occ.data.wide$J)
occ.data.wide$J <- ifelse(occ.data.wide$t1+occ.data.wide$t2+occ.data.wide$t3==0,3,occ.data.wide$J)
occ.data.wide$y <- ifelse(occ.data.wide$t1+occ.data.wide$t2+occ.data.wide$t3==0,0,1)

occ.data.wide.2 <- occ.data.wide[,c(1:7,17,18)]


# Join covariates  

BCRDataAPI::reset_api()

BCRDataAPI::set_api_server('192.168.137.180')

BCRDataAPI::add_columns(c('TransectNum',
                          'BCR',
                          'SelectionMethod',
                          'Year',
                          'Point',
                          'PointVisitStartTime',
                          'primaryHabitat',
                          'o_canopy_percent',
                          'o_mean_height',
                          'shrub_cover',
                          'shrub_mean_height',
                          'gc_woody',
                          'gc_herb',
                          'gc_grass',
                          'gc_live_grass',
                          'gc_grass_height',
                          'gc_live_grass_height',
                          'UnusableDataPrimaryHabitat',
                          'CanopyCover',
                          'CanopyHeight',
                          'ShrubCover',
                          'ShrubHeight',
                          'GroundCover',
                          'GrassHeightResidual',
                          'GrassHeightLive'
))

BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Year in 2015,2016,2017')
BCRDataAPI::filter_on('BCR in 18,19,35')

pt.cov.data <- c()
pt.cov.data <- unique(BCRDataAPI::get_data())  


occ.cov.data <- merge(occ.data.wide.2,pt.cov.data,by=c("TransectNum","Year","BCR","SelectionMethod","Point"),all.x=T)
occ.cov.data$TransectNum <- as.character(occ.cov.data$TransectNum)
occ.cov.data$SelectionMethod <- as.character(occ.cov.data$SelectionMethod)
occ.cov.data$Stratum <- as.character(occ.cov.data$Stratum)
occ.cov.data$Species <- as.character(occ.cov.data$Species)
occ.cov.data$primaryHabitat <- as.character(occ.cov.data$primaryHabitat)
occ.cov.data$UnusableDataPrimaryHabitat <- as.character(occ.cov.data$UnusableDataPrimaryHabitat)
occ.cov.data$CanopyCover <- as.character(occ.cov.data$CanopyCover)
occ.cov.data$CanopyHeight <- as.character(occ.cov.data$CanopyHeight)
occ.cov.data$ShrubCover <- as.character(occ.cov.data$ShrubCover)
occ.cov.data$ShrubHeight <- as.character(occ.cov.data$ShrubHeight)
occ.cov.data$GroundCover <- as.character(occ.cov.data$GroundCover)
occ.cov.data$GrassHeightResidual <- as.character(occ.cov.data$GrassHeightResidual)
occ.cov.data$GrassHeightLive <- as.character(occ.cov.data$GrassHeightLive)

occ.cov.data$Year <- as.integer(as.character(occ.cov.data$Year))
occ.cov.data$BCR <- as.integer(as.character(occ.cov.data$BCR))
occ.cov.data$Point <- as.integer(as.character(occ.cov.data$Point))
occ.cov.data$J <- as.integer(as.character(occ.cov.data$J))
occ.cov.data$y <- as.integer(as.character(occ.cov.data$y))

occ.cov.data$PointVisitStartTime <- as.numeric(as.character(occ.cov.data$PointVisitStartTime))
occ.cov.data$o_canopy_percent <- as.numeric(as.character(occ.cov.data$o_canopy_percent))
occ.cov.data$o_mean_height <- as.numeric(as.character(occ.cov.data$o_mean_height))
occ.cov.data$shrub_cover <- as.numeric(as.character(occ.cov.data$shrub_cover))
occ.cov.data$shrub_mean_height <- as.numeric(as.character(occ.cov.data$shrub_mean_height))
occ.cov.data$gc_woody <- as.numeric(as.character(occ.cov.data$gc_woody))
occ.cov.data$gc_herb <- as.numeric(as.character(occ.cov.data$gc_herb))
occ.cov.data$gc_grass <- as.numeric(as.character(occ.cov.data$gc_grass))
occ.cov.data$gc_live_grass <- as.numeric(as.character(occ.cov.data$gc_live_grass))
occ.cov.data$gc_grass_height <- as.numeric(as.character(occ.cov.data$gc_grass_height))
occ.cov.data$gc_live_grass_height <- as.numeric(as.character(occ.cov.data$gc_live_grass_height))

occ.cov.data$o_canopy_percent <- ifelse(occ.cov.data$o_canopy_percent==-1,0,occ.cov.data$o_canopy_percent)
occ.cov.data$o_mean_height <- ifelse(occ.cov.data$o_mean_height==-1,0,occ.cov.data$o_mean_height)
occ.cov.data$shrub_cover <- ifelse(occ.cov.data$shrub_cover==-1,0,occ.cov.data$shrub_cover)
occ.cov.data$shrub_mean_height <- ifelse(occ.cov.data$shrub_mean_height==-1,0,occ.cov.data$shrub_mean_height)
occ.cov.data$gc_woody <- ifelse(occ.cov.data$gc_woody==-1,0,occ.cov.data$gc_woody)
occ.cov.data$gc_herb <- ifelse(occ.cov.data$gc_herb==-1,0,occ.cov.data$gc_herb)
occ.cov.data$gc_grass <- ifelse(occ.cov.data$gc_grass==-1,0,occ.cov.data$gc_grass)
occ.cov.data$gc_live_grass <- ifelse(occ.cov.data$gc_live_grass==-1,0,occ.cov.data$gc_live_grass)
occ.cov.data$gc_grass_height <- ifelse(occ.cov.data$gc_grass_height==-1,0,occ.cov.data$gc_grass_height)
occ.cov.data$gc_live_grass_height <- ifelse(occ.cov.data$gc_live_grass_height==-1,0,occ.cov.data$gc_live_grass_height)

occ.cov.data$o_canopy_percent <- ifelse(occ.cov.data$CanopyCover==T,NA,occ.cov.data$o_canopy_percent)

occ.cov.data$o_canopy_percent <- ifelse(is.na(occ.cov.data$o_canopy_percent)==T,mean(occ.cov.data$o_canopy_percent, na.rm=T),occ.cov.data$o_canopy_percent)

occ.cov.data$o_mean_height <- ifelse(occ.cov.data$CanopyHeight==T,NA,occ.cov.data$o_mean_height)

occ.cov.data$o_mean_height <- ifelse(is.na(occ.cov.data$o_mean_height)==T,mean(occ.cov.data$o_mean_height, na.rm=T),occ.cov.data$o_mean_height)

occ.cov.data$shrub_cover <- ifelse(occ.cov.data$ShrubCover==T,NA,occ.cov.data$shrub_cover)

occ.cov.data$shrub_cover <- ifelse(is.na(occ.cov.data$shrub_cover)==T,mean(occ.cov.data$shrub_cover, na.rm=T),occ.cov.data$shrub_cover)

occ.cov.data$shrub_mean_height <- ifelse(occ.cov.data$ShrubHeight==T,NA,occ.cov.data$shrub_mean_height)

occ.cov.data$shrub_mean_height <- ifelse(is.na(occ.cov.data$shrub_mean_height)==T,mean(occ.cov.data$shrub_mean_height, na.rm=T),occ.cov.data$shrub_mean_height)

occ.cov.data$gc_woody <- ifelse(occ.cov.data$GroundCover==T,NA,occ.cov.data$gc_woody)

occ.cov.data$gc_woody <- ifelse(is.na(occ.cov.data$gc_woody)==T,mean(occ.cov.data$gc_woody, na.rm=T),occ.cov.data$gc_woody)

occ.cov.data$gc_herb <- ifelse(occ.cov.data$GroundCover==T,NA,occ.cov.data$gc_herb)

occ.cov.data$gc_herb <- ifelse(is.na(occ.cov.data$gc_herb)==T,mean(occ.cov.data$gc_herb, na.rm=T),occ.cov.data$gc_herb)

occ.cov.data$gc_grass <- ifelse(occ.cov.data$GroundCover==T,NA,occ.cov.data$gc_grass)

occ.cov.data$gc_grass <- ifelse(is.na(occ.cov.data$gc_grass)==T,mean(occ.cov.data$gc_grass, na.rm=T),occ.cov.data$gc_grass)

occ.cov.data$gc_live_grass <- ifelse(occ.cov.data$GroundCover==T,NA,occ.cov.data$gc_live_grass)

occ.cov.data$gc_live_grass <- ifelse(is.na(occ.cov.data$gc_live_grass)==T,mean(occ.cov.data$gc_live_grass, na.rm=T),occ.cov.data$gc_live_grass)

occ.cov.data$gc_grass_height <- ifelse(occ.cov.data$GrassHeightResidual==T,NA,occ.cov.data$gc_grass_height)

occ.cov.data$gc_grass_height <- ifelse(is.na(occ.cov.data$gc_grass_height)==T,mean(occ.cov.data$gc_grass_height, na.rm=T),occ.cov.data$gc_grass_height)

occ.cov.data$gc_live_grass_height <- ifelse(occ.cov.data$GrassHeightLive==T,NA,occ.cov.data$gc_live_grass_height)

occ.cov.data$gc_live_grass_height <- ifelse(is.na(occ.cov.data$gc_live_grass_height)==T,mean(occ.cov.data$gc_live_grass_height, na.rm=T),occ.cov.data$gc_live_grass_height)

occ.cov.data.2 <- occ.cov.data[,c(1:21)]

occ.cov.data.2$herb_cover <- occ.cov.data.2$gc_herb + occ.cov.data.2$gc_live_grass

occ.cov.data.2$max_grass_ht <- ifelse(occ.cov.data.2$gc_grass_height>occ.cov.data.2$gc_live_grass_height,occ.cov.data.2$gc_grass_height,occ.cov.data.2$gc_live_grass_height)


# Join w by species  

spp.count <- aggregate(occ.cov.data.2$y,by=list(occ.cov.data.2$Species),FUN="sum")
names(spp.count) <- c("Species","Count")
spp.count$w <- ifelse(spp.count$Count>0,1,0)

names(guild.2)[2] <- "Species"
guild.2$Common.Name <- as.character(guild.2$Common.Name)
guild.2$Guild <- as.character(guild.2$Guild)

spp.list <- merge(guild.2,spp.count,by=c("Species"))
spp.list <- spp.list[,c(1:3,5)]

occ.cov.data.3 <- merge(occ.cov.data.2,spp.list,by=c("Species"),all.x=T)  

# Format StartTime covariate (Turns time into fractional time)

occ.cov.data.5$td2 <- occ.cov.data.5$PointVisitStartTime/100
occ.cov.data.5$td3 <- floor(occ.cov.data.5$td2)
occ.cov.data.5$td4 <- occ.cov.data.5$td2 - occ.cov.data.5$td3
occ.cov.data.5$td5 <- occ.cov.data.5$td4/0.6
occ.cov.data.5$TD <- occ.cov.data.5$td3 + occ.cov.data.5$td5

occ.cov.data.6 <- occ.cov.data.5[,c(1:36,41)]

# Join grid-level covariates (Joins remotely sensed covariates)

grid.bcr18.historic <- read.delim("LEPC_Covariates_BCR18_HISTORIC.csv",header=T,sep=",")

grid.bcr18.usng <- read.delim("LEPC_Covariates_USNG_BCR18_v2.csv",header=T,sep=",")

grid.bcr19.usng <- read.delim("LEPC_Covariates_USNG_BCR19_v2.csv",header=T,sep=",")

grid.bcr35.usng <- read.delim("LEPC_Covariates_USNG_BCR35.csv",header=T,sep=",")

grid.usng <- rbind(grid.bcr18.usng,grid.bcr19.usng,grid.bcr35.usng)

grid.usng.2 <- grid.usng[,c(1:6,8:10)]

grid.cov <- rbind(grid.usng.2,grid.bcr18.historic)

grid.cov.2 <- grid.cov[,c(2:4,7:9)]
names(grid.cov.2)[6] <- "TransectNum"

occ.cov.data.7 <- merge(occ.cov.data.6,grid.cov.2,by=c("TransectNum"),all.x=T)

occ.cov.data.7$Psi.LPCI <- ifelse(occ.cov.data.7$Stratum=="LPCI18-SP"|
                                    occ.cov.data.7$Stratum=="LPCI19-MP"|
                                    occ.cov.data.7$Stratum=="LPCI19-SC"|
                                    occ.cov.data.7$Stratum=="LPCI18-SO",
                                  1,0) 

occ.cov.data.7 <- occ.cov.data.7[order(-occ.cov.data.7$w,
                                       occ.cov.data.7$Guild,
                                       occ.cov.data.7$Species,
                                       occ.cov.data.7$Year,
                                       occ.cov.data.7$TransectNum,
                                       as.numeric(as.character(occ.cov.data.7$Point)
                                       )),]

# Create indices (only needed if we have offsets for guilds) #

occ.cov.data.7$guild.id <- match(occ.cov.data.7$Guild,unique(occ.cov.data.7$Guild))

occ.cov.data.7$spp.id <- match(occ.cov.data.7$Species,unique(occ.cov.data.7$Species))

occ.cov.data.7$yr.id <- match(occ.cov.data.7$Year,unique(occ.cov.data.7$Year))

occ.cov.data.7$yr.grid <- paste(occ.cov.data.7$TransectNum,occ.cov.data.7$Year,sep="-")

occ.cov.data.7$yr.grid.id <- match(occ.cov.data.7$yr.grid,unique(occ.cov.data.7$yr.grid))

occ.cov.data.7$grid.id <- match(occ.cov.data.7$TransectNum,unique(occ.cov.data.7$TransectNum))

occ.cov.data.7$pt.id <- match(occ.cov.data.7$Point,unique(occ.cov.data.7$Point))

occ.cov.data.7$grid.pt <- paste(occ.cov.data.7$TransectNum,occ.cov.data.7$Point,sep="-")

occ.cov.data.7$grid.pt.id <- match(occ.cov.data.7$grid.pt,unique(occ.cov.data.7$grid.pt))

occ.cov.data.7$yr.grid.pt <- paste(occ.cov.data.7$TransectNum,occ.cov.data.7$Year,occ.cov.data.7$Point,sep="-")

occ.cov.data.7$yr.grid.pt.id <- match(occ.cov.data.7$yr.grid.pt,unique(occ.cov.data.7$yr.grid.pt))


# Setup number that corresponds with indices

spp.guild.frame <- unique(occ.cov.data.7[,c("Species","spp.id","Guild","guild.id")])

spp.guild.id <- as.vector(spp.guild.frame$guild.id)

n.spp <- length(unique(occ.cov.data.7$Species))

n.guild <- length(unique(occ.cov.data.7$Guild))

n.yr.grid <- length(unique(occ.cov.data.7$yr.grid.id))

n.yr.grid.pt <- length(unique(occ.cov.data.7$yr.grid.pt.id))

n.rows <- nrow(occ.cov.data.7)

unique(occ.cov.data.7[,c("Species","spp.id")])

# Inits

u.inits <- occ.cov.data.7$y

z.inits.prelim=aggregate(occ.cov.data.7$y,by=list(occ.cov.data.7$spp.id,occ.cov.data.7$yr.grid.id),max)

z.inits=matrix(NA,74,320)
for(i in 1:dim(z.inits.prelim)[1]) {z.inits[z.inits.prelim[i,1],z.inits.prelim[i,2]]=z.inits.prelim[i,3]}

w.inits.prelim <- unique(occ.cov.data.7[,c("Species","w")])

w.inits <- as.vector(w.inits.prelim$w)


save.image("~/Monitoring/LPCH/LEPC_FSA_PF_analysis_2017/Multispp_occupancy/Manuscript analysis/HB_lepc_data_3_mean_test.RData")


load("~/Monitoring/LPCH/LEPC_FSA_PF_analysis_2017/Multispp_occupancy/Manuscript analysis/HB_lepc_data_3_mean_test.RData")

pt1=proc.time()
HB.lepc.spp.rich.1 <- jags(data=list(n.spp=n.spp,n.guild=n.guild,n.yr.grid=n.yr.grid,n.rows=n.rows,y=as.vector(occ.cov.data.7$y),J=as.vector(occ.cov.data.7$J),spp.guild.id=spp.guild.id,yr.grid.id=as.vector(occ.cov.data.7$yr.grid.id),spp.id=as.vector(occ.cov.data.7$spp.id)),
                           inits=list(list("u"=u.inits,"z"=z.inits,"w"=w.inits),
                                      list("u"=u.inits,"z"=z.inits,"w"=w.inits),
                                      list("u"=u.inits,"z"=z.inits,"w"=w.inits)),
                           parameters.to.save=c('p','psi','theta','omega','w','mu.u','rho.ad','rho.ab','rho.bd','sigma.a0','sigma.d0','sigma.b0','mu0','Tau.u','mu0','eta1','tau.eta1','eta2','tau.eta2','eta3','tau.eta3'),
                           model.file='HB_lepc_model_3_mean_test_revised.R',
                           n.iter=500,n.burnin=250,n.thin=3,n.chains=3)
pt2=proc.time()
run.time=pt2-pt1
save(HB.lepc.spp.rich.1,file='HB_lepc_output_3_mean_test_revised_0.5k.rdata')



traplot(HB.lepc.spp.rich.1, parms = c('mu.u[1]'))

traplot(HB.lepc.spp.rich.1, parms = c('mu0[74,1]'))

traplot(HB.lepc.spp.rich.1, parms = c('p[10]'))

traplot(HB.lepc.spp.rich.1, parms = c('eta2[3,1]'))

traplot(HB.lepc.spp.rich.1, parms = c('tau.eta2[8]'))

traplot(HB.lepc.spp.rich.1, parms = c('rho.bd'))

traplot(HB.lepc.spp.rich.1, parms = c('sigma.d0'))

traplot(HB.lepc.spp.rich.1, parms = c('mu0[3,3]'))

traplot(HB.lepc.spp.rich.1, parms = c('omega[2]'))

traplot(HB.lepc.spp.rich.1, parms = c('Tau.u[3,3]'))


print(HB.lepc.spp.rich.1)

str(HB.lepc.spp.rich.1)

# Calculating statistics....... 
# Warning in process.output(samples, DIC = DIC, codaOnly, verbose = verbose) :
#   At least one Rhat value could not be calculated.

# Done. 
# Warning messages:
#   1: In run.model(model.file, data, inits, parameters.to.save, n.chains,  :
#                     Reached max of 10000 adaption iterations; set n.adapt to > 10000
#                   2: In run.model(model.file, data, inits, parameters.to.save, n.chains,  :
#                                     JAGS reports adaptation was # incomplete. Consider increasing n.adapt
#                                   3: In jags.samples(model, variable.names, n.iter, thin, type = "trace",  :
#                                                        Failed to set trace monitor for mu0
#                                                      Monitor already exists and cannot be duplicated


