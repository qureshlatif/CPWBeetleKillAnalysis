library(RODBC)
library(stringr)
library(BCRDataAPI)
library(dplyr)
library(stringr)

setwd("/home/RMBO.LOCAL/quresh.latif/CPW_beetle")

#### Species list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('192.168.137.180')
BCRDataAPI::add_columns(c('BCR',
                          'Year',
                          'primaryHabitat',
                          'SelectionMethod',
                          'Migrant',
                          'Stratum',
                          'BirdCode',
                          'Species'
))
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('BirdCode <> NOBI')
BCRDataAPI::filter_on('BCR = 16')
grab <- BCRDataAPI::get_data(interpolate_effort=F)
spp.list <- grab %>%
  mutate(Year = as.integer(as.character(Year))) %>%
  filter(Year %in% 2008:2017) %>%
  mutate(primaryHabitat = as.character(primaryHabitat)) %>%
  filter(primaryHabitat %in% c("LP", "SF", "PP", "MC", "AS", "II")) %>% #II?
  select(primaryHabitat, BirdCode, Species) %>%
  mutate(BirdCode = as.character(BirdCode)) %>%
  mutate(Species = as.character(Species)) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(!Species %in% c("Squirrel, Red", "Ruffed Grouse", "Turkey Vulture", "Wild Turkey",
                         "Sandhill Crane", "Bald Eagle", "American Kestrel", "Red-tailed Hawk",
                         "Great Blue Heron", "Swainson's Hawk", "Canada Goose", "Squirrel, Abert's",
                         "Northern Pygmy-Owl", "Northern Goshawk", "Sharp-shinned Hawk", "Green-winged Teal",
                         "Cooper's Hawk", "Great Horned Owl", "Pika", "Gambel's Quail", "Osprey",
                         "Common Merganser", "White-tailed Ptarmigan", "Peregrine Falcon",
                         "Boreal Owl", "Spotted Owl", "Black-crowned Night-Heron", "Ring-necked Duck",
                         "California Gull", "Northern Saw-whet Owl", "Long-eared Owl", "Flammulated Owl",
                         "Prairie Falcon", "Northern Harrier", "American White Pelican", "Western Screech-Owl",
                         "Double-crested Cormorant", "Bufflehead", "Thicket Tinamou", "Dusky Grouse", "Mallard",
                         "Golden Eagle"))
#spp.all$BirdCode[which(str_detect(spp.all$Species, "Dark-eyed Junco"))]

# Collapsing sub-species and renamed species #
spp.list <- spp.list %>%
  mutate(BirdCode = replace(BirdCode, which(BirdCode %in% c("GHJU","PSJU","RBJU","ORJU")), "DEJU"))
spp.list <- spp.list %>%
  mutate(Species = replace(Species, which(BirdCode == "DEJU"), "Dark-eyed Junco"))
spp.list <- spp.list %>%
  mutate(BirdCode = replace(BirdCode, which(BirdCode == c("WESJ")), "WOSJ"))
spp.list <- spp.list %>%
  mutate(Species = replace(Species, which(BirdCode == "WOSJ"), "Woodhouse's Scrub-Jay"))
spp.list <- spp.list %>%
  mutate(BirdCode = replace(BirdCode, which(BirdCode == c("AUWA")), "YRWA"))
spp.list <- spp.list %>%
  mutate(Species = replace(Species, which(BirdCode == "YRWA"), "Yellow-rumped Warbler"))
spp.list <- spp.list %>%
  mutate(BirdCode = replace(BirdCode, which(BirdCode == c("RSFL")), "NOFL"))
spp.list <- spp.list %>%
  mutate(Species = replace(Species, which(BirdCode == "NOFL"), "Northern Flicker"))
spp.list <- spp.list %>%
  mutate(BirdCode = replace(BirdCode, which(BirdCode == c("MWCS")), "WCSP"))
spp.list <- spp.list %>%
  mutate(Species = replace(Species, which(BirdCode == "WCSP"), "White-crowned Sparrow"))

spp.all <- spp.list %>%
  select(BirdCode, Species) %>%
  unique
spp.red <- spp.list %>%
  filter(primaryHabitat %in% c("LP", "SF", "II")) %>%
  select(BirdCode, Species) %>%
  unique
spp.additional <- spp.all %>%
  filter(!BirdCode %in% spp.red$BirdCode)
