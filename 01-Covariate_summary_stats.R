require(dplyr)
require(stringr)
setwd("F:/research stuff/BCR/projects/CPW")

dat <- read.csv("Covariates.csv", header = T, stringsAsFactors = F) %>%
  tbl_df() %>%
  mutate(Stratum = Point %>% str_sub(1, 9))

dat.LP <- 