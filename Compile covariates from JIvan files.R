# 3/30/18
# Create point-level "dead conifer" covariate for BCR analysis starting from "BCR_Dead_Conifer_Qry" in BeetleKill database
#

library(dplyr)
setwd("F:/research stuff/BCR/projects/CPW/JIvan_files")

data <- read.csv("BCR_Dead_Conifer_Qry_Headers.csv", header=T, stringsAsFactors = F) %>%
  tbl_df

plot <- unique(data$RMBOName)

newdata = NULL
#newdataframe <- data.frame(RMBOName = character(), PointNumber = double(), DeadConifPctPt = double()) 

for (i in 1:length(plot)){
  tempplot <- subset(data, RMBOName==plot[i])
  for (j in 1:16){
    temppoint <- subset(tempplot, PointNumber==j)
    temppoint$DeadConif <- temppoint$PctAbundance*0.01*(temppoint$PctRed*0.01 + temppoint$PctSilver*0.01)
    temppoint$DeadConifLP <- temppoint$PctAbundance*0.01*(temppoint$PctRed*0.01 + temppoint$PctSilver*0.01)
    tempdata <- cbind(as.character(plot[i]), j, sum(temppoint$DeadConifPct))
    newdata <- rbind(newdata, tempdata)
  }
}

DeadConiferPctPt <- as.data.frame(newdata, stringsAsFactors=FALSE) 
names(DeadConiferPctPt) <- c("RMBOName", "PointNumber", "DeadConifPct")
DeadConiferPctPt$PointNumber <- as.integer(DeadConiferPctPt$PointNumber)
DeadConiferPctPt$DeadConifPct <- as.numeric(DeadConiferPctPt$DeadConifPct)

write.csv(DeadConiferPctPt, "DeadConiferPctPt.csv", quote=TRUE, row.names=FALSE)
