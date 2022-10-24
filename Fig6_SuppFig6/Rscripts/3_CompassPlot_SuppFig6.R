# Title:  Generate Compass Plots
# Author: Emma K. Costa
# Date: code compiled on 20220822
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(sf)
library(grid)
library(dplyr)
library(knitr)
library(plyr)
library(doMC)
library(writexl)
library(readxl)
library("gridExtra")
library('tidyr')
library(ggExtra)
library(patchwork)
library(RColorBrewer)
library(philentropy)
library(splines)
library(zoo)
library(diffdf)
library("padr")
library(scales)
library('tseries')

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
setwd('/../1_ProcessingDLCOutPut/Output/filtered_and_interpolated/')
#Source: https://stackoverflow.com/questions/39024758/how-to-use-r-package-circular-to-make-rose-plot-of-histogram-data-on-360

#make one large csv file with data from all videos
f <- list.files(pattern=c("*.csv"), full.names=T, recursive=FALSE) 
tables <- lapply(f, read.csv)
df <- do.call(rbind, tables) 

#Remove animals and videos you do not want to include
df$recording_num <- df$recording_num + 1
df <- subset(df, recording_num <= 17) #only include trials 1-17
animals.to.include <- c("1", "2", "3", "4", "5", "9","13","14","15", "16", "19", "20", "26")
df <- subset(df, fishnum %in% animals.to.include) #only include these animals

#save this aggregate CSV somewhere
setwd('/../3_SuppFig6/Output/')
#write.csv(df, "allanimals_REC1-17_fullyfilteredinterpolated_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv")
df <- read.csv('allanimals_REC1-17_fullyfilteredinterpolated_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv')




# ------------------------------------------------------------------
# Get Angles of Animal Trajectory Vectors
# ------------------------------------------------------------------
#get beginning position (t = 2sec) and end position (t = 9sec) for all animals first 7 trials and last 7 trials

#extract the coordinates at 2 secs and 9 secs
beg <- subset(df, seconds == 2)
end <- subset(df, seconds == 9)

#for each, define a vector or a line
fishIDs <- unique(df$fishnum)
totrecnum <- nrow(unique(df[c("fishnum", "recording_num")]))
recs <- unique(df[c("fishnum", "recording_num")])

recs$x1 <- NA
recs$y1 <- NA
recs$x2 <- NA
recs$y2 <- NA

recs$angle <- NA

#find the angle of the vector in 2D plane (convert to polar coord (r, theta) and report theta)
for (i in 1:totrecnum){
  r <- recs[i,]
  fish <- r$fishnum
  vid <- r$recording_num
  x1 <- subset(beg, fishnum == fish & recording_num == vid)$y2
  y1 <- subset(beg, fishnum == fish & recording_num == vid)$x2
  x2 <- subset(end, fishnum == fish & recording_num == vid)$y2
  y2 <- subset(end, fishnum == fish & recording_num == vid)$x2
  
  deltay = y2 - y1
  deltax = x2 - x1
  m = deltay/deltax
  
  recs[i,]$x1 <- y1
  recs[i,]$y1 <- x1
  recs[i,]$x2 <- y2
  recs[i,]$y2 <- x2
  
  if(deltay > 0 & deltax > 0){recs[i,]$angle <- atan(m)*180/pi}
  if(deltay < 0 & deltax < 0){recs[i,]$angle <- 180 + atan(m)*180/pi}
  if(deltay > 0 & deltax < 0){recs[i,]$angle <- 180 + atan(m)*180/pi}
  if(deltay < 0 & deltax > 0){recs[i,]$angle <- 360 + atan(m)*180/pi}
} #note there are a few videos for which no angle can be calculated as t=2 or t=9 may not have values, skip those (might need to do so manually)

#save this aggregate CSV somewhere
setwd('/../3_SuppFig6/Output/')
write.csv(recs, "2secto9sec_angles.csv")

first7 <- subset(recs, recording_num %in% c(1:7))
last7 <- subset(recs, recording_num %in% c(11:17))

setwd('/../3_SuppFig6/Output/')
write.csv(first7, "2secto9sec_angles_first7.csv")
write.csv(last7, "2secto9sec_angles_last7.csv")

first7 <- read.csv('2secto9sec_angles_first7.csv')
last7 <- read.csv('2secto9sec_angles_last7.csv')




# ------------------------------------------------------------------
#  Compass Plot
# ------------------------------------------------------------------
#bin the angles and get frequency
#bin size is 20
first7.freq = hist(first7$angle, breaks=15, include.lowest=TRUE, plot=FALSE)

first7.df <- data.frame(range = first7.freq$mids, frequency = first7.freq$counts)


c1 <- ggplot(data=first7.df,aes(x=range,y=frequency))+
  geom_bar(stat="identity")+ 
  coord_polar(start = -pi/2)+ #sets where 90deg appears on the graph
  scale_x_continuous(breaks = seq(0, 360, 30))+
  ylim(c(0,20))

setwd('/../3_SuppFig6/Plots/')
pdf("220820_compassplot_1st7trials.pdf") 
c1
dev.off()


last7.freq   = hist(last7$angle, breaks=15, include.lowest=TRUE, plot=FALSE)

last7.df <- data.frame(range = last7.freq$mids, frequency = last7.freq$counts)


c2 <- ggplot(data=last7.df,aes(x=range,y=frequency))+
  geom_bar(stat="identity")+ 
  coord_polar(start = -pi/2)+ #sets where 90deg appears on the graph
  scale_x_continuous(breaks = seq(0, 360, 30))+
  ylim(c(0,20))

setwd('/../3_SuppFig6/Plots/')
pdf("220820_compassplot_last7trials.pdf") 
c2
dev.off()
