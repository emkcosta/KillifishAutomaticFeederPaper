# Title:  Conduct Pre-processing on DLC Output
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
#Take in all of the raw DLC trajectories
setwd('/../1_ProcessingDLCOutPut/Input/Training Trajectories/') #this is the unzipped Figure6-SourceData2

#make one large csv file with data from all videos
f <- list.files(path = '.', pattern=c('*csv'), full.names=FALSE, recursive=FALSE) 
tables <- lapply(f, read.csv)

df <- do.call(rbind, tables) #note the videos were processed horizontally, so the x and y coordinates are swapped.

#Remove animals and videos you do not want to include
df <- subset(df, recordingnum <= 16) #only include trials 0-16 (first 17 trials)

animals.to.include <- c("1", "2", "3", "4", "5", "9","13","14","15", "16", "19", "20", "26")
df <- subset(df, fishnum %in% animals.to.include) #only include these animals based on exclusion from manual quantification

nrow(df) #this will give you the # of total frames prefiltering

#save one large CSV with all animal data
setwd('/../1_ProcessingDLCOutPut/Output/prefiltering/')
write.csv(df, "allanimals_allREC_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv")



# ------------------------------------------------------------------
# (1) Remove low likelihood frames < 0.999 likelihood
# ------------------------------------------------------------------
df <- subset(df, fishhead_likelihood > 0.999) 
nrow(df) #this will give you the # of total frames after likelihood thresholding

setwd('/../1_ProcessingDLCOutPut/Output/prefiltering/')
write.csv(df, "allanimals_allREC_likelihood0999_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv")




# ------------------------------------------------------------------
# (2) Differential Smoothing - Eliminate single data points where position of body part moves by a large amount
# ------------------------------------------------------------------
#Find euclidean distance between contiguous frames for all videos (head only)
df$euclidean <- NA

for(a in animals.to.include){
  df.animal <- subset(df,fishnum == a)
  for(rec in 0:20){
    val <- which(df.animal$recordingnum == rec)
    if(!identical(val, integer(0))){ #for some reason, animal 6 had some videos missing, so I needed to account for that
      df.rec <- subset(df.animal,recordingnum == rec)
      nrows <- nrow(df.rec)
      if (nrows > 1){
        for(i in 2:nrows){ #unfortunately you cannot calculate euclidean distance for first one
          r <- df.rec[i,]
          r0 <- df.rec[i-1,]
          framedist <- r$bodyparts_coords - r0$bodyparts_coords
          frame <- r$bodyparts_coords
          if(frame != 0 & framedist == 1){
            deltay <- r$fishhead_y - r0$fishhead_y
            deltax <- r$fishhead_x - r0$fishhead_x
            eu <- sqrt(((deltax)^2) + ((deltay)^2))
            idx <- which((df$fishnum == a) & (df$recordingnum == rec) & (df$bodyparts_coords == frame))
            df[idx,]$euclidean <- eu
          }
        }
      }
    }
  }
}

setwd('/../1_ProcessingDLCOutPut/Output/prefiltering/')
write.csv(df, "220810_allanimals_allREC_likelihood0999_euclidean_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv")

threshold <- quantile(df$euclidean, 0.95, na.rm = T) 

length(which(df$euclidean < threshold)) #this will give you the # of total frames after euclidean distance thresholding

#set a threshold, if the change greater than this threshold, the latter frame was set to NaN
indices <- which(df$euclidean > threshold)
for (i in indices){
  df[i,]$fishhead_x <- NA
  df[i,]$fishhead_y <- NA
}



# ------------------------------------------------------------------
# (3) Interpolation w/ cubic spline
# ------------------------------------------------------------------
#For interpolation, need dataframe with all frame #s so you can do the interpolation, use the original df, but need to set some vals to NA from previous steps
#one way you can do this is by assessing which frames you removed in the previous filtering steps by comparing the OG dataframe to the current version
#if a row from the OG does not exist in the current version, you can set it to NA

setwd('/../1_ProcessingDLCOutPut/Output/prefiltering/')
df.OG <- read.csv('allanimals_allREC_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv')

df.OG$fishhead_x[df.OG$fishhead_likelihood < 0.999] <- NA #Set position values of low likelihood frames (< 0.999 likelihood) to NA
df.OG$fishhead_y[df.OG$fishhead_likelihood < 0.999] <- NA #Set position values of low likelihood frames (< 0.999 likelihood) to NA
df.OG$seconds <- (df.OG$bodyparts_coords)/20

indices <- which(df$euclidean > threshold)    #Set position values to NA in original df for where the euclidean distance is above threshold
for (i in indices){
  r <- df[i,]
  frame <- r$bodyparts_coords
  fishID <- r$fishnum
  rec <- r$recordingnum
  idx <- which((df.OG$bodyparts_coords == frame) & (df.OG$fishnum == fishID) & (df.OG$recordingnum == rec))
  df.OG[idx,]$fishhead_x <- NA
  df.OG[idx,]$fishhead_y <- NA
}

#use na.spline in the zoo package to interpolate across the NAs
fishIDs <- unique(df.OG$fishnum)
totrecnum <- nrow(unique(df.OG[c("fishnum", "recordingnum")]))
recs <- unique(df.OG[c("fishnum", "recordingnum")])

spline_x <- as.list(rep("", totrecnum)) #list of zoo objects
spline_y <- as.list(rep("", totrecnum))
spline_xy <- as.list(rep("", totrecnum))

#calculate interpolation for each recording
for (i in 1:totrecnum){
  r <- recs[i,]
  df.animal <- subset(df.OG,fishnum == r$fishnum)
  df.rec <- subset(df.animal,recordingnum == r$recordingnum)
  
  x1 <- zoo(df.rec$fishhead_x, df.rec$bodyparts_coords)
  y1 <- zoo(df.rec$fishhead_y, df.rec$bodyparts_coords)
  
  x2 <- na.spline(x1) 
  y2 <- na.spline(y1)
  
  spline_x[[i]] <- x2 #each column will be a new recording for each individual
  spline_y[[i]] <- y2 #each column will be a new recording for each individual
  
  x2y2 <- merge(x2, y2)
  spline_xy[[i]] <- x2y2
}

setwd('/../1_ProcessingDLCOutPut/Output/filtered_and_interpolated/')#convert each individual spline to df and save each individual csv, then later you can remerge into one dataframe/csv
for(i in 1:length(spline_xy)){
  if(!is.null(dim(spline_xy[[i]]))){  #some spline indices are empty bc all of the positional values were filtered out due to low likelihood or otherwise
    df <- as.data.frame(spline_xy[[i]])
    df$fishnum <- recs[,1][i]
    df$recording_num <- recs[,2][i]
    df$bodyparts_coords <- as.numeric(row.names(df))
    df$seconds <- (df$bodyparts_coords)/20
    vidInfo <- paste("t",recs[,1][i],"_", "REC00",recs[,2][i], sep = "")
    file <- paste(vidInfo,"DLC_likelihood-euclidean-filtered_splines-interpolated_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv")
    write.csv(df, file)
  }
}



# ------------------------------------------------------------------
# (4) Velocity Calculations
# ------------------------------------------------------------------
setwd('/../1_ProcessingDLCOutPut/Output/filtered_and_interpolated/')
#take in the csvs
f <- list.files(pattern=c("*.csv"), full.names=F, recursive=FALSE) 
tables <- lapply(f, read.csv)

#for each video, perform velocity calculations
for(i in 1:length(f)){
  filename <- f[i]
  df <- tables[[i]]
  df <- subset(df, select = -c(X))
  
  #add in rows with NAs for the missing frames
  frames <- c(0:600) #fill up to 600 frames = 30 seconds
  missing.frames.indices <- which(!frames %in% df$bodyparts_coords) #find the indices of the missing frame #s, make a list
  missing.frames <- missing.frames.indices - 1
  if(!identical(missing.frames, numeric(0))){
    df.toadd <- as.data.frame(missing.frames)
    df.toadd$bodyparts_coords <- missing.frames
    df.toadd$x2 <- NA
    df.toadd$y2 <- NA
    df.toadd$fishnum <- df$fishnum[1]
    df.toadd$recording_num <- df$recording_num[1]
    df.toadd$seconds <- (df.toadd$bodyparts_coords)/20
    df.toadd <- subset(df.toadd, select = -c(missing.frames))
    
    df <- rbind(df, df.toadd) #add rows with these frame numbers
    row.names(df) <- df$bodyparts_coords #rename rows 
    
    df <- df %>% arrange(bodyparts_coords) #reorder the dataframe
  }
  
  #velocity in x (top to bottom)
  df$velocity_x <- NA
  df[1,]$velocity_x <- 0
  for (j in 2:dim(df)[1]){
    df[j,]$velocity_x <- (df[j,]$x2 - df[j-1,]$x2)*20
  }
  
  #velocity in y (left to right)
  df$velocity_y <- NA
  df[1,]$velocity_y <- 0
  for (k in 2:dim(df)[1]){
    df[k,]$velocity_y <- (df[k,]$y2 - df[k-1,]$y2)*20
  }
  
  #velocity in xy   
  df$velocity_xy <- NA
  df[1,]$velocity_xy <- 0
  for (l in 2:dim(df)[1]){
    df[l,]$velocity_xy <- sqrt((df[l,]$x2 - df[l-1,]$x2)^2 + (df[l,]$y2 - df[l-1,]$y2)^2)*20
  }
  
  setwd('/../1_ProcessingDLCOutPut/Output/velocity/')
  file <- paste(strsplit(filename, "interpolated")[[1]][1], "interpolated_kinematics_02",strsplit(filename,"interpolated")[[1]][2], sep = '')
  write.csv(df, file)
}