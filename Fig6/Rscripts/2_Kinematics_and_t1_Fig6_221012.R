# Title:  Animal Trajectory Analysis and t1 Calculations
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
#Dataframe with Red Light and Food Drop values
#using the unfiltered, uninterpolated data as filtering and interpolation were conducted considering only head position
#we want to see all of the red light or food drop position values and their likelihoods
setwd('/../1_ProcessingDLCOutPut/Output/prefiltering/')
df2 <- read.csv('allanimals_allREC_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv')



# ------------------------------------------------------------------
# Figure 6D - Plotting Trajectories for animal t3
# ------------------------------------------------------------------
#Sources:
#https://rpubs.com/JoFrhwld/trajectories
#https://ggplot2.tidyverse.org/reference/geom_path.html
#https://ggplot2.tidyverse.org/reference/geom_segment.html
#https://r-spatial.github.io/sf/articles/sf1.html
#https://stackoverflow.com/questions/60359379/netcdf-display-of-trajectory-data-in-r

setwd('/../1_ProcessingDLCOutPut/Output/velocity/')
f <- list.files(pattern=c("*t3_"), full.names=F, recursive=FALSE) 

#make one large table
tables <- lapply(f, read.csv)
df <- do.call(rbind, tables)

#add one to the trial #'s so it's no longer 0 indexed
df$recording_num <- df$recording_num + 1


#for each video, find the average highest likelihood position for red light
#First 7 trials
#for highest redlight likelihood, calculate the average position and use that as your redlight position
df2.lightsubset <- subset(df2, fishnum == 3 & recordingnum == 4 & redlight_likelihood > 0.999)
avg.redlight.x <- mean(df2.lightsubset$redlight_x, na.rm = T)
avg.redlight.y <- mean(df2.lightsubset$redlight_y, na.rm = T)

df2.foodsubset <- subset(df2, fishnum == 3 & recordingnum == 3 & fooddrop_likelihood > 0.1) #use the food drop location from previous vid bc no high likelihood vals 
avg.fooddrop.x <- mean(df2.foodsubset$fooddrop_x, na.rm = T)
avg.fooddrop.y <- mean(df2.foodsubset$fooddrop_y, na.rm = T)

df.subset <- subset(df, recording_num == 4 & seconds <= 9 & seconds >= 2)

#Plot xy trajectory of a single fish
w <- ggplot(df.subset, aes(x = y2, y = x2))+
  geom_path(aes(colour=seconds), lineend = "round")+
  scale_color_gradientn( colours = c("skyblue", "blue", "red", "orange", "yellow"))+
  scale_x_continuous(breaks = seq(0, 600, by = 50))+
  scale_y_continuous(breaks = seq(0, 700, by = 50), limits = c(50,650)) +
  labs(title = paste("Trajectory of: t3 in trial 4"), x = "x position along tank width", y = "y position")+
  geom_point(aes(x=avg.redlight.y, y=avg.redlight.x), colour="red") +
  geom_rect(mapping = aes(xmin = avg.redlight.y - 10, xmax = avg.redlight.y + 10, ymin = avg.redlight.x - 10, ymax = avg.redlight.x + 10), color = 'red', alpha = 0)+
  geom_point(aes(x=avg.fooddrop.y, y=avg.fooddrop.x), colour="brown") +
  geom_rect(mapping = aes(xmin = avg.fooddrop.y - 10, xmax = avg.fooddrop.y + 10, ymin = avg.fooddrop.x - 10, ymax = avg.fooddrop.x + 10), color = 'brown', alpha = 0)+
  geom_text(aes(x=avg.redlight.y, label="\nred light", y=avg.redlight.x + 30), colour="red", text=element_text(size=11))

w

setwd('/../2_Fig6/Plots/')
pdf("Fig6D_unsuccessful_t3_trial4_trajectory.pdf") 
w
dev.off()


#Last 7 trials
#for highest redlight likelihood, calculate the average position and use that as your redlight position
df2.lightsubset <- subset(df2, fishnum == 3 & recordingnum == 17 & redlight_likelihood > 0.999)
avg.redlight.x <- mean(df2.lightsubset$redlight_x, na.rm = T)
avg.redlight.y <- mean(df2.lightsubset$redlight_y, na.rm = T)

df2.foodsubset <- subset(df2, fishnum == 3 & recordingnum == 17 & fooddrop_likelihood > 0.06)
avg.fooddrop.x <- mean(df2.foodsubset$fooddrop_x, na.rm = T)
avg.fooddrop.y <- mean(df2.foodsubset$fooddrop_y, na.rm = T)

df.subset <- subset(df, recording_num == 17 & seconds <= 9 & seconds >= 2)


#Plot xy trajectory of a single fish
w <- ggplot(df.subset, aes(x = y2, y = x2))+
  geom_path(aes(colour=seconds), lineend = "round")+
  scale_color_gradientn( colours = c("skyblue", "blue", "red", "orange", "yellow"))+
  scale_x_continuous(breaks = seq(0, 600, by = 50))+
  scale_y_continuous(breaks = seq(0, 700, by = 50)) +
  labs(title = paste("Trajectory of: t3 in trial 17"), x = "x position along tank width", y = "y position")+
  geom_point(aes(x=avg.redlight.y, y=avg.redlight.x), colour="red") +
  geom_rect(mapping = aes(xmin = avg.redlight.y - 10, xmax = avg.redlight.y + 10, ymin = avg.redlight.x - 10, ymax = avg.redlight.x + 10), color = 'red', alpha = 0)+
  geom_point(aes(x=avg.fooddrop.y, y=avg.fooddrop.x), colour="brown") +
  geom_rect(mapping = aes(xmin = avg.fooddrop.y - 10, xmax = avg.fooddrop.y + 10, ymin = avg.fooddrop.x - 10, ymax = avg.fooddrop.x + 10), color = 'brown', alpha = 0)+
  geom_text(aes(x=avg.redlight.y, label="\nred light", y=avg.redlight.x + 30), colour="red", text=element_text(size=11))

w

setwd('/../2_Fig6/Plots/')
pdf("Fig6D_successful_t3_trial17_trajectory.pdf") 
w
dev.off()




# ------------------------------------------------------------------
# Figure 6E - Vertical Velocity Heatmap
# ------------------------------------------------------------------
setwd('/../1_ProcessingDLCOutPut/Output/velocity/') #use the outputs from earlier with interpolated trajectories and velocity calculations

animalID <- 't3'            #user input required - do this one-by-one for each animal, as the thresholds may be different for each
f <- list.files(pattern=c(paste("*",animalID,"_" ,sep = "")), full.names=F, recursive=FALSE)            


#make one large table
tables <- lapply(f, read.csv)
df <- do.call(rbind, tables)

#add one to the trial #'s so it's no longer 0 indexed
df$recording_num <- df$recording_num + 1

#compute rolling averages
#NOTE: the videos were processed horizontally, so the x and y in the dataframe are swapped
df <- df %>%
  dplyr::arrange(fishnum) %>% 
  dplyr::group_by(fishnum,recording_num) %>%
  dplyr::mutate(velocity_x_5frame = zoo::rollmean(velocity_x, k = 5, fill = NA),  #5 frames
                velocity_x_10frame = zoo::rollmean(velocity_x, k = 10, fill = NA), #10 frames
                velocity_x_20frame = zoo::rollmean(velocity_x, k = 20, fill = NA), #20 frames
                velocity_x_40frame = zoo::rollmean(velocity_x, k = 40, fill = NA) #40 frames
  )

df <- subset(df, seconds <= 18)
df <- subset(df, recording_num <= 17)


###Plot Heatmap: y velocity, 20 frame rolling average
upper <- quantile(df$velocity_x_20frame, 0.95, na.rm = T)[[1]] #for better visualization on the heatmap, threshold the extremes 
lower <- quantile(df$velocity_x_20frame, 0.05, na.rm = T)[[1]]

data <- df %>% 
  mutate(velocity_x2 = case_when(
    velocity_x_20frame >= upper ~ upper,
    T ~ velocity_x_20frame
  )) %>%
  mutate(velocity_x3 = case_when(
    velocity_x2 <= lower ~ lower,
    T ~ velocity_x2
  ))

p7 <- ggplot(data,aes(recording_num,seconds,fill=velocity_x3)) +
  geom_tile(color= "white",size=0.0001) + 
  scale_fill_gradientn(colors = viridis(n=10)) #NOTE: the ggplot isnt actually taking into account the # of viridis colors, so I need to extract them
p7 <-p7 + scale_y_discrete(limits=1:18)
p7 <-p7 + scale_x_discrete(limits=1:17)
p7 <-p7 + theme_minimal(base_size = 8)
p7 <-p7 + labs(title= paste("Y Velocity - 20 frame rolling average"), x="Trial", y="Seconds")
p7 <-p7 + theme(legend.position = "bottom")+
  theme(plot.title=element_text(size = 14))+
  theme(axis.text.y=element_text(size=6)) +
  theme(strip.background = element_rect(colour="white"))+
  theme(plot.title=element_text(hjust=0))+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=7))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=6))+
  removeGrid()#ggExtra

p7


setwd('/../2_Fig6/Plots/')
pdf(paste("Fig6E_",animalID,"_verticalcomponentvelocity_20framesrollingavg_heatmap.pdf", sep = ''))
p7
dev.off()


# ------------------------------------------------------------------
# Figure 6F,G, I, J - Calculating t1
# ------------------------------------------------------------------
#Note: This part is not fully automated and will require some experimenter discretion
###Distance from Surface--------------------------------------------
rm(df)

###USER INPUT REQUIRED 
animalID <- 't1' 
at_surface <- 150 #user input required - surface threshold
fooddrop_threshold <- 0.9 #user input required - this will require manual inspection of df2
vel.threshold <- 75 #user input required - velocity threshold
heatmap.threshold <- 0.95 #user input require - surface heatmap threshold

inputDir <- c('/../1_ProcessingDLCOutPut/Output/velocity/')
outputDir <- c('/../2_Fig6/Output/')
plotsDir <- c('/../2_Fig6/Plots/')
####


setwd(paste(inputDir))
f <- list.files(pattern=c(paste("*",animalID,"_" ,sep = "")), full.names=F, recursive=FALSE)            

#make one large table
tables <- lapply(f, read.csv)
df <- do.call(rbind, tables)

#add one to the trial #'s so it's no longer 0 indexed
df$recording_num <- df$recording_num + 1

df <- df %>%
  dplyr::arrange(fishnum) %>% 
  dplyr::group_by(fishnum,recording_num) %>%
  dplyr::mutate(velocity_x_5frame = zoo::rollmean(velocity_x, k = 5, fill = NA),  #5 frames
                velocity_x_10frame = zoo::rollmean(velocity_x, k = 10, fill = NA), #10 frames
                velocity_x_20frame = zoo::rollmean(velocity_x, k = 20, fill = NA), #20 frames
                velocity_x_40frame = zoo::rollmean(velocity_x, k = 40, fill = NA) #40 frames
  )


#As before for the red light position, here for finding the average highest likelihood food drop position
#we are using the unfiltered and uninterpolated data (df2)

#for each video, find the average highest likelihood position for food drop
fishIDs <- unique(df2$fishnum)
totrecnum <- nrow(unique(df2[c("fishnum", "recordingnum")]))
recs <- unique(df2[c("fishnum", "recordingnum")])
recs$fooddrop_x <- NA
recs$fooddrop_y <- NA

for(i in 1:nrow(recs)){
  r <- recs[i,]
  df.animal <- subset(df2,fishnum == r$fishnum)
  df.rec <- subset(df.animal,recordingnum == r$recordingnum)
  
  df.rec <- subset(df.rec, fooddrop_likelihood >= fooddrop_threshold)         
  
  recs[i,]$fooddrop_x <- mean(df.rec$fooddrop_x)
  recs[i,]$fooddrop_y <- mean(df.rec$fooddrop_y)
}

#find average position of food drop for this animal --  should be in roughly same position across videos assuming camera position hasnt changed
recs.subset <- subset(recs, fishnum == paste(strsplit(animalID, split = 't')[[1]][2])) #need user input
avg_fooddrop_x <- mean(recs.subset$fooddrop_x, na.rm = T)
avg_fooddrop_y <- mean(recs.subset$fooddrop_y, na.rm = T)

df$surface_dist <- NA

#find the distance of the animal from food drop location at each frame (aka 'surface distance')
for (i in 1:nrow(df)){
  r <- df[i,]$recording_num
  idx <- which(recs.subset$recordingnum == r)
  x <- avg_fooddrop_x
  y <- avg_fooddrop_y
  
  deltay <- df[i,]$y2 - y
  deltax <- df[i,]$x2 - x
  eu <- sqrt(((deltax)^2) + ((deltay)^2))
  
  df[i,]$surface_dist <- eu
}

#save this data
setwd(paste(outputDir))
write.csv(df, paste(animalID,"_allRECs_DLC_likelihood-euc-filtered_splines-int_kinematics_surfacedistance.csv", sep = ''))

rm(df.subset)
df.subset <- subset(df, seconds <= 18 & recording_num <= 17)

#for visualization purposes, threshold these distance values
upper <- quantile(df$surface_dist, paste(heatmap.threshold), na.rm = T)[[1]]

data <- df.subset %>% 
  mutate(surface_dist_2 = case_when(
    surface_dist >= upper ~ upper,
    T ~ surface_dist
  )) 

#axes flipped heatmap of euclidean distance from surface drop location
p1 <- ggplot(data,aes(recording_num,seconds,fill=surface_dist_2)) +
  geom_tile(color= "white",size=0.0001) + 
  scale_fill_gradientn(colors = viridis(n=10), trans = 'reverse') 
p1 <-p1 + scale_y_discrete(limits=1:18)
p1 <-p1 + scale_x_discrete(limits=1:17)
p1 <-p1 + theme_minimal(base_size = 8)
p1 <-p1 + labs(title= paste("Distance from Food Drop Location on Water Surface"), x="Trial", y="Seconds")
p1 <-p1 + theme(legend.position = "bottom")+
  theme(plot.title=element_text(size = 14))+
  theme(axis.text.y=element_text(size=6)) +
  theme(strip.background = element_rect(colour="white"))+
  theme(plot.title=element_text(hjust=0))+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=7))+
  theme(legend.title=element_text(size=8))+
  theme(legend.text=element_text(size=6))+
  removeGrid()#ggExtra

p1

setwd(paste(plotsDir))
pdf(paste(animalID,"_surfacedistance_heatmap.pdf", sep = '')) #can be helpful for visualization, but is not included in the paper
p1
dev.off()

#subset the dataframe to include x2, y2, fishnum, recording_num, bodyparts_coords, seconds, surfacedist, velocity_x_20frame
df3 <- select(df, c('x2', 'y2', 'fishnum', 'recording_num', 'bodyparts_coords', 'seconds', 'surface_dist', 'velocity_x_20frame'))


#find beginning times of regions of closeness to surface - "fish arrival times at surface"
df.dist <- data.frame()

counter <- 0
all.dist <- list()
for(rec in c(1:17)){
  if (length(which(df3$recording_num == rec)) == 0){      # if the video is missing, skip it
    all.dist[[length(all.dist)+1]] <- list(NA)
    next
  }
  df.rec <- subset(df3, recording_num == rec)
  l <- c()
  for(f in 1:nrow(df.rec)){
    dist <- df.rec[f,]$surface_dist
    if(!(is.na(dist))){
      if(dist > at_surface){counter <- 0}
      else{
        if(counter == 1){next}
        else{
          counter <- 1
          l <- c(l,df.rec[f,]$seconds)
        }
      }
    }
  }
  all.dist[[length(all.dist)+1]] <- list(l)
}

all.dist

setwd(paste(outputDir))
capture.output(all.dist, file = paste(animalID,"_distanceburst_indices.txt", sep = ''))  

###Velocity Bursts----------------------------------------------------------
#find the beginning and end of a velocity bursting region

counter <- 0
all.vel.beg <- list()
all.vel.end <- list()

df3.naremove <- na.omit(df3)

for(rec in c(1:17)){
  if (length(which(df3.naremove$recording_num == rec)) == 0){      # if the video is missing, skip it
    all.vel.beg[[length(all.vel.beg)+1]] <- list(NA)
    all.vel.end[[length(all.vel.end)+1]] <- list(NA)
    next
  }
  df.rec <- subset(df3.naremove, recording_num == rec)
  beg <- c()
  end <- c()
  counter <- 0
  for(f in 1:nrow(df.rec)){
    velocity <- df.rec[f,]$velocity_x_20frame
    if(velocity < vel.threshold){
      if(!counter == 0){
        end <- c(end,df.rec[f-1,]$seconds)
        counter <- 0
      }
    }
    else{
      if(counter == 1){next}
      else{
        counter <- 1
        beg <- c(beg,df.rec[f,]$seconds)
      }
    }
  }
  all.vel.beg[[length(all.vel.beg)+1]] <- list(beg)
  all.vel.end[[length(all.vel.end)+1]] <- list(end)
}

all.vel.beg
all.vel.end

setwd(paste(outputDir))
capture.output(all.vel.beg, file = paste(animalID,"_velocityburst_begindices.txt", sep = ''))  
capture.output(all.vel.end, file = paste(animalID,"_velocityburst_endindices.txt", sep = ''))  

###Calculate t1----------------------------------------------------------
#take the index for end of a velocity burst and find the minimum time to distance value
t1s <- data.frame('recording_num' = c(1:17), 't1' =  NA)
t1s$fishnum <- paste(animalID)

for(rec in c(1:17)){
  idx.beg <- all.vel.beg[[rec]][[1]]
  idx.end <- all.vel.end[[rec]][[1]]
  idx.dist <- all.dist[[rec]][[1]]
  
  if(is.null(idx.dist) || is.na(idx.dist) || is.na(idx.beg) || is.na(idx.end) || is.null(idx.end)){
    next
  }
  
  df.rec <- subset(df3.naremove, recording_num == rec)
  
  first_velocity <- idx.beg[1] 
  y <- which(idx.dist < first_velocity)
  if(!(identical(y, integer(0)))){
    idx.dist <- idx.dist[-y]
  }
  first_surface <- idx.dist[1]
  
  #subset the distances to values greater than the minimum velocity index
  
  if(!(length(idx.beg) == length(idx.end))){
    x = min(length(idx.beg), length(idx.end))
    idx.beg <- idx.beg[1:x]
    idx.end <- idx.end[1:x]
  }
  
  
  
  closest_velocity <- NA #hold the end index
  vel_diff <- NA #hold the velocity difference
  
  if(is.null(idx.end)){
    next
  }
  
  for(i in 1:length(idx.end)){
    if(idx.beg[i] < first_surface){
      deltat <- first_surface - idx.beg[i]
      if(!(is.na(vel_diff))){
        if(deltat < vel_diff){
          vel_diff <- deltat
          closest_velocity <- idx.beg[i]
        }
      }
      else{
        vel_diff <- deltat
        closest_velocity <- idx.beg[i]
      }
    }
  }
  t1s[rec,]$t1 <- closest_velocity
}

t1s$t1[t1s$t1 > 18] <- 18  #for t1 values greater than the last second we consider in the video, reset this value to 18seconds

View(t1s)

#save this data
setwd(paste(outputDir))
write.csv(t1s, paste(animalID,"_t1values.csv", sep = ''))


