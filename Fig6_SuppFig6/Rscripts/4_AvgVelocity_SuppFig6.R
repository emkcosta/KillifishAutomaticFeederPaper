# Title:  Animal Average Velocities
# Author: Emma K. Costa
# Date: code compiled on 20220822
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022


# ------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------
inputDir <- c('/../1_ProcessingDLCOutPut/Output/velocity/')
outputDir <- c('/../3_SuppFig6/Output/')

setwd(inputDir)

f <- list.files(pattern=c("*.csv"), full.names=F, recursive=FALSE) 

#make one large table
tables <- lapply(f, read.csv)
df <- do.call(rbind, tables)

#add one to the trial #'s so it's no longer 0 indexed
df$recording_num <- df$recording_num + 1


#compute rolling averages
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

# ------------------------------------------------------------------
# Calculate the mean velocities
# ------------------------------------------------------------------
animals.to.include <- c("1", "2", "3", "4", "5", "9","13","14","15", "16", "19", "20", "26")
df <- subset(df, fishnum %in% animals.to.include)

avg.vel <- data.frame('fishnum' = animals.to.include, 'avg.velocity' =  NA)

for(i in 1:length(animals.to.include)){
  animal <- animals.to.include[i] 
  df.subset <- subset(df, fishnum <= animal) #do this on an individual animal basis
  
  upper <- quantile(df.subset$velocity_x_20frame, 0.85, na.rm = T)[[1]] #thresholding top 15%
  lower <- quantile(df.subset$velocity_x_20frame, 0.15, na.rm = T)[[1]] #thresholding bottom 15%
  
  data <- df.subset %>% 
    mutate(velocity_x2 = case_when(
      velocity_x_20frame >= upper ~ upper,
      T ~ velocity_x_20frame
    )) %>%
    mutate(velocity_x3 = case_when(
      velocity_x2 <= lower ~ lower,
      T ~ velocity_x2
    ))
  
  avg.vel[i,]$avg.velocity <- mean(data$velocity_x3, na.rm = T)
}

View(avg.vel)

#save this data
setwd(paste(outputDir))
write.csv(avg.vel, "allanimals_avgvelocities.csv")
