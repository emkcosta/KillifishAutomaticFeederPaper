library(ggplot2)

inputDir = '/../SuppFig1/Input/'
dirout = '/../SuppFig1/Output/'

readin = paste(inputDir, "Figure2-SourceData4_missed_feeding_90feeders.csv", sep = "")
df = read.table(readin, sep= ",", header = TRUE)
df = df[df$deviation >= 0,]

numfeeders = length(unique(df$feeder))


# loading in deviations from ideal (3 or 7) feedings
p <- ggplot(data = df,aes(df$deviation)) + geom_histogram() +  theme_classic(base_size=16)
saveout = paste(dirout , "SuppFig1A.pdf", sep ="")
ggsave(saveout)

# create table
total = length(df$deviation)
# getting counts
vals = summary(as.factor(df$deviation))
percent = vals/total







