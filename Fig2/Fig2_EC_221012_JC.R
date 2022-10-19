library(ggplot2)

inputDir = '/../Fig2/Input/'
dirout = "/../Fig2/Output/"


##########################################
# Plotting Figure 2A
##########################################

readin = paste(inputDir,"Figure2-SourceData1_feeding_logs_1feeder.csv", sep = "")
df30days = read.table(readin, sep= ",", header = TRUE)

pdf(file = paste(dirout,'Fig2A.pdf'), width = 6, height = 6)
p <- ggplot(df30days)
p + geom_point(aes(Day, NumberTimesFed)) + ylim(0,8) +  theme_classic(base_size = 16)
dev.off()



##########################################
# Plotting Figure 2B
##########################################

# loading in deviations from ideal (3 or 7) feedings for 41 feeders over 2281 days worth of feedings
readin = paste(inputDir,"Figure2-SourceData2_missed_feeding_41feeders.csv", sep = "")
freqd = read.table(readin, header = TRUE, sep = ",")
p <- ggplot(data = freqd,aes(freqd$MissedFeedings)) + geom_histogram() +  theme_classic(base_size=16)
saveout = paste(dirout,"Fig2B.pdf", sep ="")
ggsave(saveout)



##########################################
# Plotting Figure 2D: pulling in percentage deviation from the mean
##########################################


devf = read.table(paste0(inputDir,"rawdeviations_percent.csv"), header = TRUE, sep = ",")

devf$category = factor(devf$category, levels = c("Multiple People", "Single Person", "Automatic Feeder", "Single Automatic Feeder"))

pdf(file = paste0(dirout,'Fig2D.pdf'), width = 6, height = 6)
ggplot(devf, aes(x = category, y = deviation, color = category))+ geom_jitter( stat="identity", size = 1, width=0.15) + theme_classic() + geom_hline(yintercept = 0) + ylim(-0.8,0.8)
dev.off()






