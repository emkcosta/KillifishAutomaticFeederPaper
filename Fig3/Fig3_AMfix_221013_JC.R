# Growth rate for AL vs DR
library(ggplot2)
set.seed(1234)

###################################################################
# Functions
###################################################################
# dataset, xaxis, yaxis, colorgroup

# "global" variables
gcol = c(AL = "#0072B5", DR = "#EF932F", OE = "#010101")
pd <- position_jitter(width = 0.1)

# define boxplot function
boxit <- function(df, xa, ya, color, filename, xlabel, ylabel) {
	# getting pval
	pval = wilcox.test(df[,ya] ~ df[,xa])
	
	# real plot, fixed extra point with "outlier.shape = NA" argument
	p = ggplot(df, aes_string(x = xa, y = ya, color = color)) + 
	  geom_boxplot(width=0.5, outlier.shape = NA)  +
	  geom_point(position = pd) + 
	  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() , text=element_text(size=16,family="Helvetica", colour = "#010101"), legend.position ="none", axis.text.x = element_text(colour = "#010101"), axis.text.y = element_text(colour = "#010101")) + 
	  ggtitle(sprintf("p-value = %06f", pval$p.value)) +
	  ylab(ylabel)+ xlab(xlabel) +
	  scale_color_manual(values=gcol)

	print(p)

# filename
saveout = paste("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/GitHub/Fig3/Output/", filename, sep = "")
ggsave(saveout, p)
}


# define getstat function
getstats <- function(dataset, cond1, cond2, stat1) {
x1 = dim(dataset[dataset$FeedingScheme == cond1, ])[1]
print(paste("Number", cond1, ":", x1 ))
y1 = dim(dataset[dataset$FeedingScheme == cond2, ])[1]
print(paste("Number", cond2, ":", y1 ))
	x2 = mean(dataset[dataset$FeedingScheme == cond1, c(stat1)])
print(paste("Mean", cond1, ":", x2 ))
y2 = mean(dataset[dataset$FeedingScheme == cond2, c(stat1)])
print(paste("Mean", cond2, ":", y2 ))
x3 =  median(dataset[dataset$FeedingScheme == cond1, c(stat1)])
print(paste("Median", cond1, ":", x3 ))
y3 = median(dataset[dataset$FeedingScheme == cond2, c(stat1)])
print(paste("Median", cond2, ":", y3 ))
}


###################################################################
# Load growth rate data for DR vs AL
###################################################################

setwd('/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/GitHub/Fig3/Input/')
data = read.table("fig3growthrate_DR.csv", sep = ',', header = TRUE)

data$FeedingScheme = as.factor(data$FeedingScheme)
data$FeedingScheme = factor(data$FeedingScheme, levels = c("AL", "DR"))

# split out males and females
dfm = data[data$Sex == "m", ]
dff = data[data$Sex == "f",]

# select only fish that were observed (1 = observed)
dfm_1 = dfm[dfm$Observed == 1, ]
dff_1 = dff[dff$Observed == 1, ]

# select only fish with lifespan shorter than 120 days (4 months)
dfm_120 = dfm_1[dfm_1$Lifespan_days <= 120, ]
dff_120 = dff_1[dff_1$Lifespan_days <= 120, ]


###################################################################
# Plot all fish growth rates first DR vs AL
###################################################################

getstats(data, "AL", "DR", "GrowthRate")

boxit(data, "FeedingScheme", "GrowthRate", "FeedingScheme", "all_fishsizes_DR.pdf", "Feeding Scheme", "Growth Rate (cm/week)")

###################################################################
# Female fish only DR vs AL
###################################################################

getstats(dff_120, "AL", "DR", "GrowthRate")

boxit(dff_120, "FeedingScheme", "GrowthRate", "FeedingScheme", "female_fishsizes_DR.pdf", "Feeding Scheme", "Growth Rate (cm/week)")


###################################################################
# Male fish only DR vs AL
###################################################################

getstats(dfm_120, "AL", "DR", "GrowthRate")

boxit(dfm_120, "FeedingScheme", "GrowthRate", "FeedingScheme", "male_fishsizes_DR.pdf", "Feeding Scheme", "Growth Rate (cm/week)")

###################################################################
# Growth Rate for OE vs AL
###################################################################

data_OE = read.table("fig3growthrate_OE.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

data_OE = data[data$FeedingScheme != "manual",]

data_OE$FeedingScheme[data_OE$FeedingScheme == "2hr"] = "AL"
data_OE$FeedingScheme[data_OE$FeedingScheme == "1hr"] = "OE"

data_OE$FeedingScheme = as.factor(data_OE$FeedingScheme)

# split out males and females
dfm_OE = data_OE[data_OE$Sex == "male", ]
dff_OE = data_OE[data_OE$Sex == "female",]

# select only fish that were observed (1 = observed)
dfm_1_OE = dfm_OE[dfm_OE$Observed == 1, ]
dff_1_OE = dff_OE[dff_OE$Observed == 1, ]

# select only fish with lifespan shorter than 120 days (4 months)
dfm_120_OE = dfm_1_OE[dfm_1_OE$Lifespan_days <= 120, ]
dff_120_OE = dff_1_OE[dff_1_OE$Lifespan_days <= 120, ]

###################################################################
# Plot all fish growth rates first OE vs AL
###################################################################

#Getting stats
getstats(data_OE, "AL", "OE", "GrowthRate")

boxit(data_OE, "FeedingScheme", "GrowthRate", "FeedingScheme", "all_fishsizes_OE.pdf", "Feeding Scheme", "Growth Rate (cm/week)")

###################################################################
# Female fish only
###################################################################

getstats(dff_120_OE, "AL", "OE", "GrowthRate")

boxit(dff_120_OE, "FeedingScheme", "GrowthRate", "FeedingScheme", "female_fishsizes_OE.pdf", "Feeding Scheme", "Growth Rate (cm/week)")

###################################################################
# Male fish only
###################################################################

getstats(dfm_120_OE, "AL", "OE", "GrowthRate")

boxit(dfm_120_OE, "FeedingScheme", "GrowthRate", "FeedingScheme", "male_fishsizes_OE.pdf", "Feeding Scheme", "Growth Rate (cm/week)")

######################################################################








######################################################################
# Load fertility
######################################################################

df = read.table("fig3_DRfertility.csv", sep = ',', header = TRUE,stringsAsFactors = FALSE)

# set date uncrossed to factor, then order levels
df$Uncrossed = as.factor(df$Uncrossed) 

# swap "7times" for "AL" to match with other data (and color scheme)
df$FeedingScheme[df$FeedingScheme == "7times"] = "AL"

######################################################################
# LiveEmbryos first, DR vs AL
######################################################################

#Getting stats

df[df$Crossed == "6/29/20", c("Pair", "FeedingScheme")]

getstats(df, "AL", "DR", "LiveEmbryosHarvestedDay0")

boxit(df, "FeedingScheme", "LiveEmbryosHarvestedDay0", "FeedingScheme", "DRfertilityLiveEmbryos.pdf", "Feeding Scheme", "Fertility (Fertilized Embryos per mating)")

######################################################################
# Total Embryos, DR vs AL
######################################################################

getstats(df, "AL", "DR", "TotalEmbryosHarvested")

boxit(df, "FeedingScheme", "TotalEmbryosHarvested", "FeedingScheme", "DRfertilityTotalEmbryos.pdf", "Feeding Scheme", "Fertility (Total Embryos per mating)")


######################################################################
# Load fertility data and plot
######################################################################

df = read.table("fig3_OEfertility.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)

# subset to automatic feeding only
df = df[df$FeedingScheme != "manual", ]

df$FeedingScheme[df$FeedingScheme == "7times/AL"] = "AL"
df$FeedingScheme[df$FeedingScheme == "12times/OF"] = "OE"


######################################################################
# Fertilized embryos first
######################################################################

df[df$Crossed == "7/22/19", c("Pair", "FeedingScheme")]

getstats(df, "AL", "OE", "LiveEmbryosHarvestedDay0")

boxit(df, "FeedingScheme", "LiveEmbryosHarvestedDay0", "FeedingScheme", "fertilityOE_LiveEmbryos.pdf", "Feeding Scheme", "Fertility (Fertilized Embryos per mating)")

######################################################################
# Total embryos
######################################################################

getstats(df, "AL", "OE", "TotalEmbryosHarvested")

boxit(df, "FeedingScheme", "TotalEmbryosHarvested", "FeedingScheme", "fertilityOE_TotalEmbryos.pdf", "Feeding Scheme", "Fertility (Total Embryos per mating)")




