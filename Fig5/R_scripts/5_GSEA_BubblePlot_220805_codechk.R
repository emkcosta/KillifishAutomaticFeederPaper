# Title: Plot GSEA results 
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/GSEA/")
set.seed(1234)

# Load packages
library('ggplot2')
library('dplyr')


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------

# Load Male & Female GOGSEA csv list
# Liver
Input1 <- read.csv('Output/liver_Male_DRoverAL_DEG_220401_GOGSEA.csv', header = T)
Input2 <- read.csv('Output/liver_Female_DRoverAL_DEG_220401_GOGSEA.csv', header = T)


# ------------------------------------------------------------------
# Clean up, Filter, Sort the input data
# ------------------------------------------------------------------

# Extract only the NES and padj terms
Input1 <- Input1[c('ID', 'Description', 'NES', 'p.adjust')]
Input2 <- Input2[c('ID','NES', 'p.adjust')]

# Rename columns to keep track of Male vs. Female data
names(Input1) <- c('ID', 'Description', 'Male_NES', 'Male_p.adjust')
names(Input2) <- c('ID', 'Female_NES', 'Female_p.adjust')


# Make a large dataframe
data <- full_join(Input1, Input2, by = 'ID')


# ------------------------- Define GO terms of interest -----------------------------------------

go <- c('GO:0044391', 'GO:0006415', 'GO:0042254', 'GO:0044272', 'GO:0035384', 'GO:0071616', 'GO:0050778', 'GO:0019221','GO:0002694', 'GO:0006457', 'GO:0034976', 'GO:0000502', 'GO:0006486', 'GO:0048193')

# GO:0044391 ribosomal subunit
# GO:0006415 translational termination
# GO:0042254 ribosome biogenesis
# GO:0044272 sulfur compound biosynthetic process
# GO:0035384 thioester biosynthetic process
# GO:0071616 acyl-CoA biosynthetic process
# GO:0050778 positive regulation of immune response
# GO:0019221 cytokine-mediated signaling pathway
# GO:0002694 regulation of leukocyte activation
# GO:0006457 protein folding
# GO:0034976 response to ER stress
# GO:0000502 proteasome complex
# GO:0006486 protein glycosylation
# GO:0048193 Golgi vesicle transport


# ------------------------------------------------------------------
# Getting the data ready for plotting
# ------------------------------------------------------------------

# Keep only the terms of interest
data_go <- filter(data, data$ID %in% go)

# Make a long list of the data for plotting
data_go_f <- subset(data_go, select = -c(Male_NES, Male_p.adjust))
data_go_m <- subset(data_go, select = -c(Female_NES, Female_p.adjust))
names(data_go_f) <- c('ID', 'Description', 'NES', 'p.adjust')
names(data_go_m) <- c('ID', 'Description', 'NES', 'p.adjust')

# Keep track of the Condition
data_go_f['Condition']<- 'Female'
data_go_m['Condition']<- 'Male'

# Match the order so that ID = my GO list 
df <- data_go_f[match(go, data_go_f$ID), ]
dm <- data_go_m[match(go, data_go_m$ID), ]

# Make the final data frame for plotting
FinalData <- rbind(df, dm)

# Change the level of 'Description' based on the order I defined here.
order <- df$Description
FinalData$Description <- factor(FinalData$Description, levels=rev(order))


# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------

# Liver
pdf('Plots/Liver_GSEA_GO_MvsF_selectedGOterms_220707.pdf', width = 8, height = 10)
ggplot(data=FinalData, aes(x=Condition, y=Description)) +
  geom_point(aes(color=NES, size = -log10(p.adjust))) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  labs(title = 'Liver GSEA by GO terms for males vs. females', x='Conditions', y = '',size='-log(FDR)', color='NES')
dev.off()


sessionInfo() 
