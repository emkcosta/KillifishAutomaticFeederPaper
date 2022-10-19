# Title: Plot selected GO term genes as heatmaps
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/GSEA/")
set.seed(1234)

# load package
library(pheatmap)
library(dplyr)

# define color palette
color2<-colorRampPalette(c("blue","white","red"))

# ------------------------------------------------------------------
# Load and clean up data
# ------------------------------------------------------------------

# Load data
nCounts <- read.csv("Input/CountsNormDESeq2_Liver_ALDR_splitMF_220402.csv", header = T)

# Pick male or female GSEA results to plot 
# Female results:
dataHE <- read.csv("Output/liver_Female_DRoverAL_DEG_220401_GOGSEA_HumanName.csv", header = T)
GOterm <- read.csv("Output/liver_Female_DRoverAL_DEG_220401_GOGSEA.csv", header = T)

# Male results:
# dataHE <- read.csv("Output/liver_Male_DRoverAL_DEG_220401_GOGSEA_HumanName.csv", header = T)
# GOterm <- read.csv("Output/liver_Male_DRoverAL_DEG_220401_GOGSEA.csv", header = T)


# rename the killifish gene column to "Gene", matching the nCounts file
colnames(dataHE) <- c('human', "Gene", 'mlog10QvalxFC', 'ENTREZID')

# Join the two input files by the killifish gene name
allCounts <- left_join(dataHE, nCounts, by=c("Gene"))

# Remove columns that are unnecessary
allCounts = subset(allCounts, select = -c(mlog10QvalxFC, ENTREZID))


# ------------------------------------------------------------------
# Make heatmap for scaled Normalized counts
# ------------------------------------------------------------------

GeneGO = NULL

# Select a GO term to plot. Pick one termID at a time (keep the other ones commented out). Note that different GSEA result files are used for different termIDs:

# termID = "GO:0050778" # positive regulation of immune response (I used the male GSEA results to plot)
# termID = "GO:0034976" # response to ER stress (I used the female GSEA results to plot)
termID = "GO:0035384" # thioester biosynthetic process (I used the female GSEA results to plot)

# Extract the gene list for the specific GO term of interest as a data frame (GeneGO)
for (item in 1:length(GOterm$ID)){
  if (GOterm[item,2] == termID){
    GeneGO <- strsplit(GOterm[item,12], split = "/")
    GeneGO <- as.data.frame(GeneGO[[1]])
  }
}

# Add "human' as the column name to the data frame GeneGO, and a dummy tracker value of 1 
names(GeneGO) = "human"
GeneGO$tracker <- 1

# Join the data frames to add all the counts to the gene list, and then remove rows with tracker = NA
allCounts_GeneGO <- full_join(allCounts, GeneGO, by="human")
allCounts_GeneGO <- filter(allCounts_GeneGO, !is.na(allCounts_GeneGO$tracker))

# Remove unnecessary columns for plotting
allCounts_GeneGO_plot = subset(allCounts_GeneGO, select = -c(human, tracker))
allCounts_GeneGO_plot <- unique(allCounts_GeneGO_plot)

# Get the data frame ready for plotting
rownames(allCounts_GeneGO_plot) <- allCounts_GeneGO_plot$Gene
allCounts_GeneGO_plot = subset(allCounts_GeneGO_plot, select = -c(Gene))
allCounts_GeneGO_plot <- allCounts_GeneGO_plot[,sort(names(allCounts_GeneGO_plot))]


# ----------------------------- Plot z-score -----------------------------

# Calculate z-score
allCounts_GeneGO_plot$mean <- rowMeans(allCounts_GeneGO_plot)
allCounts_GeneGO_plot$sd <- apply(allCounts_GeneGO_plot, 1, sd)
allCounts_GeneGO_plot$z_f_AL_L_1 <- ((allCounts_GeneGO_plot$f_AL_L_1)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_f_AL_L_2 <- ((allCounts_GeneGO_plot$f_AL_L_2)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_f_AL_L_3 <- ((allCounts_GeneGO_plot$f_AL_L_3)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_f_AL_L_4 <- ((allCounts_GeneGO_plot$f_AL_L_4)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)

allCounts_GeneGO_plot$z_f_DR_L_1 <- ((allCounts_GeneGO_plot$f_DR_L_1)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_f_DR_L_2 <- ((allCounts_GeneGO_plot$f_DR_L_2)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_f_DR_L_3 <- ((allCounts_GeneGO_plot$f_DR_L_3)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_f_DR_L_4 <- ((allCounts_GeneGO_plot$f_DR_L_4)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)

allCounts_GeneGO_plot$z_m_AL_L_1 <- ((allCounts_GeneGO_plot$m_AL_L_1)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_m_AL_L_2 <- ((allCounts_GeneGO_plot$m_AL_L_2)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_m_AL_L_3 <- ((allCounts_GeneGO_plot$m_AL_L_3)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_m_AL_L_4 <- ((allCounts_GeneGO_plot$m_AL_L_4)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)

allCounts_GeneGO_plot$z_m_DR_L_1 <- ((allCounts_GeneGO_plot$m_DR_L_1)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_m_DR_L_2 <- ((allCounts_GeneGO_plot$m_DR_L_2)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_m_DR_L_3 <- ((allCounts_GeneGO_plot$m_DR_L_3)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)
allCounts_GeneGO_plot$z_m_DR_L_4 <- ((allCounts_GeneGO_plot$m_DR_L_4)-(allCounts_GeneGO_plot$mean))/(allCounts_GeneGO_plot$sd)

# Calculate z-score mean of each condition 
allCounts_GeneGO_plot$z_f_AL_ave <- rowMeans(subset(allCounts_GeneGO_plot, select = c(z_f_AL_L_1, z_f_AL_L_2, z_f_AL_L_3, z_f_AL_L_4)))
allCounts_GeneGO_plot$z_f_DR_ave <- rowMeans(subset(allCounts_GeneGO_plot, select = c(z_f_DR_L_1, z_f_DR_L_2, z_f_DR_L_3, z_f_DR_L_4)))
allCounts_GeneGO_plot$z_m_AL_ave <- rowMeans(subset(allCounts_GeneGO_plot, select = c(z_m_AL_L_1, z_m_AL_L_2, z_m_AL_L_3, z_m_AL_L_4)))
allCounts_GeneGO_plot$z_m_DR_ave <- rowMeans(subset(allCounts_GeneGO_plot, select = c(z_m_DR_L_1, z_m_DR_L_2, z_m_DR_L_3, z_m_DR_L_4)))

# Keep only the average z-scores 
aveL <- subset(allCounts_GeneGO_plot, select = c(z_f_AL_ave, z_f_DR_ave, z_m_AL_ave, z_m_DR_ave))

# Set all values to numeric. This step is needed to make CHUNK function run properly.
aveL[,1] <- as.numeric(aveL[,1])
aveL[,2] <-as.numeric(aveL[,2])
aveL[,3] <-as.numeric(aveL[,3])
aveL[,4] <-as.numeric(aveL[,4])

# create the cluster
CHUNK_aveL <-hclust(as.dist(1-cor(t(aveL))), method="ward.D2") 


# Pick the proper file names based on which termID is used:

# # # For termID = "GO:0050778": positive regulation of immune response
# pdf('Plots/HeatMap_Liver_MvF_M_PosRegImm_AveZscoreCounts_220708.pdf', width = 10, height = 10)
# pheatmap(aveL, main="Heatmap for Positve regulation of immune response (by average z score of normalized counts)", Rowv=as.dendrogram(CHUNK_aveL), clustering_method="ward.D2", treeheight_row=20, treeheight_col=1, cluster_cols = F, cluster_rows = T, fontsize_row=5, col=color2(100), show_colnames=T)
# dev.off()

# For termID = "GO:0034976": response to ER stress
# pdf('Plots/HeatMap_Liver_MvF_F_ERstressResp_AveZscoreCounts_220708.pdf', width = 10, height = 10)
# pheatmap(aveL, main="Heatmap for ER stress response (by average z score of normalized counts)", Rowv=as.dendrogram(CHUNK_aveL), clustering_method="ward.D2", treeheight_row=20, treeheight_col=1, cluster_cols = F, cluster_rows = T, fontsize_row=5, col=color2(100), show_colnames=T)
# dev.off()

# # For termID = "GO:0035384": thioester biosynthetic process 
pdf('Plots/HeatMap_Liver_MvF_F_Thioester_AveZscoreCounts_220708.pdf', width = 10, height = 10)
pheatmap(aveL, main="Heatmap for Thioester biosynthesis (by average z score of normalized counts)", Rowv=as.dendrogram(CHUNK_aveL), clustering_method="ward.D2", treeheight_row=20, treeheight_col=1, cluster_cols = F, cluster_rows = T, fontsize_row=10, col=color2(100), show_colnames=T)
dev.off()


sessionInfo() 

