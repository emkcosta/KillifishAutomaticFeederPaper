# Title: PCA and DEseq2 analysis comparing across tissues, sexes, and dietary conditions 
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
 
# Set wd to the current directory
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/DEseq2")

# Please create the following folders in the directory if not present:
# 'Plots', and then sub-folders named 'PCA', 'MAplot',
# 'Output'

# Load packages
library("DESeq2")
library("ggplot2")
library("dplyr")
library("PCAtools")

# Define color variables for plotting
c17 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075')
c20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------

# Load experiment design and count matrix
sampleTable <- read.csv("Input/DRexp20211219ExperimentDesign.csv", row.names = 1)  
lib_name <- as.character(sampleTable$lib)
sampleTable

countdata = read.csv("Input/Counts_ALDR_220325.csv", header = T, row.names = 1)
colnames(countdata) <- c(lib_name)

coldata <- DataFrame(sampleTable)
rownames(coldata)
colnames(countdata) == rownames(coldata) # Make sure that the column names are identical

# Add a pseudocount of 1 to all the Count values
countdata_1 = countdata + 1


# ------------------------------------------------------------------
# LIVER SECTION
# ------------------------------------------------------------------
# Subset for Liver samples only
sampleTable_liver <- sampleTable[which(sampleTable$tissue=='L'), ]
sample_liver <- as.character(sampleTable_liver$lib)
countdata_1_liver <- countdata_1[, c(sample_liver)]

# To make sure that the column names are identical
coldata_liver <- DataFrame(sampleTable_liver)
colnames(countdata_1_liver) == rownames(coldata_liver)

# Make dds object for liver DESEQ2
dds_liver <- DESeqDataSetFromMatrix(countData = countdata_1_liver,
                              colData = coldata_liver,
                              design = ~ Condition)

# only keep rows that have at least 15 sum count
dds_liver <- dds_liver[ rowSums(counts(dds_liver)) > 16, ]



#---------------- Plot PCA for liver samples only ----------------
# Plot PCA using vst normalization
vsd_liver <- vst(dds_liver)
head(assay(vsd_liver), 3)

# Plot other PC combinations; related Figure 5B and Supplemental figure 5B
p <- pca(assay(vsd_liver), metadata = colData(dds_liver))
pdf("Plots/PCA/ALDR_PCA_liveronly_PC1PC2_220331.pdf", width = 10, height = 10)
biplot(p,
       x = 'PC1', y = 'PC2',
       lab = NA,
       colby = 'Condition',
       shape = 'sex',
       legendPosition = 'top')
dev.off()

p <- pca(assay(vsd_liver), metadata = colData(dds_liver))
pdf("Plots/PCA/ALDR_PCA_liveronly_PC1PC3_220331.pdf", width = 10, height = 10)
biplot(p,
       x = 'PC1', y = 'PC3',
       lab = NA,
       colby = 'Condition',
       shape = 'sex',
       legendPosition = 'top')
dev.off()


# ----------- DESeq (split male & female) analysis for liver samples -----------
# Run DEseq
dds_liver <- DESeq(dds_liver)

# Normalized counts
dds_liver <- estimateSizeFactors(dds_liver)
normcounts <- counts(dds_liver, normalized=TRUE)
nrow(normcounts)
colnames(normcounts) <- rownames(sampleTable_liver)
write.table(normcounts, file="Output/CountsNormDESeq2_Liver_ALDR_splitMF_220331.csv", sep= ",")
head(normcounts)

# extract results for male, comparing DR over AL (contrast parameters = c(factor to compare, numerator, denominator))
liver_Male_DRoverAL <- results (dds_liver, contrast = c("Condition", "m_DR_L", "m_AL_L"))
write.table(liver_Male_DRoverAL, file="Output/liver_Male_DRoverAL_DEG_220401.txt", sep="\t")

# extract results for female, comparing DR over AL (contrast parameters = c(factor to compare, numerator, denominator))
liver_Female_DRoverAL <- results (dds_liver, contrast = c("Condition", "f_DR_L", "f_AL_L"))
write.table(liver_Female_DRoverAL, file="Output/liver_Female_DRoverAL_DEG_220401.txt", sep="\t")

# extract results to compare male-AL over female-AL (contrast parameters = c(factor to compare, numerator, denominator))
liver_mAL_over_fAL <- results (dds_liver, contrast = c("Condition", "m_AL_L", "f_AL_L"))
write.table(liver_mAL_over_fAL, file="Output/liver_SexDiverged_AL_MoverF_DEG_220426.txt", sep="\t")

# extract results to compare male-DR over female-DR (contrast parameters = c(factor to compare, numerator, denominator))
liver_mAL_over_fDR <- results (dds_liver, contrast = c("Condition", "m_DR_L", "f_DR_L"))
write.table(liver_mAL_over_fDR, file="Output/liver_SexDiverged_DR_MoverF_DEG_220716.txt", sep="\t")


# ----------- Plot MA Plots ----------- 
# related to Figure 5-figure supplement

# Plot MA plots for the results extracted for males, comparing DR over AL
data_liver_Male_DRoverAL <-plotMA(liver_Male_DRoverAL, alpha = 0.05, colSig = 'Red', main = 'MA plot for Male DR over AL in liver', xlab = 'mean of normalized counts', cex = 0.8, returnData = TRUE)
pdf("Plots/MAplot/liver_Male_DRoverAL_MAplot_220620_1.pdf", width = 8, height = 5)
ggplot(data_liver_Male_DRoverAL, aes(x=log10(mean), y=lfc, color=isDE)) +
  geom_point(data=subset(data_liver_Male_DRoverAL, isDE==FALSE), color="grey", alpha=0.2) +
  geom_point(data=subset(data_liver_Male_DRoverAL, isDE==TRUE), color="red") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylab("log2 fold change") +
  xlab("log10 (mean of normalized counts)")
dev.off()

# Plot MA plots for the results extracted for females, comparing DR over AL
data_liver_Female_DRoverAL <-plotMA(liver_Female_DRoverAL, alpha = 0.05, colSig = 'Red', main = 'MA plot for Female DR over AL in liver', xlab = 'mean of normalized counts', cex = 0.8, returnData = TRUE)

pdf("Plots/MAplot/liver_Female_DRoverAL_MAplot_220620_1.pdf", width = 8, height = 5)
ggplot(data_liver_Female_DRoverAL, aes(x=log10(mean), y=lfc, color=isDE)) +
  geom_point(data=subset(data_liver_Female_DRoverAL, isDE==FALSE), color="grey", alpha=0.2) +
  geom_point(data=subset(data_liver_Female_DRoverAL, isDE==TRUE), color="red") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylab("log2 fold change") +
  xlab("log10 (mean of normalized counts)")
dev.off()

# Plot MA plots for the results extracted for sex-DEGs in AL, comparing male-AL over female-AL
data_liver_mAL_over_fAL <- plotMA(liver_mAL_over_fAL, alpha = 0.05, colSig = 'Red', main = 'MA plot for Male AL over Female AL in liver', xlab = 'mean of normalized counts', cex = 0.8, returnData = TRUE)
pdf("Plots/MAplot/liver_MaleALoverFemaleAL_MAplot_220620_1.pdf", width = 8, height = 5)
ggplot(data_liver_mAL_over_fAL, aes(x=log10(mean), y=lfc, color=isDE)) +
  geom_point(data=subset(data_liver_mAL_over_fAL, isDE==FALSE), color="grey", alpha=0.2) +
  geom_point(data=subset(data_liver_mAL_over_fAL, isDE==TRUE), color="red") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylab("log2 fold change") +
  xlab("log10 (mean of normalized counts)")
dev.off()

# Plot MA plots for the results extracted for sex-DEGs in DR, comparing male-DR over female-DR
data_liver_mAL_over_fDR <- plotMA(liver_mAL_over_fDR, alpha = 0.05, colSig = 'Red', main = 'MA plot for Male DR over Female DR in liver', xlab = 'mean of normalized counts', cex = 0.8, returnData = TRUE)
pdf("Plots/MAplot/liver_MaleDRoverFemaleDR_MAplot_220724.pdf", width = 8, height = 5)
ggplot(data_liver_mAL_over_fDR, aes(x=log10(mean), y=lfc, color=isDE)) +
  geom_point(data=subset(data_liver_mAL_over_fDR, isDE==FALSE), color="grey", alpha=0.2) +
  geom_point(data=subset(data_liver_mAL_over_fDR, isDE==TRUE), color="red") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylab("log2 fold change") +
  xlab("log10 (mean of normalized counts)")
dev.off()



# ----- DESeq (Interaction) analysis for liver samples -----
# Make dds object for liver DESEQ2
dds_liverDEG <- DESeqDataSetFromMatrix(countData = countdata_1_liver,
                                    colData = coldata_liver,
                                    design = ~ sex + feeding + sex:feeding)

# Only keep rows that have at least 15 sum count
dds_liverDEG <- dds_liverDEG[ rowSums(counts(dds_liverDEG)) > 16, ]

# Relevel to specify the reference to be "female_AL"
dds_liverDEG$sex <- relevel(dds_liverDEG$sex, ref = "f")
dds_liverDEG$feeding <- relevel(dds_liverDEG$feeding, ref = "AL")

# Run DEseq
dds_liverDEG <- DESeq(dds_liverDEG)

# Normalized counts
dds_liverDEG <- estimateSizeFactors(dds_liverDEG)
normcounts <- counts(dds_liverDEG, normalized=TRUE)
nrow(normcounts)
colnames(normcounts) <- rownames(sampleTable_liver)
write.table(normcounts, file="Output/CountsNormDESeq2_Liver_ALDR_interaction_220331.csv", sep= ",")
head(normcounts)

liver_Allcomparison_DEG <- results (dds_liverDEG, name = "sexm.feedingDR")
write.table(liver_Allcomparison_DEG, file="Output/liver_ALDR_interaction_DEG_220331.txt", sep="\t")



# ------------------------------------------------------------------
# BRAIN SECTION
# ------------------------------------------------------------------
# Subset for Brain samples only
sampleTable_brain <- sampleTable[which(sampleTable$tissue=='BR'), ]
sample_brain <- as.character(sampleTable_brain$lib)
countdata_1_brain <- countdata_1[, c(sample_brain)]

# To make sure that the column names are identical
coldata_brain <- DataFrame(sampleTable_brain)
colnames(countdata_1_brain) == rownames(coldata_brain)

# Make dds object for brain DESEQ2
dds_brain <- DESeqDataSetFromMatrix(countData = countdata_1_brain,
                                    colData = coldata_brain,
                                    design = ~ Condition)

# only keep rows that have at least 15 sum count
dds_brain <- dds_brain[ rowSums(counts(dds_brain)) > 16, ]


# ---------- Remove outliers ----------
# reference: https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.Rmd

# vst normalization
vsd_brain <- vst(dds_brain)
head(assay(vsd_brain), 3)

rv <- rowVars(assay(vsd_brain))
pc <- prcomp(t(assay(vsd_brain)[head(order(-rv),1000),]))
plot(pc$x[,1:2])
pc1 <- pc$x[ ,1]
zscore_f_AL_BR_4 = (pc1[14] - mean(pc1))/sd(pc1) #2.985454, ~3 standard deviation from the mean
zscore_f_DR_BR_1 = (pc1[2] - mean(pc1))/sd(pc1) #1.231851
zscore_m_DR_BR_2 = (pc1[6] - mean(pc1))/sd(pc1) #1.301088

idx <- pc$x[,1] < 50 #2.6059 SD away from the mean
sum(idx)
plot(pc$x[,1:2], col=idx+1, pch=20, asp=1)


# Remove outlier
countdata_1_brain_noOutliner <-countdata_1_brain[idx]
dds_brain_noOutliner <- dds_brain[,idx]
sampleTable_brain_noOutlier <- sampleTable_brain[idx,]


# ------------------------------------------------------------------
# Plot PCA for brain samples without outlier
# ------------------------------------------------------------------
# related to Figure 5B

# --------------Plot PCA using vst normalization --------------------
vsd_brain_noOutliner <- vst(dds_brain_noOutliner)
head(assay(vsd_brain_noOutliner), 3)

p1 <- pca(assay(vsd_brain_noOutliner), metadata = colData(vsd_brain_noOutliner))
pdf("Plots/PCA/ALDR_Brain_noOutlier_PC1PC2_bySex_220512.pdf", width = 10, height = 10)
biplot(p1,
       x = 'PC1', y = 'PC2',
       lab = NA,
       colby = 'Condition',
       shape = 'sex',
       legendPosition = 'top')
dev.off()

p2 <- pca(assay(vsd_brain_noOutliner), metadata = colData(vsd_brain_noOutliner))
pdf("Plots/PCA/ALDR_Brain_noOutlier_PC1PC3_bySex_220512.pdf", width = 10, height = 10)
biplot(p1,
       x = 'PC1', y = 'PC3',
       lab = NA,
       colby = 'Condition',
       shape = 'sex',
       legendPosition = 'top')
dev.off()

sessionInfo() 

