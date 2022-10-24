# Title: PCA and DEseq2 analysis comparing across tissues, sexes, and dietary conditions 
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/sexDEG_enrich")
set.seed(1234)

# load package
library(dplyr)


# ------------------------------------------------------------------
# Load data and find overlaps
# ------------------------------------------------------------------
# Pick either the sexDEGs identified in AL or those identified in DR as input:
# sexDEG <- read.csv("Input/liver_SexDiverged_AL_MoverF_DEG_220426.csv", header = T)
sexDEG <- read.csv("Input/liver_SexDiverged_DR_MoverF_DEG_220714.csv", header = T)

# Load the dietDEGs identified in males or in females
dietDEG_f <- read.csv("Input/liver_Female_DRoverAL_DEG_220401.csv", header = T)
dietDEG_m <- read.csv("Input/liver_Male_DRoverAL_DEG_220401.csv", header = T)

# Find the number of genes overlapped between female dietDEG and male dietDEG
dietDEG_f_sig <- subset(dietDEG_f, padj < 0.05)
dietDEG_m_sig <- subset(dietDEG_m, padj < 0.05)
overlap <- inner_join(dietDEG_m_sig, dietDEG_f_sig, by = 'Gene')
length(overlap$Gene)

# ------------------------------------------------------------------
# Select input to run the following common script 
# ------------------------------------------------------------------
# Pick either male dietDEGs or female dietDEGs as input
# dietDEG <- dietDEG_f
dietDEG <- dietDEG_m

# ------------------------------------------------------------------
# common script for males and females - start here **************
# ------------------------------------------------------------------
# Rename the columns to reflect whether they are from the sexDEG vs. dietDEG dataframes
colnames(sexDEG) <- c('Gene', 'baseMean_s', 'log2FoldChange_s', 'lfcSE_s', 'stat_s', 'pvalue_s', 'padj_s')
colnames(dietDEG) <- c('Gene', 'baseMean_d', 'log2FoldChange_d', 'lfcSE_d', 'stat_d', 'pvalue_d', 'padj_d')

# Keep all genes that are only present in both data frames
join_sexDiet <- inner_join(sexDEG, dietDEG, by = 'Gene')

# Remove the genes with 'NA' for padj (in either sexDEG, dietDEG, or both) 
join_sexDiet_noNA <- subset(join_sexDiet, (!is.na(join_sexDiet$padj_s)) & (!is.na(join_sexDiet$padj_d)))
  # Note: dietDEG_female joined with sexDEG_DR = 16426, dietDEG_male joined with sexDEG_DR = 17732
  # Note: dietDEG_female joined with sexDEG_AL = 16426, dietDEG_male joined with sexDEG_AL = 17732

# Subset based on criteria
sexDEG_sig <- subset(join_sexDiet_noNA, padj_s < 0.05)
sexDEG_notsig <- subset(join_sexDiet_noNA, padj_s >= 0.05)
dietDEG_sig <- subset(join_sexDiet_noNA, padj_d < 0.05)
dietDEG_notsig <- subset(join_sexDiet_noNA, padj_d >= 0.05)

# Find genes that are shared in sexDEG and dietDEG lists
share_notsig <- inner_join(sexDEG_notsig, dietDEG_notsig , by = 'Gene')
share_sig <- inner_join(sexDEG_sig, dietDEG_sig , by = 'Gene')


# ------------------------------ Find the non-dietDEG 'control' gene set ------------------------------
# The control gene set shares the same expression distribution and gene number as dietDEGs

nCounts <- read.csv('Input/CountsNormDESeq2_Liver_ALDR_splitMF_220402.csv', header = T)

# Find the average of normalized counts, using all data
rownames(nCounts) <- nCounts$Gene
nCounts <- subset(nCounts, select = -Gene) # remove Gene column to do math in the next step
nCounts$ave <- rowMeans(nCounts) # find the average counts of all samples
nCounts <- cbind(Gene = rownames(nCounts), nCounts) # add the Gene column back for 'join' function later
nCounts <- subset(nCounts, select = c(Gene, ave)) # clean up dataframe

# Find the average expression for dietDEGs, using all data 
dietDEG_sig_nCounts_all <- left_join(dietDEG_sig, nCounts, by ='Gene') # find counts for dietDEG_sig
dietDEG_sig_nCounts_all <- subset(dietDEG_sig_nCounts_all, select = c(Gene, ave)) # keep only the average column
n <- length(dietDEG_sig$Gene) # calculate the number of genes

nonDietDEG_ref <- data.frame(matrix(ncol = 4, nrow = n)) # create an empty dataframe
colnames(nonDietDEG_ref) <- c('Gene', 'nonDietDEG_RefGenes', 'num_sexDEG', 'num_notSexDEG') # name each column

# Find genes whose expression are +/- 2% of each diet DEG's normalized count
for (i in (1:n)) {
  nCounts_i <- subset(nCounts, (nCounts$ave < (dietDEG_sig_nCounts_all$ave[i])*1.02) & (nCounts$ave > (dietDEG_sig_nCounts_all$ave[i]*0.98)))# within 2 % of expression (note: 1% results in 1 DEG without REF gene)
  nCounts_i <- subset(nCounts_i, Gene!= dietDEG_sig_nCounts_all$Gene[i]) # remove the dietDEG from the dataframe
  nonDietDEG_ref$Gene[i] <- dietDEG_sig_nCounts_all$Gene[i] # save the dietDEG name
  a <- list(nCounts_i$Gene)
  nonDietDEG_ref$nonDietDEG_RefGenes[i] <- a # save the non-dietDEG reference gene list
  nonDietDEG_ref$num_sexDEG[i] <- length((inner_join(nCounts_i,sexDEG_sig, by = 'Gene')$Gene)) # find the overlap between sexDEG and this non-dietDEG reference gene list
  nonDietDEG_ref$num_notSexDEG[i] <- length(nCounts_i$Gene)-nonDietDEG_ref$num_sexDEG[i] # find the non-sexDEG, non-dietDEG list
}

# -------------------------- Sanity checks --------------------------
# sexDEG enrichment in the not-diet DEG reference gene list
sexDEG_ref <- sum(nonDietDEG_ref$num_sexDEG)
notsexDEG_ref <- sum(nonDietDEG_ref$num_notSexDEG)
ratio = sexDEG_ref/(sexDEG_ref+notsexDEG_ref) #female = 2per: 6.3808%; male = 2per: 6.5584%

# find the number of dietDEGs with no reference genes within 2% normalized counts away from their expression
zero_nonDietDEG_ref <- subset(nonDietDEG_ref, (num_sexDEG=='0') &(num_notSexDEG == '0'))
percent_Non_zero <- (n - length(zero_nonDietDEG_ref$Gene))/n  # should be 100% non-zero


# -------------------------- Create non-diet DEG reference gene sets -------------------------- 
nonDietDEG_ref_set <- data.frame(matrix(ncol = 1, nrow = n))
colnames(nonDietDEG_ref_set) <- c('Gene')


# -------Run Fisher's exact test using the non-Diet DEG control gene set ---------
# Simulate the non-dietDEG control gene set for 1000 times, and run Fisher's exact test each time.

# Set up the data table
refset <- data.frame(matrix(ncol=1, nrow=n)) # create a new dataframe
colnames(refset) <- c('Gene') # set column names
ref_table <- data.frame(matrix(ncol=2, nrow=1000)) # create a new dataframe 
colnames(ref_table) <- c('enrichment','pvalue') # set column names

# 1000-times Loop to make the control set, test, and store enrichment & p-value each time. 
for (k in (1:1000)){
  for (i in (1:n)) {
    refset_i <- data.frame(nonDietDEG_ref$nonDietDEG_RefGenes[i]) # for each gene, convert its control gene set from being a list to dataframe
    colnames(refset_i) <- c('Gene') # set column name
    refGene <- sample(refset_i$Gene, 1) # randomly select 1 member from the control gene set 
    refset$Gene[i] <- refGene # save the selected gene. 
    } # At the end of this loop, each gene's control gene constitutes a member of the non-Diet DEG control gene set.
  a <- nrow(inner_join(refset, sexDEG_sig, by = 'Gene'))
  c <- nrow(refset) - a 
  b <- nrow(sexDEG_sig) - a 
  d <- nrow(share_notsig)
  dat <- data.frame ('dietRefDEG_sig' = c (a, c), 'dietRefDEG_notsig' = c(b, d), row.names = c ('sexDEG_sig', 'sexDEG_notsig'), stringsAsFactors = FALSE)
  colnames(dat) <- c('dietRefDEG_sig', 'dietRefDEG_notsig')
  test <- fisher.test(dat)
  ref_table$pvalue[k] <- test$p.value
  ref_table$enrichment[k] <- a/n
}


# ------------------------------------------------------------------
# common script for males and females - end here **************
# ------------------------------------------------------------------



# Below are the codes for exporting data & plotting graphs. Select the proper file name for each run.

# ----------------------------------------------------------------------------------------------
# Export data  
# ----------------------------------------------------------------------------------------------

# Export p-values and enrichment

# write.csv(ref_table, 'Output/FisherExactTest_2per_NonDietDEG_REF_f_220706.csv')
# write.csv(ref_table, 'Output/FisherExactTest_2per_NonDietDEG_REF_m_220706.csv')
# write.csv(ref_table, 'Output/FisherExactTest_2per_DRsexDEG_NonDietDEG_REF_f_220723.csv')
write.csv(ref_table, 'Output/FisherExactTest_2per_DRsexDEG_NonDietDEG_REF_m_220723.csv')



# ---------------------------------------------------------------------------------------------
# Plotting Pi charts
# ---------------------------------------------------------------------------------------------

# *************** PART 1: Observed data *************** 
# Run each section separately for each input chosen above.

# # ------------ The percentage of AL-sexDEGs for the female dietDEGs ------------ 
# both <- nrow(share_sig)
# dietOnly <- nrow(dietDEG_sig) - both 
# slices <- c(both, dietOnly)
# lbls <- c("sex-DEG", "not sex-DEG")
# pct <- round(slices/sum(slices)*100)
# lbls <- paste(lbls, pct) # add percents to labels
# lbls <- paste(lbls,"%",sep="") # ad % to labels
# 
# pdf('Plots/PiChart_Liver_ALsexDEG_dieDEG_enrich_f_220724.pdf')
# pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of AL-sex-DEGs among the liver female diet-DEGs")
# dev.off()

# ------------ The percentage of AL-sexDEGs for the male dietDEGs ------------ 
# both <- nrow(share_sig)
# dietOnly <- nrow(dietDEG_sig) - both
# slices <- c(both, dietOnly)
# lbls <- c("sex-DEG", "not sex-DEG")
# pct <- round(slices/sum(slices)*100)
# lbls <- paste(lbls, pct) # add percents to labels
# lbls <- paste(lbls,"%",sep="") # ad % to labels
# 
# pdf('Plots/PiChart_Liver_ALsexDEG_dieDEG_enrich_m_220724.pdf')
# pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of AL-sex-DEGs among the liver male diet-DEGs")
# dev.off()

# # ------------ The percentage of DR-sexDEGs for the female dietDEGs ------------ 
# both <- nrow(share_sig)
# dietOnly <- nrow(dietDEG_sig) - both
# slices <- c(both, dietOnly)
# lbls <- c("sex-DEG", "not sex-DEG")
# pct <- round(slices/sum(slices)*100)
# lbls <- paste(lbls, pct) # add percents to labels
# lbls <- paste(lbls,"%",sep="") # ad % to labels
# 
# pdf('Plots/PiChart_Liver_DRsexDEG_dieDEG_enrich_f_220724.pdf')
# pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of DR-sex-DEGs among the liver female diet-DEGs")
# dev.off()
# 
# # ------------ The percentage of DR-sexDEGs for the male dietDEGs ------------ 
both <- nrow(share_sig)
dietOnly <- nrow(dietDEG_sig) - both
slices <- c(both, dietOnly)
lbls <- c("sex-DEG", "not sex-DEG")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels

pdf('Plots/PiChart_Liver_DRsexDEG_dieDEG_enrich_m_220724.pdf')
pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of DR-sex-DEGs among the liver male diet-DEGs")
dev.off()
# 

# *************** PART 2: Simulated control data *************** 

# load data
# ref_table_al_f <- read.csv('Output/FisherExactTest_2per_NonDietDEG_REF_f_220706.csv', row.names = 1)
# ref_table_al_m <- read.csv('Output/FisherExactTest_2per_NonDietDEG_REF_m_220706.csv', row.names = 1)
# ref_table_dr_f <- read.csv('Output/FisherExactTest_2per_DRsexDEG_NonDietDEG_REF_f_220723.csv', row.names = 1)
ref_table_dr_m <- read.csv('Output/FisherExactTest_2per_DRsexDEG_NonDietDEG_REF_m_220723.csv', row.names = 1)

# select data to run the common script in the next step
# table <- ref_table_al_f
# table <- ref_table_al_m
# table <- ref_table_dr_f
table <- ref_table_dr_m

# ---- common script starts here **** -----

both <- median(table$enrichment)
dietOnly <- 1 - both
slices <- c(both, dietOnly)
lbls <- c("sex-DEG", "not sex-DEG")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels

# ---- common script ends here **** -----

# Pick the section to run based on which input 'table' was selected in Line 248.

# Plotting
# pdf('Plots/PiChart_Liver_ALsexDEG_dieDEG_enrich_ctrl_f_220724.pdf')
# pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of AL-sex-DEGs among the liver female non-diet-DEGs controls")
# dev.off()

# pdf('Plots/PiChart_Liver_ALsexDEG_dieDEG_enrich_ctrl_m_220724.pdf')
# pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of AL-sex-DEGs among the liver male non-diet-DEGs controls")
# dev.off()

# pdf('Plots/PiChart_Liver_DRsexDEG_dieDEG_enrich_ctrl_f_220724.pdf')
# pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of DR-sex-DEGs among the liver female non-diet-DEGs controls")
# dev.off()
# 
pdf('Plots/PiChart_Liver_DRsexDEG_dieDEG_enrich_ctrl_m_220724.pdf')
pie(slices,labels = lbls, col=c("gray","white"), main="Percentage of DR-sex-DEGs among the liver male non-diet-DEGs controls")
dev.off()
# 

sessionInfo() 
