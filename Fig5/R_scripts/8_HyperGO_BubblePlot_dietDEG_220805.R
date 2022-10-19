# Title: Plot Hypergeometric GO analysis results
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/HyperGO")
set.seed(1234)

# Load packages
library('ggplot2')
library('dplyr')


# ------------------------------------------------------------------
# Load data, clean up
# ------------------------------------------------------------------
# dietDEG in males
m_al <- read.csv('Input/GOBP_dietDEG_m_AL_UP_220714.csv', header = TRUE) # male
m_dr <- read.csv('Input/GOBP_dietDEG_m_DR_UP_220714.csv', header = TRUE) # male

# dietDEG in females
f_al <- read.csv('Input/GOBP_dietDEG_f_AL_UP_220714.csv', header = TRUE) # female
f_dr <- read.csv('Input/GOBP_dietDEG_f_DR_UP_220714.csv', header = TRUE) # female

# Extract only the enrichment and padj terms
m_al <- m_al[c('GOBPID', 'Term', 'enrichment', 'padj')]
m_dr <- m_dr[c('GOBPID', 'Term', 'enrichment', 'padj')]
f_al <- f_al[c('GOBPID', 'Term', 'enrichment', 'padj')]
f_dr <- f_dr[c('GOBPID', 'Term', 'enrichment', 'padj')]

# Make the GO term names as row names
rownames(m_al) <- m_al$GOBPID 
rownames(m_dr) <- m_dr$GOBPID
rownames(f_al) <- f_al$GOBPID 
rownames(f_dr) <- f_dr$GOBPID

# ------------------------------------------------------------------
# Subset specific GO terms, define plotting order
# ------------------------------------------------------------------
# Subset specific GO terms
m_al_go <- m_al[c('GO:0015850', 'GO:0046486', 'GO:0006644', 'GO:0008610', 'GO:0009987'), ]
m_dr_go <- m_dr[c('GO:0050708', 'GO:0010817', 'GO:0044283', 'GO:0043066', 'GO:0006725'), ]
f_al_go <- f_al[c('GO:0019752', 'GO:0006633', 'GO:0044272', 'GO:0042445', 'GO:0030198', 'GO:0009987'), ]
f_dr_go <- f_dr[c('GO:0036503', 'GO:0006486', 'GO:0006457', 'GO:0006986', 'GO:0055085', 'GO:1903506'), ]
# note: the last GO term of each list is used for calibrating the 4 graphs so that their scales are in a more similar range

# order the dataframe based on decreasing padj values
m_al_go <- m_al_go[order(m_al_go$padj), ]
m_dr_go <- m_dr_go[order(m_dr_go$padj), ]
f_al_go <- f_al_go[order(f_al_go$padj), ]
f_dr_go <- f_dr_go[order(f_dr_go$padj), ]

# Change the level of 'Term' based on the order I defined here.
m_al_go$Term <- factor(m_al_go$Term, levels=rev(unique(m_al_go$Term)))
m_al_go['Condition']<- 'Male DR downregulated'

m_dr_go$Term <- factor(m_dr_go$Term, levels=rev(unique(m_dr_go$Term)))
m_dr_go['Condition']<- 'Male DR upregulated'

f_al_go$Term <- factor(f_al_go$Term, levels=rev(unique(f_al_go$Term)))
f_al_go['Condition']<- 'Female DR downregulated'

f_dr_go$Term <- factor(f_dr_go$Term, levels=rev(unique(f_dr_go$Term)))
f_dr_go['Condition']<- 'Female DR upregulated'

# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------

pdf('Plots/Liver_DietDEG_mDown_GO_bubble_220723.pdf', width = 8, height = 10)
ggplot(data=m_al_go, aes(x=Condition, y=Term)) + 
  geom_point(aes(color=enrichment, size = -log10(padj))) + 
  scale_colour_gradient2(high = "dark blue", space = "Lab", na.value = "grey50") +
  labs(title = 'GO terms for diet-DEGs downregulated in male livers', x='Conditions', y = '',size='-log(FDR)', color='Enrichment') 
dev.off() 

pdf('Plots/Liver_DietDEG_mUp_GO_bubble_220723.pdf', width = 8, height = 10)
ggplot(data=m_dr_go, aes(x=Condition, y=Term)) + 
  geom_point(aes(color=enrichment, size = -log10(padj))) + 
  scale_colour_gradient2(high = "red", space = "Lab", na.value = "grey50") +
  labs(title = 'GO terms for diet-DEGs upregulated in male livers', x='Conditions', y = '',size='-log(FDR)', color='Enrichment') 
dev.off() 

pdf('Plots/Liver_DietDEG_fDown_GO_bubble_220723.pdf', width = 8, height = 10)
ggplot(data=f_al_go, aes(x=Condition, y=Term)) + 
  geom_point(aes(color=enrichment, size = -log10(padj))) + 
  scale_colour_gradient2(high = "dark blue", space = "Lab", na.value = "grey50") +
  labs(title = 'GO terms for diet-DEGs downregulated in female livers', x='Conditions', y = '',size='-log(FDR)', color='Enrichment') 
dev.off() 

pdf('Plots/Liver_DietDEG_fUp_GO_bubble_220723.pdf', width = 8, height = 10)
ggplot(data=f_dr_go, aes(x=Condition, y=Term)) + 
  geom_point(aes(color=enrichment, size = -log10(padj))) + 
  scale_colour_gradient2(high = "red", space = "Lab", na.value = "grey50") +
  labs(title = 'GO terms for diet-DEGs upregulated in female livers', x='Conditions', y = '',size='-log(FDR)', color='Enrichment') 
dev.off() 


sessionInfo() 
