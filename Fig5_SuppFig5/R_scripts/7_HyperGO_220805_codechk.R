# Title: PCA and DEseq2 analysis comparing across tissues, sexes, and dietary conditions 
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/HyperGO")
set.seed(1234)

## This is to do GO enrichment analysis based on BH list from zebrafish/human
# *** IMPORTANT: Run it one by one for each of BP, MF and CC in the function #####
library("GOstats")
library('dplyr')



# --------------------------------------------- 
# Load data
# ---------------------------------------------
# Pick a dataset to run the script Leave the rest labeled with '#' sign. 

# dietDEG data
# data <- read.csv("Input/liver_Female_DRoverAL_DEG_220401.csv", header = TRUE)
# data <- read.csv("Input/liver_Male_DRoverAL_DEG_220401.csv", header = TRUE)

# sexDEG data
data <- read.csv("Input/liver_SexDiverged_AL_MoverF_DEG_220426.csv", header = TRUE)
# data <- read.csv("Input/liver_SexDiverged_DR_MoverF_DEG_220714.csv", header = TRUE)


# --------------------- Process the data ------------------------

# Find the significant differential-expressed genes, then categorize them based on up-regulation or down-regulation 

# For dietDEG data, run these two lines of code:
# data_dr_UP <- subset(data, padj < 0.05 & log2FoldChange > 0)
# data_al_UP <- subset(data, padj < 0.05 & log2FoldChange < 0)

# For sexDEG data, run these two lines of code:
data_m_UP <- subset(data, padj < 0.05 & log2FoldChange > 0)
data_f_UP <- subset(data, padj < 0.05 & log2FoldChange < 0)


# --------------------------------------------- 
# Define background and genes of interest
# ---------------------------------------------
# List of universe genes (removed padj = NA). Background.
universe <- subset(data, !is.na(data$padj))
universe <- subset(universe, select = 'Gene')
colnames(universe) <- 'id'

# ---------------------- Select Gene list to be tested -----------------------
# Select one of the 4 options, given the specific data set and the condition of interest (up-regulation or down-regulation)  

# For dietDEG data, pick one of these two options:
# genes <- subset(data_dr_UP, select = 'Gene')
# genes <- subset(data_al_UP, select = 'Gene')

# For sexDEG data, pick one of these two options:
# genes <- subset(data_m_UP, select = 'Gene')
genes <- subset(data_f_UP, select = 'Gene')


# clean up dataframe
colnames(genes) <- 'id'


# ---------------------- Select Output files ----------------------
# Select the proper output file name given the specific data set, as well as upregulation or downregulation (condition of interest)  

# For dietDEG data, pick one of these options:
# outfilename = "Output/GOBP_dietDEG_f_DR_UP_220714.txt" #enrichment for the dietDEGs upregulated in DR females
# outfilename = "Output/GOBP_dietDEG_f_AL_UP_220714.txt" #enrichment for the dietDEGs downregulated in DR females
# outfilename = "Output/GOBP_dietDEG_m_DR_UP_220714.txt" #enrichment for the dietDEGs upregulated in DR males
# outfilename = "Output/GOBP_dietDEG_m_AL_UP_220714.txt" #enrichment for the dietDEGs downregulated in DR males

# For sexDEG data, pick one of these options: 
# outfilename = "Output/GOBP_sexDEG_AL_mUP_220705.txt" #enrichment for the sexDEGs identified in AL and higher expression in males
outfilename = "Output/GOBP_sexDEG_AL_fUP_220705.txt" #enrichment for the sexDEGs identified in AL and higher expression in females
# outfilename = "Output/GOBP_sexDEG_DR_mUP_220716.txt" #enrichment for the sexDEGs identified in DR and higher expression in males
# outfilename = "Output/GOBP_sexDEG_DR_fUP_220716.txt" #enrichment for the sexDEGs identified in DR and higher expression in females

# ---------------------- Select Output file with genes in GO terms ----------------------
# Select the proper output file name to record the genes associated with each GO term.

# For dietDEG data, pick one of these options:
# gotermlist = "Output/GeneGOBP_dietDEG_f_DR_UP_220714.txt"
# gotermlist = "Output/GeneGOBP_dietDEG_f_AL_UP_220714.txt"
# gotermlist = "Output/GeneGOBP_dietDEG_m_DR_UP_220714.txt"
# gotermlist = "Output/GeneGOBP_dietDEG_m_AL_UP_220714.txt"

# For sexDEG data, pick one of these options: 
# gotermlist = "Output/GeneGOBP_sexDEG_AL_mUP_220705.txt"
gotermlist = "Output/GeneGOBP_sexDEG_AL_fUP_220705.txt"
# gotermlist = "Output/GeneGOBP_sexDEG_DR_mUP_220716.txt"
# gotermlist = "Output/GeneGOBP_sexDEG_DR_fUP_220716.txt"

# --------------------------------------------- 
# Hypergeometric test for GO analysis
# ---------------------------------------------

# --------------------------------------------- INPUT FILES --------------------------------------------------------
# Go terms list. Use either zebrafish or human. Human works a bit better due to better annotations.
frame = read.table(file ="GO_terms/GO_killifish-human_best-hits.txt", header = T, colClasses=c(rep("factor",3)))
# ontology MF, BP, CC
ontolg = "BP"

# Minimum number of genes for a term to filter
mingenes = 5 # Bare minimum is 2. More will get more general terms, less more specific. 5-10 is a good number.
# Relative enrichment filter
relenrich = 0 # I generally use 0 to get all terms

# ------------------------------------------------------------------------------------------------------------------

# This is just to get the 3 column. I already have these so I don't need it, still.
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)

# put your data into a GOFrame object
goFrame=GOFrame(goframeData,organism="Human")
head(goframeData)

# cast this object to a GOAllFrame object will tap into the GO.db package and populate this object with the implicated GO2All mappings for you
goAllFrame=GOAllFrame(goFrame)

# generate geneSetCollection objects
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# Process the universe list
universe = universe$id
universe = lapply(universe, as.character)
universe = unlist(universe)
head(universe)

# Process the gene list of interest
genes = genes$id
genes = lapply(genes, as.character)
genes = unlist(genes)
head(genes)

params <- GSEAGOHyperGParams(name="My Custom GSEA based annotation parameters", 
                             geneSetCollection=gsc, 
                             geneIds = genes, 
                             universeGeneIds = universe, 
                             ontology = ontolg,
                             pvalueCutoff = 1, # 1 will get all terms, and then we can filter later
                             conditional = F, # To consider GO DAG structure or not. Doesn't affect much, but we can always try
                             testDirection = "over") # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over". Fo depleted terms use under.

# call hyperGTest to do the test
Over <- hyperGTest(params)
head(summary(Over))

# calculate enrichment and add it to data frame.
# Relative enrichment factor (E-value) for a GO term = (count/size)/(size/universe)
enrichment = (summary(Over)[5]$Count / summary(Over)[6]$Size) / (summary(Over)[6]$Size / length(universe))

# create a new frame
SummaryOver = data.frame(summary(Over), enrichment)
head(SummaryOver)

# Filter the Over variable on parameters other than P-value
# Filter the summary of OVER with size of the term, at least 2 genes for a go term
FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$enrichment >= relenrich),]
head(FilteredSummaryOver)

# adjust p value for multiple correction
padj = p.adjust(FilteredSummaryOver$Pvalue, "BH")

# Add padj to the data frame
FinalSummaryOver = data.frame(FilteredSummaryOver, padj)

# write to a file
write.table(FinalSummaryOver, outfilename, quote = F, row.names = F, sep = "\t") # Final result file


# --------------------- To get the genes for each GO terms ---------------------------------------------

# isolate indexes for the go terms in final results
ind.GO <- is.element(names(Over@goDag@nodeData@data), eval(parse(text=paste("FinalSummaryOver$", "GO",ontolg,"ID", sep=''))))
selected.GO <- Over@goDag@nodeData@data[which(ind.GO)]

# get a go terms and genes in a new variable for all the terms in the results of enrichment
goTerms <- lapply(selected.GO, 
                  function(x) x$geneIds)
names(goTerms) <- names(Over@goDag@nodeData@data)[ind.GO]

# This will create a new file that will have GO terms and gene names in each GO terms.
# Number of GO terms or lines should be equal to the enriched go terms as in the other file generated by this script.
# This needs to be further processed to generate the desired files.

for (i in 1:length(goTerms)){
  
  test = as.data.frame(do.call(rbind, goTerms[i]))
  write.table(test, file = gotermlist, quote = F, col.names = F, append = T) # append each line - so make sure the file is empty for each run, or renamed after each run
  rm(test)
}


rm(list = ls())

sessionInfo() 
