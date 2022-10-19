# Title:  Get a large count table to be used as input for DEseq2 analysis
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
# Set wd to the current directory
setwd("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/GetCounts")

# Please create the following folder if not present: 'Output'

# ---------------  Input variables ---------------------------
genename = "gene_name"
homeDir = as.character("/Users/jingxun/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/GetCounts")
sampleTable <- read.csv(paste0(homeDir, "/Input/DRexp20211219ExperimentDesign.csv"), row.names = 1, stringsAsFactors = F)


# ------------------------------------------------------------------
# Define functions for rpkm and tpm
# ------------------------------------------------------------------
## reference: https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## reference: https://www.biostars.org/p/171766/

rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}


                    
# ------------------------------------------------------------------
# Loop to get the counts
# ------------------------------------------------------------------
for (i in 1:length(sampleTable$files)){
  
  outfilename <- file.path(homeDir, "/Output/", paste0(sampleTable$lib[i], "_counts.csv")) # Change this to the path of output file
  print (sampleTable$files[i])
  outfilename
  cts <- read.table(paste0(homeDir, "/featureCount/", sampleTable$files[i], ".featureCounts"), header = TRUE)
  colnames(cts)[7] <- paste("Counts")
  colnames(cts)[1] <- paste("GeneID")
  x = cts[c("GeneID", "Length", "Counts")]
  
  # get tpm, this is not normalized across samples
  x$TPM = tpm(x$Counts, x$Length)  
  write.table(x,file=outfilename,quote=FALSE,sep=",",row.names=FALSE)

  rm(outfilename, cts, x) # remove the temp variables for the next run

}

sessionInfo() 

