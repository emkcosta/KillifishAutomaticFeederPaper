Script to do PCA and DEseq2 analysis
=============================================================================================

This script performs PCA and DEseq2 analysis, output differentially expressed genes, and generate MA plots and PCA plots. 

Related to Figure 5B, as well as Figure 5-figure supplements
 

Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
1. The Count supermatrix .csv file generated using the '2_getCombinedCountMatrix_220805_codechk.pl' script
2. Experimental design .csv file that contains library name, sex, feeding, batch, and condition (a complex variable with all the info). 

------------------------------+
Plots   |
------------------------------+
This is where the PCA plots and MA plots will be stored.


-----------------+
Output          |
-----------------+
This is where all the DEseq2 comparisons and normalized counts will be stored. 
