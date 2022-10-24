Script for Gene Ontology enrichment
==========================

This script can be run with Biological Process (BP), Molecular Function (MF) and Cellular Component (CC). Change the variable in the script and rerun. See the comments inside the script for more details.


Directory Structure:
==========================

-----------------+
Input       |
-----------------+
The input folder contains the output files of the DEseq2 analysis - (1) comparing AL vs DR conditions, for males and females separately (dietDEGs); and (2) comparing males vs females, for AL and DR conditions separately (sexDEGs).

The script contains codes to separate the different groups of genes (e.g., dietDEGs unregulated in DR females). See the comments inside the script for selecting the specific group of genes to run the script with.  


-----------------+
GO_terms         |     |
-----------------+
These have the latest list of GO and KEGG terms from Ensembl release 100. 


-----------------+
Output         |
-----------------+
Results directory will have the output lists. 
Two per run: one for enriched terms and one for associated gene lists.
Gene list and GO terms can be combined to get a single file with gene ids.
