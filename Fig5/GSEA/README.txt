Scripts to do enrichments based on Gene Set Enrichment Analysis (GSEA) using ClusterProfiler:
=============================================================================================

GSEA works slightly better with low numbers because it uses the whole lists into account. 
The input must be a ranked list (ranked based on anything e.g. FDR, Fold Change etc.). I 
rank based on: <-log10(qvalue) * FC>. This takes both FDR and fold change into account.
The input file is the DEseq2 output file that contains p-values and log2FoldChange.

4_GSEA_220805_codechk.R: the script that takes the input list and ranks and runs GSEA
using Gene Ontology

5_GSEA_BubblePlot_220805_codechk.R: the script that plots the GSEA results for males and females together

6_GSEA_Heatmap_220805_codechk.R: the script that plots heatmaps using the scaled normalized counts of the genes in selected GO terms.   

Refer to ClusterProfiler vignette for more details on tests and plots.

Directory Structure:
==========================

-----------------+
Input          |
-----------------+
The input files are the output files of the DEseq2 analysis comparing AL vs DR conditions, for males and females separately. 

-----------------+
Output          |
-----------------+
4_GSEA_220805_codechk.R: output lists in .csv format for all different enrichments. Gene names - human orthologs will be in the last column.

-----------------+
Plots          |
-----------------+
5_GSEA_BubblePlot_220805_codechk.R: Dot plot. Figure 5C
6_GSEA_Heatmap_220805_codechk.R: Heatmaps. Figure 5D



