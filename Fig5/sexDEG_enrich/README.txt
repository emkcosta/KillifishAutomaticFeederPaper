Script for calculating sexDEG enrichment in dietDEGs
==========================

This script calculates the observed enrichments of sexDEGs (identified in AL or DR conditions) in the dietDEGs found in males or in females. Fisher's exact test was used to test for statistical significance of the enrichment. 

It also generates non-dietDEG control gene sets for male dietDEGs or female dietDEGs, respectively. These non-dietDEG control gene sets share the same number of genes as the dietDEGs, and each control gene is randomly selected from the control gene group that has an expression value within 2% of a given gene in the dietDEG list. The enrichment of sexDEGs in the non-dietDEG control group was simulated 1000 times. Fisher's exact test was run for each bootstrap iteration. 

Lastly, the script plots Pi-charts to represent the enrichments.


Directory Structure:
==========================

-----------------+
Input       |
-----------------+
The input folder contains the output files of the DEseq2 analysis - (1) comparing AL vs DR conditions, for males and females separately (dietDEGs); (2) comparing males vs females, for AL and DR conditions separately (sexDEGs); and (3) normalized counts for all liver samples. 


-----------------+
Output         |
-----------------+
The simulated data will be stored here. 


-----------------+
Plots         |
-----------------+
The Pi-charts will be stored here. 
