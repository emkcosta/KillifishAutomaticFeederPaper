Script to Process Outputs of DeepLabCut
=============================================================================================
1_DLCOutputPreProcessing.R

This script takes the CSV outputs of DeepLabCut, and performs likelihood and euclidean distance filtering, interpolation, and ycomponent of velocity calculations.

Related to Figure 6, as well as Supplemental Figure 6.
 

Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
The CSV file outputs from DLC for each of the animals. These are stored in "Figure6-SourceData2_unfiltered_raw_data_DeepLabCut.zip".



------------------------------+
Output   |
------------------------------+
This is where the files with filtered and interpolated trajectories and velocity calculations will be stored. The files that are generated are included in "Figure6-SourceData3_filtered_interpolated_data.zip". 