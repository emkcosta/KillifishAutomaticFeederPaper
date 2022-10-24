Script to Generate Velocity and Trajectory Plots and Calculate t1
=============================================================================================
2_Kinematics_and_t1_Fig6.R

This script takes in the CSVs with velocity calculations (Output from 1_DLCOutputPreprocessing.R) and generates velocity heat map plots, fish trajectory plots,  and calculations for surface distance and t1. 

Related to Figure 6, as well as Supplemental Figure 6.
 

Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
1. The CSV file outputs with instantaneous velocity data from 1_DLCOutputPreprocessing.R for all animals. The files are included in "Figure6-SourceData3_filtered_interpolated_data.zip". This zip file was too big to upload here.

2. The CSV file output with prefiltered 'allanimals_allREC_resnet50_TrainAllTheThingsJun28shuffle1_1030000.csv' from 1_DLCOutputPreprocessing.R

------------------------------+
Plots   |
------------------------------+
1. Velocity heat maps

2. Fish trajectory plots 

3. Surface distance heat map. This is not a figure in the paper but can be useful for visualization and sanity checking purposes.


------------------------------+
Output   |
------------------------------+
1. CSVs containing surface distance data

2. CSVs containing t1 data
