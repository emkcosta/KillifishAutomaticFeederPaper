Scripts to Generate Compass Plots & to Calculate Average Velocity Across all Trials
=============================================================================================
3_CompassPlot_SuppFig6.R

This script takes the CSV files with prefiltered and interpolated positional data from 1_DLCOutputPreprocessing.R and generates compass plots for the first 7 and last 7 trials for the animal of interest.

Related to Supplemental Figure 6.
 

Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
The CSV file outputs with prefiltered and interpolated data from 1_DLCOutputPreprocessing.R.


------------------------------+
Plots   |
------------------------------+
Compass plots for the first and last 7 trials will be stored here.


------------------------------+
Output   |
------------------------------+
1. Aggregate CSV of filtered and interpolated fish trajectories will be stored here.

2. CSVs containing information about the angles used for input into the compass plot will also be stored here.

====================================================
4_AvgVelocity_SuppFig6.R

This script takes the CSV files with velocity calculations output from 2_Kinematics_and_t1_Fig6.R and calculates the average 20-frame rolling average velocity for each animal across all trials.

Related to Supplemental Figure 6.


Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
The CSV file outputs with velocity data from 2_Kinematics_and_t1_Fig6.R.


------------------------------+
Plots   |
------------------------------+
NA

------------------------------+
Output   |
------------------------------+
Aggregate CSV of average fish velocities will be stored here.




 
