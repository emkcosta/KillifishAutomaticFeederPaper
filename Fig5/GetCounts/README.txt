Script to get a gene count table for each sample based on the feature counts:
=============================================================================================

Feature counts were generated using the Subread package (not allowing multiple mapping). Non-fractional counts are required for DEseq2 to run properly. Be sure to check the output of this script to ensure that the result does not contain fractions. 
 

Directory Structure:
==========================

------------------------------+
Input   |
------------------------------+
An experimental design file that contains a 'file' column that matches the library name of Subread output (this is usually the library name provided by Novogene). 

------------------------------+
featureCount   |
------------------------------+
This is where all the feature count output files from Subread are stored.


-----------------+
Output          |
-----------------+
Output directory will have the output lists in .csv format for each sample. To get all the counts combined by gene name for all samples, use these files as input for the perl script called '2_getCombinedCountMatrix_220805_codechk.pl' and run the perl script on Terminal. 
