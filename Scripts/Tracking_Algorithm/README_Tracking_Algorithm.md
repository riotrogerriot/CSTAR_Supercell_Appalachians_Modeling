# README_Tracking_Algorithm

# Summary

This folder contains a series of interconnected python3 scripts that work together to perform analysis on each individual CM1 model output file within a single simulation. The required packages were managed using Anaconda with the CSTAR_Supercell_Appalachians_Modeling.yml file being included to replicate the environment used to run these scripts at the time of publication. An algorithm flowchart has been included in this folder for visual reference as to how these scripts work together as well.

run_cm1_analysis_v2.py is the main driver code that pulls the individual scripts together to create visualizations and a summary.csv file used to analyze each simulation. This "should" be the only script where the user would need to make any modifications to use these analysis scripts for replication/additional applications.

Note that run_cm1_analysis_mesoanalysis_plotter.py does have a Cython dependancy which is included as srh_cy.pyx and setup_srh.py (provided by Dr. Branden Katona). Users will need to run the following command first to compile the associated C code to allow this module to work as intended: python setup_srh.py build_ext --inplace

run_cm1_analysis_UH_tracker.py is the main driver for the storm-tracking algorithm used in this study. We developed a 2-5 km Updraft Helicity-centric (UH) algorithm that locks-on to a specific updraft after a set model integration time and checks for spatial overlap to ensure consistent tracking over time. Convective mode is also determined by an area-weighted UH average, with 150 m^2 s^-2 being used as the threshold for supercells in our 250 m simulations. Note, this algorithm has not been tested on model ouput with purposes outside the scope of this study so caution is advised if adpoting it for other purposes.

The additional python scripts individually build the associated plots they are named after (i.e., radar_plotter, sounding_plotter, etc.). Again these scripts should not have to be modified as the main driver code was written to allow for flexibilities regarding these plotting functions. 

Feel free to e-mail the corresponding author with any questions/concerns regarding these scripts.

# Folder Contents

1). Storm_Tracking_Algorithm_Flowchart.jpeg:
A visual aid to help understand how these interconnected scripts work together to generate our visual and numeric model output analyses.

2). run_cm1_analysis_UH_Tracker.py
A python3 script containing a user-defined function that controls the consistent tracking and evaluation of a particular updraft after a set model integration time. Built-in areal and temporal overlap checks are included to minimize tracking errors.

3). run_cm1_analysis_cross-section_plotter.py
A python3 script containing a user-defined function that controls the generation and contents of an axis object to visualize a zonal cross-section along the center of the maximum UH value for a respective model output file. Several variables of interest (e.g., radar reflectivity > 30 dBz, relative humidity, parameterized updraft rotation, parameterized cloud contours, 0-3 km cold pool boundaries, and storm-relative winds are included in the resulting plot.

# References

Riggin IV, R. R., C. E. Davenport, M. D. Eastin, K. E. McKeown, S. M. Purpura, and B. T. Katona, 202X: Idealized Simulations of Supercell Thunderstorms Interacting with the Appalachian Mountains. In Prep.
