# README_Tracking_Algorithm

# Summary

This folder contains a series of interconnected python3 scripts (Van Rossum 2009) that work together to perform analysis on each individual CM1 model output file stemming from a single simulation. The required package dependencies were managed using Anaconda with the CSTAR_Supercell_Appalachians_Modeling.yml file being included for replication purposes. Additionally, an algorithm flowchart is included in this folder for visual reference as to how these scripts work together as well.

run_cm1_analysis_v2.py is the main driver code that pulls the individual scripts together to create visualizations and a summary.csv file used to analyze each simulation. This "should" be the only script where the user would need to make any modifications to use these analysis scripts for replication/additional applications.

run_cm1_analysis_UH_tracker.py is the primary script responsible for storm-tracking. We developed a 2-5 km Updraft Helicity-centric (UH) algorithm that locks-on to a specific updraft after a set model integration time and checks for spatial overlap to ensure consistent tracking over time. Convective mode is also determined by an area-weighted UH average, with 150 m^2 s^-2 being used as the threshold for supercells in our 250 m simulations. Note, this algorithm has not been tested on model ouput with purposes outside the scope of this study so caution is advised if adpoting it for other purposes.

The additional python scripts individually build the associated plots they are named after (i.e., radar_plotter, sounding_plotter, etc.). Again these scripts should not have to be modified as the main driver code was written to allow for flexibilities regarding these analysis and plotting functions. 

These scripts require CF-Compliance NetCDF formatted ouput and are highly-dependent on the MetPy (May et al. 2022), Matplotlib (Hunter 2007), and Xarray (Hoyer and Hamman 2017) python packages.

Feel free to e-mail the corresponding author with any questions/concerns regarding these scripts.

# Folder Contents

1). Storm_Tracking_Algorithm_Flowchart.jpeg:
A visual aid to help understand how these interconnected scripts work together to generate our visual and numeric model output analyses.

2). run_cm1_analysis_UH_Tracker.py: 
A python3 script containing a user-defined function that controls the consistent tracking and evaluation of a particular updraft after a set model integration time. Built-in areal and temporal overlap checks are included to minimize tracking errors. Mesocyclone location and areal-averaged characteristics and calculated and passed back to the main driver code to be written out in the Model_Output.CSV. Convective mode is also determined through a boolean comparison based on the areal-averaged mesocyclone intensity (i.e., Areal-Averaged 2-5 km UH) and a user-defined threshold for supercellular convection.

3). run_cm1_analysis_cross-section_plotter.py: 
A python3 script containing a user-defined function that controls the generation and contents of an axis object to visualize a zonal cross-section along the center of the maximum UH value for a respective model output file. Several variables of interest (e.g., radar reflectivity > 30 dBz, relative humidity, parameterized updraft rotation, parameterized cloud contours, 0-3 km cold pool boundaries, and storm-relative winds) are included in the resulting plot.

4). run_cm1_analysis_get_meso_area.py: 
A python3 script containing a user-defined function that evaluates mesocyclone areas at various levels and surface cold pool area within a defined box centered on the local-UH max (i.e., within the overlap checks from UH_Tracker) for any respective model output file. Mesocyclone depth is also estimated by computing the difference between the highest and lowest gridpoint satisfying a set vertical vorticity value within the grid column corresponding to the local-UH max location. Cold pool intensity is also estimated by the minimum theta_pert value within the user-defined box. Resultant calculations are passed back to the main driver code and written out to the Model_Output.csv. 

5). run_cm1_analysis_initial_conditions.py: 
A python3 script containing a user-defined function to compute the initial 0-6 km mean wind and Bunker's estimated storm motion vectors within a given model grid column. Needed to compute

6). run_cm1_analysis_mesoanalysis_plotter.py: 
A python3 script containing a user-defined function that controls the generation and contents of an axis object to visualize a mesoscale environmental analysis of the near-storm environment around a simulated supercell. Several variables of interest (e.g., 1 km radar reflectivity > 10 dBz, terrain, CAPE > 500 J/kg, CIN > 50 J/kg, UH Swaths, and multi-level wind vectors) are included in the resulting plot.

NOTE: This script does have a Cython (i.e., C integrated with Python) dependency to compute 0-1 km SRH across the entire model domain; included as srh_cy.pyx and setup_srh.py (provided by Dr. Branden Katona; branden.katona@noaa.gov). Users will need to run the following command first to compile the associated C code to allow this module to work as intended: python setup_srh.py build_ext --inplace. 

7). run_cm1_analysis_radar_plotter.py: 
A python3 script containing a user-defined function that controls the generation and contents of an axis object to visualize a simulated supercell. Several variables of interest (e.g., 1 km radar reflectivity > 10 dBz, surface wind vectors, infow sounding location, UH max location, cold pool boundaries, and the user-defined window in which meso_area parameters were computed within) are included in the resulting plot.

8). run_cm1_analysis_sounding_analysis.py: 
A python3 script containing a user-defined function that computes all relevant inflow sounding parameters and passes them back to the main driver code to be written out to the Model_Output.csv.

9). run_cm1_analysis_sounding_plotter.py: 
A python3 script containing a user-defined function that controls the generation and contents of an axis object to visualize the near-inflow sounding along with annotations of relevant parameters computed in the sounding_analysis script.

10). run_cml_analysis_v2.py: 
A python3 script that acts as the main driver code to link the aformentioned scripts together to perform the associated analyses for each output file within a single simulations. The script is written so that all user-defined parameters can be adjusted within this script. The provided flowchart provides a visual aid for the fundamental layout of this code. This script also generates and writes out visual and numeric data for each output file after all analyses are complete.

11). setup_srh.py: 
A python 3 script that is used to build and compile srh_cy.c (Cython code to compute SRH throughout model domain).

12). srh_cy.c:
A copy of the compiled C code after using the "python setup_srh.py build_ext --inplace" command as required to allow the mesoanalysis plotter to function as intended.

13). srh_cy.pyx:
The Cython script developed by Dr. Branden Katona to compute SRH across the model domain.

# References

Hoyer, S. & Hamman, J., (2017). xarray: N-D labeled Arrays and Datasets in Python. Journal of Open Research Software. 5(1), p.10. DOI: https://doi.org/10.5334/jors.148

Hunter J. D., (2007). "Matplotlib: A 2D Graphics Environment," in Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, doi: 10.1109/MCSE.2007.55.

May, R. M., Goebbert, K. H., Thielen, J. E., Leeman, J. R., Camron, M. D., Bruick, Z., Bruning, E. C., Manser, R. P., Arms, S. C., and Marsh, P. T., 2022: MetPy: A Meteorological Python Library for Data Analysis and Visualization. Bull. Amer. Meteor. Soc., 103, E2273-E2284, https://doi.org/10.1175/BAMS-D-21-0125.1.

Riggin IV, R. R., C. E. Davenport, M. D. Eastin, K. E. McKeown, S. M. Purpura, and B. T. Katona, 202X: Idealized Simulations of Supercell Thunderstorms Interacting with the Appalachian Mountains. In Prep.

Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.
