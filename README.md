# CSTAR_Supercell_Appalachians_Modeling README

 Hello,
  The following repository contains all necessary files and scripts to replicate the numerical modeling and analysis
 work from my thesis and upcoming publication, Riggin et al. (202X). The following document will outline what files 
 are located where and how to utilize them in order to replicate and/or advance the current work.

Folder 1): Scripts
  This folder contains a series of interconnected python3 scripts that work together to perform analysis on each
 individual CM1 model output file within a single simulation. The required packages were managed using Anaconda
 with the CSTAR_Supercell_Appalachians_Modeling.yml file being included to replicate the environment used to run
 these scripts at the time of publication.
 
 run_cm1_analysis_v2.py is the main driver code that pulls the individual scripts together to create visualizations
 and a summary.csv file used to analyze each simulation. This "should" be the only script where the user would need
 to make any modifications to use these analysis scripts for replication/additional applications.
  
 Do note that run_cm1_analysis_mesoanalysis_plotter.py does have a Cython dependancy which is included as srh_cy.pyx 
 and setup_srh.py (provided by Dr. Branden Katona). Users will need to run the following command first to compile the
 associated C code to allow this module to work as intended: python setup_srh.py build_ext --inplace
  
 run_cm1_analysis_UH_tracker.py is the main driver for the storm-tracking algorithm used in this study. We developed
 a UH-centric algorithm that hons in on a specific updraft after a set model integration time and checks for spatial
 overlap to ensure consistent tracking over time. Do note that this algorithm has not been tested on other model ouputs
 outside of the scope of this study so caution is advised if anyone wishes to adopt if for other applications. See the
 paper for additional details.
  
 The additional python scripts individually build the associated plots they are named after (i.e., radar_plotter, 
 sounding_plotter, etc.). Again these scripts should not have to be modified as the main driver code was written
 to allow for flexibilities regarding these plotting functions. Feel free to e-mail me with any questions that may
 arise when working with these scripts.
  
 Folder 2): Steady-State
   This folder contains the necessary files to re-run any of our steady-state (Upstream Only) simulations. Within both
  the Idealized Terrain and Realistic Terrain folders, there are individual folders containing the cm1.exe, onefile.F,
  namelist, input_sounding, and perts.dat (Realistic Terrain only) files needed to re-run the simulations exactly as they 
  were designed during this study. I have also included a sample Slurm submit script that we used to submit these jobs to the 
  UNCC HPC clusters. Note that the slurm scripts contain the module dependencies that the cm1.exe files were compiled with so 
  be sure to check them out even if you are using a different job manager than Slurm.
   
 Folder 3): Variable-State
   This folder contains the necessary files to re-run any of our variable-state (BSS) simulations. Within both
  the Idealized Terrain and Realistic Terrain folders, there are individual folders containing the cm1.exe, onefile.F,
  namelist, input_sounding, and perts.dat (Realistic Terrain only) files needed to re-run the simulations exactly as they 
  were designed during this study. Additionally, utlizing BSS requires additional sounding files (i.e., thermo_replace1, 
  wind_replace1, etc.) I have also included a sample Slurm submit script that we used to submit these jobs to the 
  UNCC HPC clusters. Note that the slurm scripts contain the module dependencies that the cm1.exe files were compiled 
  with so be sure to check them out even if you are using a different job manager than Slurm.
   
Please feel free to reach out if there are any specific questions (rriggin@uncc.edu or chaos93677@me.com).

Riggin , Roger R., IV. (2022). Idealized simulations of supercell thunderstorms interacting with the appalachian mountains (Order No. 29397578). Available from Dissertations & Theses @ University of North Carolina Charlotte. (2726338467). Retrieved from https://www.proquest.com/dissertations-theses/idealized-simulations-supercell-thunderstorms/docview/2726338467/se-2

A link to the publication will be provided once it has been accepted into an AMS Journal. 
