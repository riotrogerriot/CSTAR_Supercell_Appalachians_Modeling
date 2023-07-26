#---------------------------------------------------------------------------------
# CSTAR Supercell Appalachians Modeling Project README
# Author: Roger R. Riggin IV, M.S. (PhD Student at UNCC)
# Latest Update: 07/26/2023
#---------------------------------------------------------------------------------

Hello fellow scientists, reviewers, and/or weather-nerds :)

The following repository contains the necessary data, scripts, and compiled executables to review and replicate the results of Riggin et al. (202X; In Prep for MWR Submission), which focused on isolated supercell interactions with the complex terrain of the Appalachian Mountains via an idealized numerical model. This README will outline the layout of the repository and summarize what you can find in each individual folder. Please feel free to reach out if you have any questions or concerns while utilizing anything within this repository.

Contact Info:
Academic E-mail: rriggin@uncc.edu
Personal E-mail: chaos93677@me.com

Publication Link: TBA. Feel free to contact me for a non peer-reviewed draft 

Publication Citation
#-------------------
TBA

Thesis Link: https://www.proquest.com/dissertations-theses/idealized-simulations-supercell-thunderstorms/docview/2726338467/se-2

Thesis Citation
#--------------
Riggin , Roger R., IV. (2022). Idealized simulations of supercell thunderstorms   interacting with the appalachian mountains (Order No. 29397578). Available from Dissertations & Theses @ University of North Carolina Charlotte. (2726338467). 

#---------------------------------------------------------------------------------
Folder 1): BSS_cm1r20.3_and_CI_Delay_SRC_Files
#---------------------------------------------------------------------------------

This folder contains all necessary files that are required to incorporate both Base-State Substitution (BSS; Letkewicz et al. 2013; Davenport et al. 2019) and Convective Initiation (CI) delay functionality into CM1. The BSS code was developed and graciously provided by Dr. Casey Davenport (UNCC; casey.davenport@charlotte.edu). See README within folder for additional details. 

Note, this functionality is only confirmed to work with cm1r20.3 (June 2021) and has not been tested with the latest update (cm1r21.0; April 2022). 

Included fortran files must be placed in the /cm1r20.3/src/ directory. Additionally, the Makefile will need to be modified for your specific HPC environment.

CM1 is an open-source idealized numerical model maintained by Dr. George Bryan (Bryan and Fritsch 2002). A clean build of cm1r20.3 can be found at https://www2.mmm.ucar.edu/people/bryan/cm1/

#---------------------------------------------------------------------------------
Folder 2): Conceptual_Models
#---------------------------------------------------------------------------------

This folder contains full-resolution images of the Conceptual Models that describe four key terrain-induced processes that acted to modulate supercells in our simulations (i.e., Figs. 16-19). The primary goal of these models is to improve situational awarness and aleviate short-term challenges associated with supercell events around regions of complex terrain. See the publication for a complete walkthrough of each model. A brief summary is included below.

The conceptual models summarize the following terrain-induced processes, how they peturb the near-storm environment, and expected short-term impacts on an isolated supercell passing through an area under such a regime:

1) Terrain-Blocking (Fig. 16)
2) Terrain-Channeling (Fig. 17)
3) Upslope Flow (Fig. 18)
4) Downslope Flow (Fig. 19)

These figures were manually drawn using Adobe Illustrator.

#---------------------------------------------------------------------------------
Folder 3): Scripts
#---------------------------------------------------------------------------------

This folder contains the anaconda environment file (i.e, CSTAR_Supercell_Appalachians_Modeling.yml) and python3 scripts used to track and analyze individual supercells in model output (i.e., Tracking_Algorithm), t generate a CM1-ready realistic terrain field (i.e, Generate_RLTRN_Field; graciously provided by Dr. Branden Katona (CIWRO; branden.katona@noaa.gov), and construct the figures used for our analyses in the publication (i.e., Figure_Generation). See indidivudal READMEs for more details on each.

Necessary data to re-run most of these script are included (e.g., Cleaned 5-min Output CSVs; Initial Conditions Output). Please contact corresponding author for any data not provided (unable to upload due to file sizes). 

Additionally, full-resolution figures used in the publication (including annotations and cleaning added post-script) were included for reference.

Lastly, we included the model-executables, necessary input files, analysis scripts, and 5-min output imagery from our Dry_Simulations (i.e., no moisture or microphysics), which were discussed in the publication but not shown for brevity. These simulations utilized fixed-layer Froude Numbers computed at various locations around the modeled terrain fields to confirm the presence of terrain-blocking influences within the included background environments (i.e., upstream crossing and noncrossing composites; Fig. 3a, d). Our approach was rather simple and served only to confirm blocking potential within the environment. Peak and downstream composites (e.g., Fig. 3b, c, e, f) were not tested due to inferred terrain-induced perturbations responsible for realizing such spatiotemporal heterogenities (see Purpura et al. 2023 for more details).

#---------------------------------------------------------------------------------
Folder 4): Steady-State
#---------------------------------------------------------------------------------

This folder contains all necessary files (e.g., CM1 executables, input_soundings, and namelists) required to re-run each steady-state (Base-State stems from Upstream Composites only; Fig. 3a,d) simulation. Each simulation folder also includes a onefile.F for reference and 5_Min_Output folder which contains visual references generated by our tracking algoritm (/Scripts/Tracking_Algoritm) that were discussed in the Publication but again not shown for brevity. 

For the Idealized Terrain simulations, the terrain configuration (e.g. TRN and MOD) are built into the executables (see Methods Section in Publication for configuration details) and controlled by namelist settings (i.e., itern = 1).

Realistic Terrain simulations also include a perts.dat file, which is the cm1-ready terrain file generated from USGS DEM of the south-central Appalachians (See Methods Section in Publication or /Scripts/Generate_RLTRN_Field/ for more details).

All steady-state executables were built with the same dependencies (e.g., openmpi/4.0.3, hdf5/1.10.5-mpi, and netcdf/4.7.2-mpi) and should be able to be queued in an HPC environment using the submit_steady_state.slurm script if using Slurm as a job manager. Users utilizing other job managers will need to develop their own submit scripts. We suggest a minimum wall time of 72:00 hrs to allow simulations to complete a full integration based on available resources used original via the UNCC HPC Clusters (e.g., 8 Nodes, 64gb Memory, 32 tasks per node, MPI, etc).

#---------------------------------------------------------------------------------
Folder 5): Variable-State
#---------------------------------------------------------------------------------

This folder contains all necessary files (e.g., CM1 executables, input_soundings, and namelists) required to re-run each variable-state (Base-State nudging via BSS from Upstream to Peak and then Downstream Composites; Fig. 3a-f) simulation. Each simulation folder also includes a onefile.F for reference and 5_Min_Output folder which contains visual references generated by our tracking algoritm (/Scripts/Tracking_Algoritm) that were discussed in the Publication but again not shown for brevity. 

For the Idealized Terrain simulations, the terrain configurations (e.g., BSS_TRN) are built into the executable (see Methods Section in Publication for configuration details) and controlled by namelist settings (i.e., itern = 1).

All variable-state idealized terrain executables (e.g., BSS_IDTRN) were built with the same dependencies (e.g., intel/2020, intel-rtl/2020, openmpi/4.1.0-intel, hdf5/1.10.7-intel-mpi, and netcdf/4.8.1-intel-mpi) and should be able to be queued in an HPC environment using the submit_BSS_IDTRN.slurm script if using Slurm as a job manager. Users utilizing other job managers will need to develop their own submit scripts. Again, we suggest a minimum wall time of 72:00 hrs to allow simulations to complete a full integration based on available resources used original via the UNCC HPC Clusters (e.g., 8 Nodes, 64gb Memory, 32 tasks per node, MPI, etc).

Realistic Terrain simulations also include a perts.dat file, which is the cm1-ready terrain file generated from USGS DEM of the south-central Appalachians (See Methods Section in Publication or /Scripts/Generate_RLTRN_Field/ for more details).

Each variable-state realistic terrain executable (e.g., BSS_CS_RLTRN_ALL and BSS_NC_RLTRN_ALL) were built with the different dependencies and should be able to be queued in an HPC environment using unique submit scripts in each respecitve simulations folder when using Slurm as a job manager. Users utilizing other job managers will need to develop their own submit scripts. We suggest a minimum wall time of 96:00 hrs to allow simulations to complete a full integration based on available resources used original via the UNCC HPC Clusters (e.g., 8 Nodes, 64gb Memory, 32 tasks per node, MPI, etc).

Feel free to contact the corresponding author regarding re-runs and dependencies.
