# CSTAR Supercell Appalachians Modeling Project README
# Corresponding Author: Roger R. Riggin IV, M.S. (rriggin@uncc.edu)
# Latest Update: 07/26/2023

Hello fellow scientists, reviewers, and/or weather-nerds :)

The following repository contains the necessary data, scripts, and compiled executables to review and replicate the results of Riggin IV et al. (202X), which focused on isolated supercell interactions with the complex terrain of the Appalachian Mountains via an idealized numerical model. This study is the third part of a multi-scale analysis of numerous isolated supercells interacting with the complex terrain of the Appalachians (see Purpura et al. 2023 and McKeown et. al 202x for environmental and radar-based analyses). This README will outline the layout of the repository and summarize what you can find in each individual folder. Please feel free to reach out if you have any questions or concerns while utilizing anything within this repository.

Preliminary results are already published in my Master's Thesis which is available via Proquest (Riggin 2021).

The publication (i.e., Riggin IV et al. 202X) is in preparation for submission to Monthly Weather Review. The repository will be updated with appropriate links and references once the peer-review process is complete. A non peer-reviewed draft is available upon request.

# Folder 1): BSS_cm1r20.3_and_CI_Delay_SRC_Files

This folder contains all necessary files that are required to incorporate both Base-State Substitution (BSS; Letkewicz et al. 2013; Davenport et al. 2019) and Convective Initiation (CI) delay functionality into CM1. The BSS code was developed and graciously provided by Dr. Casey Davenport (UNCC; casey.davenport@charlotte.edu). See README within this folder for additional details. 

Note, this functionality is only confirmed to work with cm1r20.3 (June 2021) and has not been tested with the latest update (cm1r21.0; April 2022). 

Included fortran files must be placed in the /cm1r20.3/src/ directory. Additionally, the Makefile will need to be modified for your specific HPC environment.

CM1 is an open-source idealized numerical model maintained by Dr. George Bryan (Bryan and Fritsch 2002). A clean build of cm1r20.3 can be found at https://www2.mmm.ucar.edu/people/bryan/cm1/

# Folder 2): Conceptual_Models

This folder contains full-resolution images of the Conceptual Models that describe four key terrain-induced processes that acted to modulate supercells in our simulations (i.e., Figs. 16-19; Riggin IV et al. 202X). The primary goal of these models is to improve situational awareness and alleviate short-term challenges associated with supercell events around regions of complex terrain. See the publication for a complete walkthrough of each model. A brief summary is included below.

The conceptual models summarize the following terrain-induced processes, how they perturb the near-storm environment, and expected short-term impacts on an isolated supercell passing through an area under such a regime:

1) Terrain-Blocking (Fig. 16)
2) Terrain-Channeling (Fig. 17)
3) Upslope Flow (Fig. 18)
4) Downslope Flow (Fig. 19)

These figures were manually drawn using Adobe Illustrator.


# Folder 3): Scripts

This folder contains the anaconda environment file (i.e, CSTAR_Supercell_Appalachians_Modeling.yml) and python3 scripts used to track and analyze individual supercells in model output (i.e., Tracking_Algorithm), to generate a CM1-ready realistic terrain field (i.e, Generate_RLTRN_Field; graciously provided by Dr. Branden Katona (CIWRO; branden.katona@noaa.gov), and construct the figures used for our analyses in the publication (i.e., Figure_Generation). See individual READMEs for more details on each.

Necessary data to re-run most of these scripts are included (e.g., Cleaned 5-min Output CSVs; Initial Conditions Output). Please contact the corresponding author for any data not provided (unable to upload due to file sizes). 

Additionally, full-resolution figures used in the publication (including annotations and cleaning added post-script) were included for reference.

Lastly, we included the model-executables, necessary input files, analysis scripts, and 5-min output imagery from our Dry_Simulations (i.e., no moisture or microphysics), which were discussed in the publication but not shown for brevity. These simulations utilized fixed-layer Froude Numbers computed at various locations around the modeled terrain fields to confirm the presence of terrain-blocking influences within the included background environments (i.e., upstream crossing and noncrossing composites; Fig. 3a, d; Riggin IV et al. 202X). Our approach was rather simple and served only to confirm blocking potential within the environment. Peak and downstream composites (e.g., Fig. 3b, c, e, f; Riggin IV et al. 202X) were not tested due to inferred terrain-induced perturbations responsible for realizing such spatiotemporal heterogeneities (see Purpura et al. 2023).

# Folder 4): Steady-State

This folder contains all necessary files (e.g., CM1 executables, input_soundings, and namelists) required to re-run each steady-state (Base-State stems from Upstream Composites only; Fig. 3a,d; Riggin IV et al. 202X) simulation. Each simulation folder also includes a onefile.F for reference and 5_Min_Output folder which contains visual references generated by our tracking algorithm (/Scripts/Tracking_Algorithm) that were discussed in the Publication but again not shown for brevity. 

For the Idealized Terrain simulations, the terrain configuration (e.g. TRN and MOD) are built into the executables (see Methods Section in Publication for configuration details) and controlled by namelist settings (i.e., itern = 1).

Realistic Terrain simulations also include a perts.dat file, which is the cm1-ready terrain file generated from USGS DEM of the south-central Appalachians (See Methods Section Riggin IV et al. 202X or /Scripts/Generate_RLTRN_Field/ for more details).

All steady-state executables were built with the same dependencies (e.g., openmpi/4.0.3, hdf5/1.10.5-mpi, and netcdf/4.7.2-mpi) and should be able to be queued in an HPC environment using the submit_steady_state.slurm script if using Slurm as a job manager. Users utilizing other job managers will need to develop their own submit scripts. We suggest a minimum wall time of 72:00 hrs to allow simulations to complete a full integration based on available resources used originally via the UNCC HPC Clusters (e.g., 8 Nodes, 64gb Memory, 32 tasks per node, MPI, etc).

# Folder 5): Variable-State

This folder contains all necessary files (e.g., CM1 executables, input_soundings, and namelists) required to re-run each variable-state (Base-State nudging via BSS from Upstream to Peak and then Downstream Composites; Fig. 3a-f; Riggin IV et al. 202X) simulation. Each simulation folder also includes a onefile.F for reference and 5_Min_Output folder which contains visual references generated by our tracking algorithm (/Scripts/Tracking_Algorithm) that were discussed in the Publication but again not shown for brevity. 

For the Idealized Terrain simulations, the terrain configurations (e.g., BSS_TRN) are built into the executable (see Methods Section in Publication for configuration details) and controlled by namelist settings (i.e., itern = 1).

All variable-state idealized terrain executables (e.g., BSS_IDTRN) were built with the same dependencies (e.g., intel/2020, intel-rtl/2020, openmpi/4.1.0-intel, hdf5/1.10.7-intel-mpi, and netcdf/4.8.1-intel-mpi) and should be able to be queued in an HPC environment using the submit_BSS_IDTRN.slurm script if using Slurm as a job manager. Users utilizing other job managers will need to develop their own submit scripts. Again, we suggest a minimum wall time of 72:00 hrs to allow simulations to complete a full integration based on available resources used original via the UNCC HPC Clusters (e.g., 8 Nodes, 64gb Memory, 32 tasks per node, MPI, etc).

Realistic Terrain simulations also include a perts.dat file, which is the cm1-ready terrain file generated from USGS DEM of the south-central Appalachians (See Methods Section of Riggin IV et al. 202X or /Scripts/Generate_RLTRN_Field/ for more details).

Each variable-state realistic terrain executable (e.g., BSS_CS_RLTRN_ALL and BSS_NC_RLTRN_ALL) were built with the different dependencies and should be able to be queued in an HPC environment using unique submit scripts in each respective simulations folder when using Slurm as a job manager. Users utilizing other job managers will need to develop their own submit scripts. We suggest a minimum wall time of 96:00 hrs to allow simulations to complete a full integration based on available resources used originally via the UNCC HPC Clusters (e.g., 8 Nodes, 64gb Memory, 32 tasks per node, MPI, etc).

Feel free to contact the corresponding author regarding re-runs and dependencies.

# References

Bryan, G. H., and J. M. Fritsch, 2002: A benchmark simulation for moist nonhydrostatic numerical models. Mon. Wea. Rev., 130, 2917–2928, https://doi.org/10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2.

Davenport, C. E., C. L. Ziegler, and M. I. Biggerstaff, 2019: Creating a More Realistic Idealized Supercell Thunderstorm Evolution via Incorporation of Base-State Environmental Variability. Mon. Wea. Rev., 147, 4177–4198, https://doi.org/10.1175/MWR-D-18-0447.1.

Letkewicz, C. E., A. J. French, and M. D. Parker, 2013: Base-state substitution: An idealized modeling technique for approximating environmental variability. Mon. Wea. Rev., 141, 3062–3086, https://doi.org/10.1175/MWR-D-12-00200.1.

McKeown, K. E., C. E. Davenport, M. D. Eastin, S. M. Purpura, and R. R. Riggin IV, 202X: Radar Characteristics of Supercell Thunderstorms Interacting with the Appalachian Mountains. Wea. Forecasting. Manuscript submitted for review.

Purpura, S. M., C. E. Davenport, M. D. Eastin, K. E. McKeown, and R. R. Riggin, 2023: Environmental Evolution of Supercell Thunderstorms Interacting with the Appalachian Mountains. Wea. Forecasting, 38, 179–198, https://doi.org/10.1175/WAF-D-22-0115.1.

Riggin , Roger R., IV. (2022). Idealized simulations of supercell thunderstorms interacting with the appalachian mountains (Order No. 29397578). Available from Dissertations & Theses @ University of North Carolina Charlotte. (2726338467). https://www.proquest.com/dissertations-theses/idealized-simulations-supercell-thunderstorms/docview/2726338467/se-2

Riggin IV, R. R., C. E. Davenport, M. D. Eastin, K. E. McKeown, S. M. Purpura, and B. T. Katona, 202X: Idealized Simulations of Supercell Thunderstorms Interacting with the Appalachian Mountains. In Prep.
