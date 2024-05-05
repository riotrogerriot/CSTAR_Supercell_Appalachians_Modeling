# README_Dry_Simulations

README file discussing contents of /Scripts/Dry_Simulations folder
Author: R.R. Riggin, July 2023 (Updated 5/5/2024)

# Summary 

We've included the model-executables, necessary input files, analysis scripts, and 5-min output imagery from our Dry_Simulations (i.e., no moisture or microphysics), which were discussed in the publication but not shown for brevity. These simulations utilized fixed-layer Froude Numbers computed at various locations around the modeled terrain fields to confirm the presence of terrain-blocking influences within the included background environments (i.e., upstream crossing and noncrossing composites; Fig. 3a, d; Riggin et al. In Review). Our approach was rather simple and served only to confirm blocking potential within the environment. Peak and downstream composites (e.g., Fig. 3b, c, e, f; Riggin et al. In Review) were not tested due to inferred terrain-induced perturbations responsible for realizing such spatiotemporal heterogeneities (see Purpura et al. 2023).

Note, these analyses were now used in conjunction with the refined dry run environmental analyses (Figs. 7-8) to affirm blocking influences. The individual panels in the respective simulation folders were referenced, but given less weight to our inferences than Figs. 7-8, which are now included in the paper after completing our first round of peer-reviewed revisions.

# Folder Contents

1). CS_IDTRN_Dry_Run: 
Contains 5-min output visual analyses of a dry simulation utilizing the upstream crossing composite (e.g., Fig. 3a; Riggin et al. In Review) as the model background field when idealized terrain is included. Results support terrain-blocking impacts modifying the local environment along the windward slope of the infinitely-long idealized ridge.

2). CS_RLTRN_Dry_Run: 
Contains 5-min output visual analyses of a dry simulation utilizing the upstream crossing composite (e.g., Fig. 3a; Riggin et al. In Review) as the model background field when realistic terrain is included. Results support terrain-blocking impacts modifying the local environment along the windward slopes of the western fringe of the realistic terrain field.

3). Dry_Run_Executable: 
Contains a pre-complied executable, onefile.F and namelist for dry run simulations. If one wishes to re-run these dry simulations they must collect the appropriate input_sounding file (available in Steady-State or Variable-State folders) before running the simulation. The namelist is natively configured for a realistic terrain simulation (i.e., itern == 4), which will also require the perts.dat file to run. This can easily be changed back to the study-specific idealized terrain profile by setting itern == 1. A submit.slurm script with required dependencies is included as well.

4). NC_IDTRN_Dry_Run: 
Contains 5-min output visual analyses of a dry simulation utilizing the upstream noncrossing composite (e.g., Fig. 3d; Riggin et al. In Review) as the model background field when idealized terrain is included. Results support terrain-blocking impacts modifying the local environment along the windward slope of the infinitely-long idealized ridge.

5). CS_RLTRN_Dry_Run: 
Contains 5-min output visual analyses of a dry simulation utilizing the upstream noncrossing composite (e.g., Fig. 3d; Riggin et al. In Review) as the model background field when realistic terrain is included. Results support terrain-blocking impacts modifying the local environment along the windward slopes of the western fringe of the realistic terrain field.

6). CM1_Dry_Run_Analysis.py: 
A post-processing python3 script that calculates and generates the analyses included in each respective IDTRN dry simulation folder. The included CSTAR_Supercell_Appalachians_Modeling.yml environment file in /Scripts/ should contain all package dependencies required to run this script.

7). CM1_Dry_Run_Analysis_RLTRN.py: 
Same as in Folder 6, but for the respective RLTRN dry simulations.

# References

Purpura, S. M., C. E. Davenport, M. D. Eastin, K. E. McKeown, and R. R. Riggin, 2023: Environmental Evolution of Supercell Thunderstorms Interacting with the Appalachian Mountains. Wea. Forecasting, 38, 179–198, https://doi.org/10.1175/WAF-D-22-0115.1.

Chen, S., and Y. Lin, 2005: Effects of Moist Froude Number and CAPE on a Conditionally Unstable Flow over a Mesoscale Mountain Ridge. J. Atmos. Sci., 62, 331–350, https://doi.org/10.1175/JAS-3380.1.

Riggin IV, R. R., C. E. Davenport, M. D. Eastin, K. E. McKeown, S. M. Purpura, and B. T. Katona, In Review: Idealized Simulations of Supercell Thunderstorms Interacting with the Appalachian Mountains. Mon. Wea. Review.
