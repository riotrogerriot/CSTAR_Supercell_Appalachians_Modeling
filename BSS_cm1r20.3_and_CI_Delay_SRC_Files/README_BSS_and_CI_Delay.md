# README_CI_DELAY

README file for how to incorporate CI Delay in CM1 when using Updraft Nudging Technique (Naylor and Gilmore 2012)
Author: R.R. Riggin, July 2023

NOTE: Code modifications only designed to delay CI with the updraft nudging technique (i.e., iinit = 12) and has only been verified with CM1r20.3. I am unsure if these changes translate over to CM1r.21.0 (Latest) smoothly.

The package for incorporating CI Delay contains several files of interest
1). init3d.F including modifications to set start time for updraft nudging (i.e., t0_wnudge)
2). input.F including modification to declare nudging start time var (i.e., t0_wnudge)
3). Sound.F including modifications using if/else logic to delay CI via updraft nudging until integration time == t0_wnudge

All modifications can be found by searching "RR added" within the respective fortran file

How to incorporate CI_Delay files:
1. Save copies of any duplicate files listed above you may already have for CM1 (i..e, cm1.F and Makefile)
   The native CM1r20.3 code can be found at https://www2.mmm.ucar.edu/people/bryan/cm1/
2. Move the files noted above into the cm1r20.3/src folder
3. Define updraft nudging parameters, including start time via t0_wnudge, in init3d.F
4. When ready to compile, return to src folder, type "make"
5. Configure namelist as usual and run executable to start simulation

In Riggin IV et al. (202X), this modification required to run the steady-state RLTRN simulations (e.g., CS_RLTRN_ALL, NC_RLTRN_ALL, CS_RLTRN_NO_OUT, NC_RLTRN_NO_OUT, CS_RLTRN_OUT, NC_RLTRN_OUT). This modification was used in conjunction with BSS (see below) to run the variable-state RLTRN simulations (e.g., BSS_CS_RLTRN_ALL, BSS_NC_RLTRN_ALL).


# README_BSS

README file for how to run continuous BSS in CM1 (updated for r20.3)
Author: C.E. Davenport, June 2023

Refer to Letkewicz and Parker (2013) and Davenport et al. (2019) for details regarding BSS methodology

The package for incorporating continous BSS contains several files of interest:
1. Makefile to compile; note that compiler selection is specific for UNC Charlotte HPC system
2. base_bss.F;
3. bss_tend.F;
4. cm1.F; updated to include variable declarations and calls to bss subroutines


How to incorporate BSS files: 
1. Save copies of any duplicate files listed above you may already have for CM1 (i..e, cm1.F and Makefile)
   The native CM1r20.3 code can be found at https://www2.mmm.ucar.edu/people/bryan/cm1/
2. Move all files to the cm1r20.3/src folder; the submit script should go in cm1r20.3/run
3. Load mpi compiler (I was able to successfully compile with "module load openmpi/3.1.2-pgi" and settings found in Makefile)
4. In src folder, type "make"


How to run CM1 with BSS:
1. Set up model as you desire (grid spacings, microphysics, initiation type, etc) within the namelist file
2. Within your namelist, go to the "param8" section, which contains flex variables. Change settings according to how you want to run BSS. Currently hardwired to apply tendencies between initial sounding (found in the "input_sounding" file in run folder) and two additional soundings. In other words, tendencies are computed between 1st (original input) and 2nd soundings, then between 2nd and 3rd soundings. It's a bit hacked, yes, but it works.

var1: BSS start time (in sec)

var2: BSS end time (in sec)

var3: Time to implement second sounding (model run time, since start of model run)

var4: Time to implement third sounding (in current formulation, should = BSS end time)

var5: Change wind profile? (1=yes,0=no)

var6: Change thermo profile? (1=yes,0=no)

3. Modify the "replace" files found in the run directory, using the format provided. There are separate ones for the thermo and wind profiles. "Replace1" corresponds to the first sounding that replaces the original "input_sounding", then "replace2" corresponds to the second sounding that replaces the prior one.
4. Update the submit script (runcm1.sh) with appropriate processor, walltime, etc. info
5. Submit the run!

Note: There are a couple of common errors I've done that are good to check for before running: 1) Make sure all soundings (input and replacement soundings) are stable (e.g., see checks in base.F that occur). If there are weird layers where N<0 or similar, that'll cause issues. 2) Make sure all soundings have depths equal to or greater than the depth of the grid. Model needs sufficient data points!

# References

Davenport, Casey E., Ziegler, Conrad L., and Biggerstaff, Michael I., 2019, "Creating a More Realistic Idealized Supercell Thunderstorm Evolution via Incorporation of Base-State Environmental Variability" Monthly Weather Review Vol. 147, No. 11, pp 4177, 1520-0493, https://doi.org/10.1175/MWR-D-18-0447.1.

Letkewicz, C. E., A. J. French, and M. D. Parker, 2013: Base-State Substitution: An Idealized Modeling Technique for Approximating Environmental Variability. Mon. Wea. Rev., 141, 3062–3086, https://doi.org/10.1175/MWR-D-12-00200.1.

Naylor, J., and M. S. Gilmore, 2012: Convective Initiation in an Idealized Cloud Model Using an Updraft Nudging Technique. Mon. Wea. Rev., 140, 3699–3705, https://doi.org/10.1175/MWR-D-12-00163.1.

Riggin IV, R. R., C. E. Davenport, M. D. Eastin, K. E. McKeown, S. M. Purpura, and B. T. Katona, 202X: Idealized Simulations of Supercell Thunderstorms Interacting with the Appalachian Mountains. In Prep.
