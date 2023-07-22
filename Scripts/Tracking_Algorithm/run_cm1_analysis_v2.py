#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 14:13:37 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Program: 'run_cm1_analysis_v2.py'
#
#   Update Records:
#       (3/20/22) Script created
#
#-----------------------------------------------------------------------------
#
# Summary:  
#
#   This script acts as the main driver code for data collection and analysis 
#   of a CM1 simulation. The user can defined necessary variables for tracking
#   and plotting model output. Also defines and generates an output statistic
#   CSV from the inflow soundings. It assumes that the output data is netCDF 
#   format and written in CF-Compliance which is now standard as of cm1r20.1.
#
# Notes:
#
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Import External Libraries
#-----------------------------------------------------------------------------
import os
import numpy as np
import xarray as xr
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl


#-----------------------------------------------------------------------------
# Import External User-Defined Functions
#-----------------------------------------------------------------------------
from run_cm1_analysis_initial_conditions import cm1_initial_conditions
from run_cm1_analysis_UH_tracker import cm1_uh_tracker
from run_cm1_analysis_radar_plotter import cm1_radar_plotter
from run_cm1_analysis_cross_section_plotter import cm1_cross_section_plotter
from run_cm1_analysis_sounding_plotter import cm1_sounding_plotter
from run_cm1_analysis_sounding_analysis import cm1_sounding_analysis
from run_cm1_analysis_mesoanalysis_plotter import cm1_mesoanalysis_plotter
from run_cm1_analysis_get_meso_area import cm1_get_meso_area

#-----------------------------------------------------------------------------
# User-Defined Variables 
#-----------------------------------------------------------------------------


#--------------------------------------
# Description Strings
#--------------------------------------

# Script Name
program = 'run_cm1_analysis_v2.py'

# Simulation Shortname
shortname = 'BSS_CS_RLTRN_ALL'

# Model Output Directory
# main_dir = '/scratch/rriggin/cm1r20.3/original_sims_corrected_snd/non_crosser/'+shortname.lower()+'/'
# main_dir = '/scratch/rriggin/cm1r20.3_init3d_mod/simulations/non_crosser_all/'
main_dir = '/scratch/rriggin/bss_cm1r20.3_fixed/RLTRN/' + shortname + '/'
# main_dir = '/users/rriggin/physical_proj/' + shortname + '/'

# Terrain Included?
terrain = True

# Use peak-relative coordinates?
id_terrain = False

# Horizontal Gridspacing (m)
dx  = 250.0
dy = dx
 
#--------------------------------------
# Model Configurations
#--------------------------------------

if( dx == 250.0 ):
    # Number of gridpoints in x-direction
    nx = 2400
    
    # Number of gridpoints in y-direction
    ny = 1600
    
    # Number of gridpoints in z-direction
    nz = 48
    
    skip_meso = 100


elif( dx == 500.0 ):
    # Number of gridpoints in x-direction
    nx = 1000

    # Number of gridpoints in y-direction
    ny = 1000

    # Number of gridpoints in z-direction
    nz = 48
    
    skip_meso = 50


elif( dx == 1000.0 ):
    # Number of gridpoints in x-direction
    nx = 500

    # Number of gridpoints in y-direction
    ny = 500

    # Number of gridpoints in z-direction
    nz = 48
    
    skip_meso = 25

    
# Model output time increment (min)
dt = 5.0

# Max integration time (hrs)
tot_time = 6.0

# CI start time (hrs)
ci_time = 0.0

# Turn on storm tracking algorithm time (min)
track_time = 100.0

# Output file CI begins
ci_on =  int( (ci_time * 60.0)  / dt )

# When would you like to turn on tracking algorithm? (Preferably 90 min after CI)
track_on = int( ci_on  + (track_time / dt ) )


#--------------------------------------
# Tracking Thresholds
#--------------------------------------

# Initial bounding_box indicies (Should be centered on CI location)
x1 = 0
x2 = nx-1
y1 = 0
y2 = ny-1

# Starting sounding x-index
sxi = 281 

# Starting sounding y-index
syj = 91

# Acceptable plotting distance for future time-steps (km)
bb_box_dist = 45.0  

# Distance to pull inflow sounding from updraft location (km)
prox_dist = 30.0

# Updraft Helicity threshold for identifying supercellular convection (m^2/s^2)
uh_thresh = 150.0 

# Radial distance across UH max to compute updraft statistics (km)
up_rad = 5.0

# Max distance away from UH max between timesteps (km)
overlap_dist = 15.0

# Mesocyclone area parameter threshold (m/s, 1/s, m/s^2, dBZ)
w_thresh = 10.0
zeta_thresh = 0.01 
w_zeta_thresh = 0.1
dbz_thresh = 30.0

# Moving window to compute mesocyclone area over (km^2)
meso_window = 35.0


#--------------------------------------
# Mapping Variables
#--------------------------------------

# Peak Elevation location (km)
peak_pos = 350.0 

# Peak Elevation of Terrain (m)
if( id_terrain == True ):
    peak_elv = 750.0
    
    # Terrain Contour Interval (m)
    cint_trn = 75.0
    
else:
    peak_elv = 2000.0
    
    # Terrain Contour Interval (m)
    cint_trn = 100.0

# Common Height Indices (Must change if z-stretch setting are modified)
m500 = 4
km1 = 6
km3 = 13
km5 = 18
km8 = 24
km10 = 28

#--------------------------------------
# Main Looping Variables
#--------------------------------------

# Which file would you like to start the output loop on? (Must be second output time minimum to work!)
i_start = 20

# Which file would you like to end the output loop on?
i_end = 97


#-----------------------------------------------------------------------------
# End: User-Defined Variables 
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# DO NOT CHANGE ANYTHING BELOW HERE!!!!!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# Begin Main Script
#-----------------------------------------------------------------------------


# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( "\nBegin {}...".format( program ) )



#-----------------------------------------------------------------------------
# Set-Up Working and Plotting Directories
#-----------------------------------------------------------------------------

# Change the working directory to the provided path
os.chdir( main_dir )

# Report to terminal
print( "\nChanging working directory to model output location... \n\t{}".format( main_dir ) )

plot_dir = main_dir + "/plots/"

# Create a new directory for plots if it does not exist
#------------------------------------------------------
if( os.path.isdir( plot_dir ) == False ):
    
    # Report to terminal1
    print( "\nCreating new directory to store plots...\n\t{}".format( plot_dir ) )
    
    # Create directory
    os.mkdir( plot_dir )
#------------------------------------------------------

#-----------------------------------------------------------------------------
# End Set-Up Working and Plotting Directories
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# Create and write header to output statistics CSV
#-----------------------------------------------------------------------------

# Create output stats CSV as writable file 
stats_file = open( "model_output_stats.csv", "w+" )

# Write out simulation type at top of the stats file
stats_file.write( 
                 '{}: UH_Thresh: {}, u_rad: {}, track_on: {} min'.format( shortname, uh_thresh, up_rad, track_time )
                )

# Write out column headers 
stats_file.write( 
                 '\nMode, Time, x1, y1, x2, y2, zs_surface,' + 
                 ' w500m, w1km, w3km, w5km,' + 
                 'zvort_suf, zvort500m, zvort1km, zvort3km, zvort5km,' + 
                 'UH_Max, UH_Area, UH_X, UH_Y, UHi, UHj,' +
                 'SX, SY, SXi, SYj, meanWnd, meanWnd_dir, rm, rm_dir, lm, lm_dir,' + 
                 'shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,' +
                 'srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin,' +
                 '3CAPE, LCL, LFC, LFC-LCL, EL,' +
                 'wA500m, wA1km, wA3km, wA5km, wA8km,' +
                 'zA500m, zA1km, zA3km, zA5km, zA8km,' +
                 'wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth, CP_Min, CP_area, dBZ_Area'
                )

# Write out column units   
stats_file.write( 
                 '\nbool, min, idx, idx, idx, idx, km,' +
                 'ms^-1, ms^-1, ms^-1, ms^-1,' + 
                 's^-1, s^-1, s^-1, s^-1, s^-1,' +
                 'm^2s^-2, m^2s^-2, km, km, idx, idx,' +
                 'km, km, idx, idx, ms^-1, deg, ms^-1, deg, ms^-1, deg,' + 
                 'ms^-1, deg, ms^-1, deg, ms^-1, deg,' +
                 'm^2s^-2, m^2s^-2, m^2s^-2, Jkg^-1, Jkg^-1, Jkg^-1, Jkg^-1, Jkg^-1, Jkg^-1,' +
                 'Jkg^-1, hPa, hPa, hPa, hPa,' +
                 'km^2, km^2, km^2, km^2, km^2,' +
                 'km^2, km^2, km^2, km^2, km^2,' +
                 'km^2, km^2, km^2, km^2, km^2, km, K, km^2, km^2' 
                )   

#-----------------------------------------------------------------------------
# End Create and write header to output statistics CSV
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin netCDF file loop
#------------------------------------------------------------------------------ 

# Report program status to terminal
print( "\nBegin model output loop...\n" )

# Loop through specified model output files
for i in np.arange( i_start, i_end + 1 ):
    
    # Get current integration time
    current_model_time = (i-1) * dt 
    
    #------------------------------------------------------------------------------
    # Begin Opening Current File 
    #------------------------------------------------------------------------------ 
    
    
    # Logic to open standard cm1out netCDF
    #-------------------------------------
    
    # Define the current filename as a string
    if ( i < 10 ):
        filename = 'cm1out_00000' + str(i) + '.nc'
        print( "\tOpening " + filename + "..." )
        
    elif( i < 100 and i >= 10 ):
        filename = 'cm1out_0000' + str(i) + '.nc'
        print( "\nOpening " + filename + "...")  
    
    else:
        filename = 'cm1out_000' + str(i) + '.nc'
        print( "\nOpening " + filename + "...")

    # Open the current netCDF file with xarray 
    DS = xr.open_dataset( filename, engine = "netcdf4", decode_cf = True )    
    
    #-------------------------------------

    
    # Logic to open terrain-interpolated netCDF
    #------------------------------------------    
   
    # If simulation includes terrain
    if( terrain ==  True ):
                    
        # Define the terrain-interpolated filename as a string
        if ( i < 10 ):
            filename2 = 'cm1out_00000' + str(i) + '_i.nc'
            
            # Report program status to terminal
            print( "\n\tOpening " + filename2 + "..." )
            
        elif( i < 100 and i >= 10 ):
            filename2 = 'cm1out_0000' + str(i) + '_i.nc'
            
            # Report program status to terminal
            print( "\tOpening " + filename2 + "..."  )
            
        else:
            filename = 'cm1out_000' + str(i) + '_i.nc'
            print( "\nOpening " + filename + "...")
    
        # Open the current terrain-interpolated netCDF file with xarray 
        DS_interp = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
    
    #------------------------------------------
    
    #------------------------------------------------------------------------------
    # End Opening Current File
    #------------------------------------------------------------------------------  
    
    
    #------------------------------------------------------------------------------
    # Begin Closing Current File
    #------------------------------------------------------------------------------ 
    
    # Close the current file
    DS.close()
    
    # Report program status to terminal
    print( "\n\t" + filename + " Successfully closed!"  )
    
    # If simulation includes terrain
    if( terrain == True ):
    
        # Close terrain-interpolated file    
        DS_interp.close()
        
        # Report program status to terminal
        print( "\t" + filename2 + " Successfully closed!"  )
    
    #------------------------------------------------------------------------------    
    # End Closing Current File
    #------------------------------------------------------------------------------


    #------------------------------------------------------------------------------    
    # Begin Module 1: run_cm1_analysis_get_initial_conditions
    #------------------------------------------------------------------------------ 
         
    # # Only run during first iteration    
    # if( i == i_start ):
       
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_get_initial_conditions.py"  )

    # Get storm-motion vectors from intital contions
    mean_uv, rm, lm = cm1_initial_conditions( DS, sxi, syj )

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_get_initial_conditions.py"  )
    
    #------------------------------------------------------------------------------    
    # End Module 1: run_cm1_analysis_get_initial_conditions
    #------------------------------------------------------------------------------

        
    #------------------------------------------------------------------------------    
    # Begin Module 2: run_cm1_analysis_storm_tracking
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_uh_tracker.py"  )

    ( x1, x2, y1, y2, uhx, uhy, uh_x, uh_y, sx, sy, sxi, syj,
     uh_max, uh_area, w_500m, w_1km, w_3km, w_5km, zvort_surface,
     zvort_500m, zvort_1km, zvort_3km, zvort_5km,
     classification ) = cm1_uh_tracker( 
                                       idx = i, ds = DS, nx = nx, ny = ny, nz = nz,
                                       dx = dx, dy = dy, x1 = x1, x2 = x2, y1 = y1, y2 = y2, 
                                       uh_thresh = uh_thresh, overlap_dist = overlap_dist, track_on = track_on,
                                       m500 = m500, km1 = km1, km3 = km3, km5 = km5, peak_pos = peak_pos,
                                       updraft_radius = up_rad, prox_dist = prox_dist, bb_box_dist = bb_box_dist,
                                       terrain = terrain, id_terrain = id_terrain
                                      )

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_uh_tracker.py"  )

    #------------------------------------------------------------------------------    
    # End Module 2: run_cm1_analysis_storm_tracking
    #------------------------------------------------------------------------------
    
    
    #------------------------------------------------------------------------------    
    # Begin Module 3: run_cm1_analysis_get_meso_area
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_get_meso_area.py"  )

    ( mx1, mx2, my1, my2,
      w_area_500m_total, w_area_1km_total, w_area_3km_total, w_area_5km_total, w_area_8km_total,
      zeta_area_500m_total, zeta_area_1km_total, zeta_area_3km_total, zeta_area_5km_total, zeta_area_8km_total,
      w_zeta_area_500m_total, w_zeta_area_1km_total, w_zeta_area_3km_total,
      w_zeta_area_5km_total, w_zeta_area_8km_total, meso_depth, 
      thpert_surface_area, thpert_surface_min, dbz_total ) = cm1_get_meso_area( 
                                                                    idx = i, ds = DS, nx = nx, ny = ny, nz = nx, dx = dx, dy = dy, z = km1,
                                                                    uh_x = uh_x, uh_y= uh_y, meso_window = meso_window, track_on = track_on, 
                                                                    w_thresh = w_thresh, zeta_thresh = zeta_thresh, w_zeta_thresh = w_zeta_thresh, dbz_thresh = dbz_thresh,
                                                                    m500 = m500, km1 = km1, km3 = km3, km5 = km5, km8 = km8, km10 = km10,
                                                                    terrain = terrain
                                                                   )

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_get_meso_area.py"  )

    #------------------------------------------------------------------------------    
    # End Module 2: run_cm1_analysis_storm_tracking
    #------------------------------------------------------------------------------
    
     # Report program status to terminal
    print( "\n\tBegin " + filename + " analysis!"  )
    
    #------------------------------------------------------------------------------    
    #  Begin: Set-Up Figure and Axes
    #------------------------------------------------------------------------------
    
    # Representative of the individual plots that compose the figure (Each letter is a separate subplot and . skips that position)
    figure_mosaic = """
                     AABBE
                     CCDDE
                    """ 
    
    # Create figure and axes objects
    fig, axes = plt.subplot_mosaic( 
                                   mosaic = figure_mosaic, 
                                   figsize = ( 28, 18 ), 
                                   tight_layout = True
                                  )
    
    #------------------------------------------------------------------------------    
    #  End: Set-Up Figure
    #------------------------------------------------------------------------------
    
    
    
    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_radar_plotter
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_radar_plotter.py"  )

    # Construct axis objects for the radar plot
    ax_A = fig.add_subplot( axes['A'] )

    # Run the Radar Plotter function
    cm1_radar_plotter( 
                      ds = DS, ax = ax_A, fig = fig, nx = nx, ny = ny,
                      x1 = x1, x2 = x2 , y1 = y1, y2 = y2, 
                      z = km1, sx = sx, sy = sy, uh_x = uhx, uh_y = uhy,
                      skip_val = 20, shortname = shortname, current_model_time = current_model_time,
                      terrain = terrain, peak_elv = peak_elv, 
                      cint = cint_trn, id_terrain = id_terrain,
                      meso_box = True, mx = mx1, my = my1,
                      mdx = meso_window, mdy = meso_window
                     )

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_radar_plotter.py"  )

    #------------------------------------------------------------------------------    
    # End Module X: run_cm1_analysis_radar_plotter
    #------------------------------------------------------------------------------



    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_sounding_plotter
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_sounding_plotter.py"  )
    
    # Get the original specs for subplot D
    ss = axes['B'].get_subplotspec()
    
    # Remove the original instance of subplot D
    axes['B'].remove()
    
    # Re-construct the subplot using the original specs with the skewed x-axis projection 
    axes['B'] = fig.add_subplot( ss, projection = 'skewx' )
    
    # Run the Sounding Plotter function
    cm1_sounding_plotter( 
                         ds = DS, fig = fig, ax = axes['B'],
                         sx = sxi, sy = syj, 
                         shortname = shortname, current_model_time = current_model_time,
                         knots = True, terrain = terrain, peak_rel = id_terrain
                        )

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_sounding_plotter.py"  )

    #------------------------------------------------------------------------------    
    # End Module X: run_cm1_analysis_radar_plotter
    #------------------------------------------------------------------------------    


    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_sounding_analysis
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_sounding_analysis.py"  )
    
    ( zs, lcl_pressure, lfc_pressure, el_pressure, 
      sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3km,
      shear1km, shear_1km_dir, shear3km, shear_3km_dir, shear6km, shear_6km_dir,
      rm_mag, rm_dir, lm_mag, lm_dir, mean_wind, mean_uv_dir,
      srh_500m, srh_1km, srh_3km, crit_angle ) = cm1_sounding_analysis( 
                                                                       DS, sxi, syj, shortname, current_model_time, 
                                                                       rm, lm, mean_uv, km3, terrain = terrain
                                                                      )
          
    # Report program status to terminal
    print( "\t\tEnter run_cm1_analysis_sounding_analysis.py"  )

    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_sounding_analysis
    #------------------------------------------------------------------------------



    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_cross_section_plotter
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_cross_section_plotter.py"  )

    # Construct axis objects for the radar plot
    ax_C = fig.add_subplot( axes['C'] )

    # (UH_X1 - 50 UH_X2 + 200 in original script)
    if( terrain == True ):
        cm1_cross_section_plotter( 
                                  ds_interp = DS_interp, fig = fig, ax = ax_C,
                                  x1 = x1, x2 = x2, z1 = 0, z2 = 43, nx = nx, nz = nz, uh_y = uh_y,
                                  sx = sx, sy = sy, rm_x = rm[0], rm_y  = rm[1],
                                  shortname = shortname, current_model_time = current_model_time, 
                                  skip_val = 20, terrain = terrain, knots = True,
                                  id_terrain = id_terrain, peak_pos = peak_pos
                                 )
    else:
        cm1_cross_section_plotter( 
                                  ds_interp = DS, fig = fig, ax = ax_C,
                                  x1 = x1, x2 = x2, z1 = 0, z2 = 43, nx = nx, nz = nz, uh_y = uh_y,
                                  sx = sx, sy = sy, rm_x = rm[0], rm_y  = rm[1],
                                  shortname = shortname, current_model_time = current_model_time, 
                                  skip_val = 20, terrain = terrain, knots = True,
                                  id_terrain = id_terrain, peak_pos = peak_pos
                                 )
        

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_cross_section_plotter.py"  )

    #------------------------------------------------------------------------------    
    # End Module X: run_cm1_analysis_cross_section_plotter
    #------------------------------------------------------------------------------


    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_mesoanalysis_plotter
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_mesoanalysis_plotter.py"  )

    # Construct axis objects for the radar plot
    ax_D = fig.add_subplot( axes['D'] )
    
    cm1_mesoanalysis_plotter( 
                              ds = DS, ax = ax_D, x1 = x1, x2 = x2, y1 = y1, y2 = y2, 
                              z1 = 0, z2 = nz, shortname = shortname, current_model_time =  current_model_time,
                              sx = sx, sy = sy, storm_motion = rm,
                              uh_x = uh_x, uh_y = uh_y, bb_box_dist = bb_box_dist, dx = dx, dy = dy,
                              terrain = terrain, id_terrain = id_terrain, knots = True,
                              km1 = km1, km3 = km3, km5 = km5, skip_val = skip_meso,
                              nx = nx, ny = ny, nz = nz, smooth = 10,
                              peak_pos = peak_pos, peak_elv = peak_elv, cint = cint_trn 
                            )

    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_mesoanalysis_plotter.py"  )

    #------------------------------------------------------------------------------    
    # End Module X: run_cm1_analysis_cross_section_plotter
    #-----------------------------------------------------------------------------

    # Report program status to terminal
    print( "\n\t\tWritting data to stats file..."  )
    
    # Write out data to stats file
    stats_file.write( 
                     '\n{}, {}, {}, {}, {}, {}, {},'.format( classification, current_model_time, x1, y1, x2, y2, zs.values ) + 
                       '{}, {}, {}, {},'.format( w_500m, w_1km, w_3km, w_5km ) + 
                       '{}, {}, {}, {}, {},'.format( zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km ) + 
                       '{}, {}, {}, {}, {}, {},'.format( uh_max, uh_area, uhx, uhy, uh_x, uh_y )+
                       '{}, {}, {}, {}, {}, {}, {}, {}, {}, {},'.format( sx, sy, sxi, syj, mean_wind.m, mean_uv_dir.m, rm_mag.m, rm_dir.m, lm_mag.m, lm_dir.m ) + 
                       '{}, {}, {}, {}, {}, {},'.format( shear1km.m, shear_1km_dir.m, shear3km.m, shear_3km_dir.m, shear6km.m, shear_6km_dir.m ) +
                       '{}, {}, {}, {}, {}, {}, {}, {}, {},'.format( srh_500m[0].m, srh_1km[0].m, srh_3km[0].m, sbcape.m, sbcin.m, mlcape.m, mlcin.m, mucape.m, mucin.m ) +
                       '{}, {}, {}, {}, {},'.format( cape3km.m, lcl_pressure.m, lfc_pressure.m, abs(lfc_pressure.m-lcl_pressure.m), el_pressure.m ) +
                       '{}, {}, {}, {}, {},'.format( w_area_500m_total, w_area_1km_total, w_area_3km_total, w_area_5km_total, w_area_8km_total ) +
                       '{}, {}, {}, {}, {},'.format( zeta_area_500m_total, zeta_area_1km_total, zeta_area_3km_total, zeta_area_5km_total, zeta_area_8km_total ) +
                       '{}, {}, {}, {}, {}, {},'.format( w_zeta_area_500m_total, w_zeta_area_1km_total, w_zeta_area_3km_total, w_zeta_area_5km_total, w_zeta_area_8km_total, meso_depth ) +
                       '{}, {}, {}'.format( thpert_surface_min, thpert_surface_area, dbz_total )
                    )

    # Report program status to terminal
    print( "\n\t\tFinished writting data to stats file..."  )
    
    
    #------------------------------------------------------------------------------    
    # Begin Module X: run_cm1_analysis_cross_section_plotter
    #------------------------------------------------------------------------------
    
    # Report program status to terminal
    print( "\n\t\tEnter run_cm1_analysis_stats_plotter.py"  )
    
    if( classification == 1 ):
        mode = 'Supercell'
    else:
        mode = 'Not Supercell'
        
    # Construct axis objects for the radar plot
    ax_E = fig.add_subplot( axes['E'] )
    ax_E.spines['top'].set_visible( False )
    ax_E.spines['right'].set_visible( False )
    ax_E.spines['left'].set_visible( False )
    ax_E.spines['bottom'].set_visible( False )
    ax_E.get_xaxis().set_ticks([])
    ax_E.get_yaxis().set_ticks([])
    
    ax_E.text( x = 0.05, y= 1.0, s = 'Storm Characteristics', fontsize = 18, fontweight = 'bold', transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.98, s = 'Objective Classification: ' + mode, fontsize = 16, transform = ax_E.transAxes ) 
    ax_E.text( x = 0.06, y= 0.96, s = '2-5km $UH_{max}$: ' + str( round(uh_max,2) ) +' m$^{2}$s$^{-2}$', fontsize = 16, transform = ax_E.transAxes ) 
    ax_E.text( x = 0.06, y= 0.94, s = '2-5km $UH_{avg}$: ' + str( round(uh_area,2) ) +' m$^{2}$s$^{-2}$', fontsize = 16, transform = ax_E.transAxes ) 
    
    ax_E.text( x = 0.05, y= 0.90, s = 'Inflow Thermodynamics', fontsize = 18, fontweight = 'bold', transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.88, s = 'SBCAPE: ' + str( round(sbcape.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.86, s = 'SBCIN: ' + str( round(sbcin.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.84, s = 'MLCAPE: ' + str( round(mlcape.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.82, s = 'MLCIN: ' + str( round(mlcin.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.80, s = 'MUCAPE: ' + str( round(mucape.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.78, s = 'MUCIN: ' + str( round(mucin.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.76, s = '3CAPE: ' + str( round(cape3km.m,2) ) +' Jkg$^{-1}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.74, s = 'LCL: ' + str( round(lcl_pressure.m,2) ) +' hPa', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.72, s = 'LFC: ' + str( round(lfc_pressure.m,2) ) +' hPa', fontsize = 16, transform = ax_E.transAxes )        
    ax_E.text( x = 0.06, y= 0.70, s = 'LFC-LCL: ' + str( round( abs( lfc_pressure.m - lcl_pressure.m ), 2) ) +' hPa', fontsize = 16, transform = ax_E.transAxes )
    
    ax_E.text( x = 0.05, y= 0.66, s = 'Inflow Kinematics', fontsize = 18, fontweight = 'bold', transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.64, s = '0-1 km Shear: ' + str( np.round(shear1km.m, 2 ) ) +' ms$^{-1}$ ' + str( np.round(shear_1km_dir.m, 2 ) ) + '\u00b0', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.62, s = '0-3 km Shear: ' + str( np.round(shear3km.m, 2 ) ) +' ms$^{-1}$ ' + str( np.round(shear_3km_dir.m, 2 ) ) + '\u00b0', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.60, s = '0-6 km Shear: ' + str( np.round(shear6km.m, 2 ) ) +' ms$^{-1}$ ' + str( np.round(shear_6km_dir.m, 2 ) ) + '\u00b0', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.58, s = 'Mean Wind: ' + str( np.round(mean_wind.m, 2 ) ) +' ms$^{-1}$ ' + str( np.round(mean_uv_dir.m, 2 ) ) + '\u00b0', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.56, s = "Bunker's Right Motion: " + str( np.round(rm_mag.m, 2 ) ) +' ms$^{-1}$ ' + str( np.round(rm_dir.m, 2 ) ) + '\u00b0', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.54, s = "Bunker's Left Motion: " + str( np.round(lm_mag.m, 2 ) ) +' ms$^{-1}$ ' + str( np.round(lm_dir.m, 2 ) ) + '\u00b0', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.52, s = "0-500 m SRH: " + str( round(srh_500m[0].m, 2 ) ) +' m$^{2}$s$^{-2}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.50, s = "0-1 km SRH: " + str( round(srh_1km[0].m, 2 ) ) +' m$^{2}$s$^{-2}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.48, s = "0-3 km SRH: " + str( round(srh_3km[0].m, 2 ) ) +' m$^{2}$s$^{-2}$', fontsize = 16, transform = ax_E.transAxes )    
    
    ax_E.text( x = 0.05, y= 0.44, s = '{}-km UH-Centric Areal Averages'.format( up_rad ), fontsize = 18, fontweight = 'bold', transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.42, s = '$W_{500m}$: ' + str( round(w_500m, 2 ) )  + ' ms$^{-1}$', fontsize = 16, transform = ax_E.transAxes )     
    ax_E.text( x = 0.06, y= 0.40, s = '$W_{1km}$: ' + str( round(w_1km, 2 ) ) + ' ms$^{-1}$', fontsize = 16, transform = ax_E.transAxes )   
    ax_E.text( x = 0.06, y= 0.38, s = '$W_{3km}$: ' + str( round(w_3km, 2 ) ) + ' ms$^{-1}$', fontsize = 16, transform = ax_E.transAxes )    
    ax_E.text( x = 0.06, y= 0.36, s = '$W_{5km}$: ' + str( round(w_5km, 2 ) ) + ' ms$^{-1}$', fontsize = 16, transform = ax_E.transAxes ) 
    ax_E.text( x = 0.06, y= 0.34, s = '$\u03b6_{surface}$: ' + str( round(zvort_surface, 6 ) ) + ' s$^{-1}$', fontsize = 16, transform = ax_E.transAxes )   
    ax_E.text( x = 0.06, y= 0.32, s = '$\u03b6_{500m}$: ' + str( round(zvort_500m, 5 ) ) + ' s$^{-1}$', fontsize = 16, transform = ax_E.transAxes )   
    ax_E.text( x = 0.06, y= 0.30, s = '$\u03b6_{1km}$: ' + str( round(zvort_1km, 5 ) ) + ' s$^{-1}$', fontsize = 16, transform = ax_E.transAxes )   
    ax_E.text( x = 0.06, y= 0.28, s = '$\u03b6_{3km}$: ' + str( round(zvort_3km, 5 ) ) + ' s$^{-1}$', fontsize = 16, transform = ax_E.transAxes )   
    ax_E.text( x = 0.06, y= 0.26, s = '$\u03b6_{5km}$: ' + str( round(zvort_5km, 5 ) ) + ' s$^{-1}$', fontsize = 16, transform = ax_E.transAxes ) 
    
    ax_E.text( x = 0.05, y= 0.22, s = 'Mesocyclone Area Metrics', fontsize = 18, fontweight = 'bold', transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.20, s = '$W_{500m}$: ' + str( round( w_area_500m_total, 4 ) ) + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.18, s = '$W_{1km}$: ' + str( round( w_area_1km_total, 4 ) ) + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.16, s = '$W_{3km}$: ' + str( round( w_area_3km_total, 4 ) ) + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.14, s = '$W_{5km}$: ' + str( round( w_area_5km_total, 4 ) )  + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.12, s = '$W_{8km}$: ' + str( round( w_area_8km_total, 4 ) )  + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.10, s = '$W\u03b6_{500m}$: ' + str( round( w_zeta_area_500m_total, 4 ) ) + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.08, s = '$W\u03b6_{1km}$: ' + str( round( w_zeta_area_1km_total, 4 ) ) + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.06, s = '$W\u03b6_{3km}$: ' + str( round( w_zeta_area_3km_total, 4 ) ) + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.04, s = '$W\u03b6_{5km}$: ' + str( round( w_zeta_area_5km_total, 4 ) )  + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.02, s = '$W\u03b6_{8km}$: ' + str( round( w_zeta_area_8km_total, 4 ) )  + ' km$^{2}$', fontsize = 16, transform = ax_E.transAxes )
    ax_E.text( x = 0.06, y= 0.00, s = 'Mesocyclone Depth: ' + str( round( meso_depth, 2 ) )  + ' km', fontsize = 16, transform = ax_E.transAxes )
    
    # Report program status to terminal
    print( "\t\tExit run_cm1_analysis_stats_plotter.py"  )

    #------------------------------------------------------------------------------    
    # End Module X: run_cm1_analysis_cross_section_plotter
    #------------------------------------------------------------------------------
    
    
    #------------------------------------------------------------------------------    
    # Begin: Create Terrain Colorbar
    #------------------------------------------------------------------------------ 
    if( terrain == True ):
        cb_ax = fig.add_axes( [0.015, -0.05, 0.975, 0.025] )
        cmap = mpl.cm.copper.reversed()
        bounds = np.arange( 0, peak_elv + cint_trn , cint_trn ) 
        norm = mpl.colors.BoundaryNorm( bounds, cmap.N )
        cbar = fig.colorbar( mpl.cm.ScalarMappable( norm = norm, cmap = cmap ),
                             cax = cb_ax, orientation = 'horizontal', alpha = 0.35 )
        cbar.set_label( label = 'Elevation (m)', fontsize = 20, fontweight = 'bold' )
        cb_ax.tick_params( axis = 'x', labelsize = 16 )
    #------------------------------------------------------------------------------    
    # End: Create Terrain Colorbar
    #------------------------------------------------------------------------------ 
            
    
    #------------------------------------------------------------------------------
    # Begin Save Figure
    #------------------------------------------------------------------------------
    
    plt.show()
    
    # Change the working directory to the plot folder 
    os.chdir( plot_dir )
    
    # Save the current RH plot
    fig.savefig(
                fname = '{}_4panel_{}_min_{}.jpeg'.format( shortname, current_model_time, i ),
                dpi = 300,
                bbox_inches = "tight"
               )

    # # Clear Current Figure for re-use on next iteration
    # fig.clf()
    
    # Change the working directory back to the model output folder
    os.chdir( main_dir )

    #------------------------------------------------------------------------------
    # End Save Figure
    #------------------------------------------------------------------------------

    # Report program status to terminal
    print( "\n\t" + filename + " analysis complete!\n"  )

#-----------------------------------------------------------------------------
# End netCDF file loop
#------------------------------------------------------------------------------

# Report program status to terminal
print( "\nEnd of model output loop..." )

# Close the output stats CSV
stats_file.close()

#-----------------------------------------------------------------------------
# End Main Script
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# DO NOT CHANGE ANYTHING ABOVE HERE!!!!!
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
   


# Confirm that script successfully ran
print( "\n{} successfully completed!".format( program ) )

# Report the time required to run the function
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )
