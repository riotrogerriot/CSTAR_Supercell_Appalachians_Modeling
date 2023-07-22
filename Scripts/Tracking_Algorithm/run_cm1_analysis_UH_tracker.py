#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 13:09:43 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Function: cm1_UH_tracker.py
#
#   Update Records:
#       (4/2/22) Script Created
#
# Notes:
#
#
#-----------------------------------------------------------------------------
#
# Summary:  
#
# 
# Arguments:
#
#
# Returns:
#
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Load external libraries
#-----------------------------------------------------------------------------
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import scipy
from matplotlib.colors import Normalize
import metpy.plots as plots
import matplotlib.patheffects as PathEffects
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import pint
import pint_xarray
import cf_xarray
# from datetime import datetime
# from run_cm1_analysis_radar_plotter import cm1_radar_plotter
# from run_cm1_analysis_sounding_plotter import cm1_sounding_plotter
# from run_cm1_analysis_initial_conditions import cm1_initial_conditions
# from run_cm1_analysis_cross_section_plotter import cm1_cross_section_plotter
# from run_cm1_analysis_mesoanalysis_plotter import cm1_mesoanalysis_plotter


#-----------------------------------------------------------------------------
# Begin: Define cm1_uh_tracker
#-----------------------------------------------------------------------------
def cm1_uh_tracker( 
                   idx, ds, nx, ny, nz, dx, dy,  x1, x2, y1, y2, 
                   uh_thresh = 200.0, overlap_dist = 15.0, track_on = 19,
                   m500 = 4, km1 = 6, km3 = 13, km5 = 18, peak_pos = 350.0,
                   updraft_radius = 5.0, prox_dist = 30.0, bb_box_dist = 40.0,
                   terrain = False, id_terrain = False
                  ):

    # Store prior window for overlap checking
    x1_past = x1
    x2_past = x2
    y1_past = y1
    y2_past = y2
    overlap_dist = overlap_dist * units.km
    
    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the 2-5 km Updraft Helicity (m^2/s^2)
    uh = ds.metpy.parse_cf( 'uh' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the w-component of the wind (m/s)   
    w = ds.metpy.parse_cf( 'winterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the vertical vorticity (1/s)   
    zvort = ds.metpy.parse_cf( 'zvort' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Drop units to speed up looping
    uh = uh.values
    w = w.values
    zvort = zvort.values
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get the x-axis coordinates (km)
    xh = ds.coords[ 'xh' ].isel( xh = slice( 0, nx ) ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)      
    yh = ds.coords[ 'yh' ].isel( yh = slice( 0, ny ) ).values * units( 'km' )
    
    # Get the z-axis coordinates (km)      
    zh = ds.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' )
    
    #-----------------------------------------------------------------------------
    # End: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Identify UH Maxima Location
    #-----------------------------------------------------------------------------

    # Find UH abs maxima
    uh_max = np.amax( uh )
    
    # Get location of UH max
    uh_max_loc = np.where( uh == uh_max )
    uh_y = int( uh_max_loc[0][0] )
    uh_x = int( uh_max_loc[1][0] )
    
    # Store locational indices
    uhx = xh[uh_x].m
    uhy = yh[uh_y].m
    
    # Report results back to terminal
    print( '\t\tUHmax: {} \tUHx: {} \tUHy: {}'.format( uh_max, xh[uh_x], yh[uh_y] ) )
    
    #-----------------------------------------------------------------------------
    # End: Identify UH Maxima Location
    #-----------------------------------------------------------------------------   
    
    
    #-----------------------------------------------------------------------------
    # Begin: Generate BB_Box Coordinates Centered on Storm
    #-----------------------------------------------------------------------------
    
    x1 = int( uh_x - ( ( 0.50 * bb_box_dist * 1000.0 ) / dx ) )
    x2 = int( uh_x + ( ( 1.50 * bb_box_dist * 1000.0 ) / dx ) )
    y1 = int( uh_y - ( ( 0.75 * bb_box_dist * 1000.0 ) / dy ) )
    y2 = int( uh_y + ( ( 1.25 * bb_box_dist * 1000.0 ) / dy ) )
    
    x1_dist = 0
    x2_dist = 0
    y1_dist = 0
    y2_dist = 0
    
    #-----------------------------------------------------------------------------
    # End: Generate BB_Box Coordinates Centered on Storm
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Overlap Check
    #-----------------------------------------------------------------------------
    
    # Turn on supercell tracking at
    if( idx >= track_on ):

	#print( '\nTracking Algorithm Activated' )
        
        # Overlap Boolean
        overlap = 0
        
        #-----------------------------------------------------------------------------
        # Overlap logic
        while( overlap == 0 ):
        
            # Report to terminal
            print( '\nChecking for overlap...' )
            
            # Is the new indices within a reasonable distance of the old ones?
            if( x1 < nx and x2 < nx and y1 < ny and y2 < ny ):
                
                x1_dist = abs( xh[x1] - xh[x1_past] ) 
                x2_dist = abs( xh[x2] - xh[x2_past] ) 
                y1_dist = abs( yh[y1] - yh[y1_past] )
                y2_dist = abs( yh[y2] - yh[y2_past] )
                
                print( 'x1_dist: {} \tx2_dist: {}'.format( x1_dist, x2_dist ) )
                print( 'y1_dist: {} \ty2_dist: {}'.format( y1_dist, y2_dist ) )
            
            else:
                      
                x1 = x1
                x2 = nx - 1
                y1 = y1
                y2 = ny - 1
                
                x1_dist = abs( xh[x1] - xh[x1_past] ) 
                x2_dist = abs( xh[x2] - xh[x2_past] ) 
                y1_dist = abs( yh[y1] - yh[y1_past] )
                y2_dist = abs( yh[y2] - yh[y2_past] )
                
            # Check to ensure the new indices are within the previous window
            if( x1_dist > overlap_dist or x2_dist > overlap_dist 
                or y1_dist > overlap_dist or y2_dist > overlap_dist ):
                
                # Report to terminal
                print( '\n\tOverlap check failed finding UH max within prior BB_Box ...' )
                
                # Find maximum UH within provided overlap dist
                uh_max = np.amax( uh[y1_past:y2_past, x1_past:x2_past] )
            
                # Get indices of the max UH location
                uh_max_loc = np.where( uh == uh_max )
                uh_y = int( uh_max_loc[0][0] )
                uh_x = int( uh_max_loc[1][0] )  
                
                # Report results back to terminal
                print( '\t\tAdjusted UHmax: {} \tUHx: {} \tUHy: {}'.format( uh_max, xh[uh_x], yh[uh_y] ) )
                
                #-----------------------------------------------------------------------------
                # Begin: Correct BB_Box Coordinates Centered on Storm
                #-----------------------------------------------------------------------------
                
                x1 = int( uh_x - ( ( 0.50 * bb_box_dist * 1000.0 ) / dx ) )
                x2 = int( uh_x + ( ( 1.50 * bb_box_dist * 1000.0 ) / dx ) )
                y1 = int( uh_y - ( ( 0.75 * bb_box_dist * 1000.0 ) / dy ) )
                y2 = int( uh_y + ( ( 1.25 * bb_box_dist * 1000.0 ) / dy ) )
                
                # Is the new indices within a reasonable distance of the old ones?
                if( x1 < nx and x2 < nx and y1 < ny and y2 < ny ):
                
                    x1 = x1
                    y1 = y1
                    x2 = x2
                    y2 = y2
                    
                    # print( '\nmod_x1: {} \tmod_x2: {}'.format( x1, x2 ) )
                    # print( 'mod_y1: {} \tmod_y2: {}'.format( y1, y2) )
                
                else:
                    # # Report to terminal
                    # print( '\n Attempting to move out of domain auto-correction being applied...' )
                    x1 = x1
                    x2 = nx - 1
                    y1 = y1
                    y2 = ny - 1
                
                #-----------------------------------------------------------------------------
                # End: Correct BB_Box Coordinates Centered on Storm
                #-----------------------------------------------------------------------------
                
                # Break the loop now that the UH tracking has been corrected
                overlap = 1
                
            else:
                overlap = 1
            
        # Store locational indices
        uhx = xh[uh_x].m
        uhy = yh[uh_y].m
        
        #-----------------------------------------------------------------------------
        # End: Overlap Check
        #-----------------------------------------------------------------------------
        
    
    #-----------------------------------------------------------------------------
    # Begin: Generate Proximity Sounding Location
    #-----------------------------------------------------------------------------
    
    # Pull Sounding 30 km SE of UH max (prox_dist should be 1/2 of requested distance [a^2 + b^2 = c^2] )
    sx_dist = ( ( prox_dist / 3.0 ) * 1000.0 ) / dx
    sy_dist = ( ( prox_dist / 1.0 ) * 1000.0 ) / dx   
    
    if( sx_dist < nx and sy_dist < ny ):
        sx = uh_x + sx_dist
        sy = uh_y - sy_dist
        
        # Ensure indices are in integer format
        sx = int(sx)
        sy = int(sy)
        
        # Save the indices
        sxi = sx   
        syj = sy 
        
        # Get the grid-relative coordinates for proximity sounding (Required for clean plotting w/ Xarray)
        sx = xh[ sx ].m
        sy = yh[ sy ].m
    
    #-----------------------------------------------------------------------------
    # End: Generate Proximity Sounding Location
    #-----------------------------------------------------------------------------
    

    #-----------------------------------------------------------------------------
    # Begin: Compute areal averaging about the UH maxima
    #-----------------------------------------------------------------------------
    
    # Define the search radius about the uh max as a callable integer value for looping purposes (Units: m)
    uh_radius = int( ( updraft_radius * 1000.0 ) / dx )
    
    # Initialize areal averaged values 
    num_gridpoints = 0
    uh_area = 0.0
    w_500m = 0.0
    w_1km = 0.0
    w_3km = 0.0
    w_5km = 0.0
    zvort_surface = 0.0
    zvort_500m = 0.0
    zvort_1km = 0.0
    zvort_3km = 0.0
    zvort_5km = 0.0
    
    # Get the looping indices 
    i_start = uh_x - uh_radius
    i_end = uh_x + uh_radius
    j_start = uh_y - uh_radius
    j_end = uh_y + uh_radius
    
    # Begin: Loop through pre-defined updraft diameter in the x-direction
    #-----------------------------------------------------------------------------
    for i in range( i_start, i_end ):
    
        # Begin: Loop through pre-defined updraft diameter in the y-direction
        #-----------------------------------------------------------------------------
        for j in range( j_start, j_end ):
            
            # If the radius lies within the domain
            if( i < nx and j < ny):
                
                # Sum up the areal values of variables within the updraft
                num_gridpoints = num_gridpoints +  1 * (dx * dy)
                uh_area = uh_area + uh[ j, i ] * (dx * dy)
                w_500m = w_500m + w[ m500, j, i ] * (dx * dy)
                w_1km = w_1km + w[ km1, j, i ] * (dx * dy)
                w_3km = w_3km + w[ km3, j, i ] * (dx * dy)
                w_5km = w_5km + w[ km5, j, i ] * (dx * dy)
                zvort_surface = zvort_surface + zvort[ 0, j, i ] * (dx * dy)
                zvort_500m = zvort_500m + zvort[ m500, j, i ] * (dx * dy)
                zvort_1km = zvort_1km + zvort[ km1, j, i ] * (dx * dy)
                zvort_3km = zvort_3km + zvort[ km3, j, i ] * (dx * dy)
                zvort_5km = zvort_5km +  zvort[ km5, j, i ] * (dx * dy)
            
            # End the loop if the updraft moves outside of the domain (Not likely but needs to be explicitly declared to prevent code from breaking)    
            else:
                break
            
        # End: Loop through pre-defined updraft diameter in the y-direction
        #-----------------------------------------------------------------------------
    
    # End: Loop through pre-defined updraft diameter in the x-direction
    #-----------------------------------------------------------------------------
    
    # Compute the areal average of updraft variables 
    uh_area = uh_area / float( num_gridpoints )
    w_500m = w_500m / float( num_gridpoints )
    w_1km = w_1km / float( num_gridpoints )
    w_3km = w_3km / float( num_gridpoints )
    w_5km = w_5km / float( num_gridpoints )
    zvort_surface = zvort_surface / float( num_gridpoints )
    zvort_500m = zvort_500m / float( num_gridpoints )
    zvort_1km = zvort_1km / float( num_gridpoints )
    zvort_3km = zvort_3km / float( num_gridpoints )
    zvort_5km = zvort_5km / float( num_gridpoints )    
    
    #-----------------------------------------------------------------------------
    # End: Compute areal averaging about the UH maxima
    #-----------------------------------------------------------------------------

        
    #-----------------------------------------------------------------------------
    # Begin: Supercell Classification
    #-----------------------------------------------------------------------------
    
    # Check to see if areal UH exceeds classification threshold
    if( uh_area >= uh_thresh ):
        
        # Boolean classification (Supercell)
        storm_class = 1
        
    else:
        
        # Boolean classification (Not Supercell)
        storm_class = 0
    
    # Store classification into array for later usage
    classification = storm_class
    
    #-----------------------------------------------------------------------------
    # End: Supercell Classification
    #-----------------------------------------------------------------------------
    
    
    # Return...
    return(
           x1, x2, y1, y2, uhx, uhy, uh_x, uh_y, sx, sy, sxi, syj,
           uh_max, uh_area, w_500m, w_1km, w_3km, w_5km, zvort_surface,
           zvort_500m, zvort_1km, zvort_3km, zvort_5km, classification
          )

#-----------------------------------------------------------------------------
# End: Define cm1_uh_tracker
#-----------------------------------------------------------------------------

















# # Record Script start time
# startTime = datetime.now()

# nx = 2400
# ny = 1000
# nz = 48

# terrain = True 
# id_terrain = True
# peak_pos = 350.0

# # Initial bounding_box indicies (Should be centered on CI location)
# x1 = 0
# x2 = nx-1
# y1 = 0
# y2 = ny-1

# dx = 250.0
# dy = 250.0

# sx = 250
# sy = 500

# # Common Height Indices (Must change if z-stretch setting are modified)
# m500 = 4
# km1 = 6
# km3 = 13
# km5 = 18

# updraft_radius = 5.0    # km ( Should range between 1.5-2.5 based on McKeown 2021 w/ CS diameter medians ~ 5 and NC ~ 3 km)
# uh_thresh = 200.0       # m^2/s^2
# prox_dist = 30.0        # km

# overlap_dist = 15.0*units.km     # TEST OTHER SETTINGS BE SURE TO NOTE FINAL CHOICE IN PAPER
# bb_box_dist = 40.0      # km

# # Output file Tracking begins
# track_on =  int( 90.0 / 5.0 ) + 1


# # Empty arrays to store storm loc
# uhx_idx = np.zeros( 74 )
# uhy_idx = np.zeros( 74 )
# uhx_grid = np.zeros( 74 )
# uhy_grid = np.zeros( 74 )

# uhwx_idx = np.zeros( 74 )
# uhwy_idx = np.zeros( 74 )
# uhwx_grid = np.zeros( 74 )
# uhwy_grid = np.zeros( 74 )

# classification = np.zeros( 74 )

# x1_past = 0
# x2_past = nx-1
# y1_past = 0
# y2_past = ny-1

# shortname = 'NC_TRN'

# # Loop through 90-360 min of integration (AKA when tracking == on)
# for k in range( 55, 74 ):
    

#     # Model output directory
#     directory = '/scratch/rriggin/cm1r20.3/original_sims_corrected_snd/non_crosser/nc_trn/'
    
#     # Get current integration time
#     current_model_time = (k-1) * 5.0
    
#     # Define the current filename as a string
#     if ( k < 10 ):
#         filename = directory + 'cm1out_00000' + str(k) + '.nc'
        
#         # Report program status to terminal
#         print( "\tOpening " + filename + "..." )
        
#     else:
#         filename = directory + 'cm1out_0000' + str(k) + '.nc'
        
#         # Report program status to terminal
#         print( "\tOpening " + filename + "..." )
    
#     # Open the current terrain-interpolated netCDF file with xarray
#     ds = xr.open_dataset( filename, engine = "netcdf4", decode_cf = True )
    
#     # If simulation includes terrain
#     if( terrain ==  True ):
                    
#         # Define the terrain-interpolated filename as a string
#         if ( k < 10 ):
#             filename2 = directory + 'cm1out_00000' + str(k) + '_i.nc'
            
#             # Report program status to terminal
#             print( "\n\tOpening " + filename2 + "..." )
            
#         else:
#             filename2 = directory +  'cm1out_0000' + str(k) + '_i.nc'
            
#             # Report program status to terminal
#             print( "\tOpening " + filename2 + "..."  )
    
#         # Open the current terrain-interpolated netCDF file with xarray 
#         ds_interp = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
    
    
#     #-----------------------------------------------------------------------------
#     # Begin: Get variables from xarray dataset
#     #-----------------------------------------------------------------------------
    
#     # Data Variables
#     #-----------------------------------------------------------------------------
    
#     # Get the 2-5 km Updraft Helicity (m^2/s^2)
#     uh = ds.metpy.parse_cf( 'uh' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
#     # Get the w-component of the wind (m/s)   
#     w = ds.metpy.parse_cf( 'winterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
#     # Get the vertical vorticity (1/s)   
#     zvort = ds.metpy.parse_cf( 'zvort' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
#     # Drop units to speed up looping
#     uh = uh.values
#     w = w.values
#     zvort = zvort.values
    
#     # Coordinate Variables
#     #-----------------------------------------------------------------------------
    
#     # Get the x-axis coordinates (km)
#     xh = ds.coords[ 'xh' ].isel( xh = slice( 0, nx ) ).values * units( 'km' )
    
#     # Get the y-axis coordinates (km)      
#     yh = ds.coords[ 'yh' ].isel( yh = slice( 0, ny ) ).values * units( 'km' )
    
#     # Get the z-axis coordinates (km)      
#     zh = ds.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' )
    
#     # Logic to determine if x-axis should be peak-relative or not
#     #------------------------------------------------------------
#     if( id_terrain == True ):
        
#         # Convert x-coordinates from grid-relative to peak-relative coordinates
#         xh = xh - peak_pos * units( 'km' )
    
#     #-----------------------------------------------------------------------------
#     # End: Get variables from xarray dataset
#     #-----------------------------------------------------------------------------
    
    
#     #-----------------------------------------------------------------------------
#     # Begin: Identify UH Maxima Location
#     #-----------------------------------------------------------------------------

#     # Find UH abs maxima
#     uh_max = np.amax( uh )
    
#     # Get location of UH max
#     uh_max_loc = np.where( uh == uh_max )
#     uh_y = int( uh_max_loc[0][0] )
#     uh_x = int( uh_max_loc[1][0] )
    
#     # Report results back to terminal
#     print( '\tUHmax: {} \tUHx: {} \tUHy: {}'.format( uh_max, xh[uh_x], yh[uh_y] ) )
    
#     #-----------------------------------------------------------------------------
#     # End: Identify UH Maxima Location
#     #-----------------------------------------------------------------------------    
    
    
#     #-----------------------------------------------------------------------------
#     # Begin: Generate BB_Box Coordinates Centered on Storm
#     #-----------------------------------------------------------------------------
    
#     x1 = int( uh_x - ( ( 0.50 * bb_box_dist * 1000.0 ) / dx ) )
#     x2 = int( uh_x + ( ( 1.50 * bb_box_dist * 1000.0 ) / dx ) )
#     y1 = int( uh_y - ( ( 0.75 * bb_box_dist * 1000.0 ) / dy ) )
#     y2 = int( uh_y + ( ( 1.25 * bb_box_dist * 1000.0 ) / dy ) )
    
#     print( '\nx1: {} \tx2: {}'.format( x1, x2 ) )
#     print( '\ny1: {} \ty2: {}'.format( y1, y2) )
    
#     #-----------------------------------------------------------------------------
#     # End: Generate BB_Box Coordinates Centered on Storm
#     #-----------------------------------------------------------------------------
    
    
#     #-----------------------------------------------------------------------------
#     # Begin: Overlap Check
#     #-----------------------------------------------------------------------------
    
#     # Turn on supercell tracking at 
#     if( k >= track_on ):
        
#         # Overlap Boolean
#         overlap = 0
        
#         #-----------------------------------------------------------------------------
#         # Overlap logic
#         while( overlap == 0 ):
        
#             # Report to terminal
#             print( '\nChecking for overlap...' )
            
#             # Is the new indices within a reasonable distance of the old ones?
#             if( x1 < nx and x2 < nx and y1 < ny and y2 < ny ):
                
#                 x1_dist = abs( xh[x1] - xh[x1_past] ) 
#                 x2_dist = abs( xh[x2] - xh[x2_past] ) 
#                 y1_dist = abs( yh[y1] - yh[y1_past] )
#                 y2_dist = abs( yh[y2] - yh[y2_past] )
                
#                 print( 'x1_dist: {} \tx2_dist: {}'.format( x1_dist, x2_dist ) )
#                 print( 'y1_dist: {} \ty2_dist: {}'.format( y1_dist, y2_dist ) )
            
#             else:
#                 # Report to terminal
#                 print( 'Attempting to move out of domain auto-correction being applied...' )
#                 x1 = x1
#                 x2 = nx - 1
#                 y1 = y1
#                 y2 = ny - 1
                
#             # Check to ensure the new indices are within the previous window
#             if( x1_dist > overlap_dist or x2_dist > overlap_dist 
#                 or y1_dist > overlap_dist or y2_dist > overlap_dist ):
                
#                 # Report to terminal
#                 print( '\Overlap check failed finding UH max within prior BB_Box ...' )
                
#                 # Find maximum UH within provided overlap dist
#                 uh_max = np.amax( uh[ y1_past:y2_past, x1_past:x2_past ] )
            
#                 # Get indices of the max UH location
#                 uh_max_loc = np.where( uh == uh_max )
#                 uh_y = int( uh_max_loc[0][0] )
#                 uh_x = int( uh_max_loc[1][0] )  
                
#                 # Report results back to terminal
#                 print( '\tAdjusted UHmax: {} \tUHx: {} \tUHy: {}'.format( uh_max, xh[uh_x], yh[uh_y] ) )
                
#                 #-----------------------------------------------------------------------------
#                 # Begin: Correct BB_Box Coordinates Centered on Storm
#                 #-----------------------------------------------------------------------------
                
#                 x1 = int( uh_x - ( ( 0.50 * bb_box_dist * 1000.0 ) / dx ) )
#                 x2 = int( uh_x + ( ( 1.50 * bb_box_dist * 1000.0 ) / dx ) )
#                 y1 = int( uh_y - ( ( 0.75 * bb_box_dist * 1000.0 ) / dy ) )
#                 y2 = int( uh_y + ( ( 1.25 * bb_box_dist * 1000.0 ) / dy ) )
                
#                 # Is the new indices within a reasonable distance of the old ones?
#                 if( x1 < nx and x2 < nx and y1 < ny and y2 < ny ):
                
#                     print( '\nmod_x1: {} \tmod_x2: {}'.format( x1, x2 ) )
#                     print( 'mod_y1: {} \tmod_y2: {}'.format( y1, y2) )
                
#                 else:
#                     # Report to terminal
#                     print( '\n Attempting to move out of domain auto-correction being applied...' )
#                     x1 = x1
#                     x2 = nx - 1
#                     y1 = y1
#                     y2 = ny - 1
                
#                 #-----------------------------------------------------------------------------
#                 # End: Correct BB_Box Coordinates Centered on Storm
#                 #-----------------------------------------------------------------------------
                
#                 # Break the loop now that the UH tracking has been corrected
#                 overlap = 1
                
#             else:
#                 overlap = 1
            
#         #-----------------------------------------------------------------------------
        
#         x1_past = x1
#         x2_past = x2
#         y1_past = y1
#         y2_past = y2
    
#         #-----------------------------------------------------------------------------
#         # Begin: Generate Proximity Sounding Location
#         #-----------------------------------------------------------------------------
        
#         # Pull Sounding 30 km SE of UH max (prox_dist should be 1/2 of requested distance [a^2 + b^2 = c^2] )
#         sx_dist = ( ( prox_dist / 2.0 ) * 1000.0 ) / dx
#         sy_dist = ( ( prox_dist / 2.0 ) * 1000.0 ) / dx   
        
#         if( sx_dist < nx and sy_dist < ny ):
#             sx = uh_x + sx_dist
#             sy = uh_y - sy_dist
            
#             # Ensure indices are in integer format
#             sx = int(sx)
#             sy = int(sy)
            
#             # Save the indices
#             sxi = sx   
#             syj = sy 
            
#             # Get the grid-relative coordinates for proximity sounding (Required for clean plotting w/ Xarray)
#             sx = xh[ sx ].m
#             sy = yh[ sy ].m
        
#         # Get the grid-relative coordinates for UH max (Required for clean plotting w/ Xarray)
#         uhx_idx[k] = uh_x 
#         uhy_idx[k] = uh_y
#         uhx_grid[k] = xh[uh_x].m
#         uhy_grid[k] = yh[uh_y].m
        
#         uhx = xh[uh_x].m
#         uhy = yh[uh_y].m
        
#         #-----------------------------------------------------------------------------
#         # End: Generate Proximity Sounding Location
#         #-----------------------------------------------------------------------------
        
    
#         #-----------------------------------------------------------------------------
#         # Begin: Compute areal averaging about the UH maxima
#         #-----------------------------------------------------------------------------
        
#         # Define the search radius about the uh max as a callable integer value for looping purposes (Units: m)
#         uh_radius = int( ( updraft_radius * 1000.0 ) / dx )
        
#         # Initialize areal averaged values 
#         num_gridpoints = 0
#         uh_area = 0.0
#         w_500m = 0.0
#         w_1km = 0.0
#         w_3km = 0.0
#         w_5km = 0.0
#         zvort_surface = 0.0
#         zvort_500m = 0.0
#         zvort_1km = 0.0
#         zvort_3km = 0.0
#         zvort_5km = 0.0
        
#         # Get the looping indices 
#         i_start = uh_x - uh_radius
#         i_end = uh_x + uh_radius
#         j_start = uh_y - uh_radius
#         j_end = uh_y + uh_radius
        
#         # Begin: Loop through pre-defined updraft diameter in the x-direction
#         #-----------------------------------------------------------------------------
#         for i in range( i_start, i_end ):
        
#             # Begin: Loop through pre-defined updraft diameter in the y-direction
#             #-----------------------------------------------------------------------------
#             for j in range( j_start, j_end ):
                
#                 # If the radius lies within the domain
#                 if( i < nx and j < ny):
                    
#                     # Sum up the areal values of variables within the updraft
#                     num_gridpoints = num_gridpoints +  1 * (dx * dy)
#                     uh_area = uh_area + uh[ j, i ] * (dx * dy)
#                     w_500m = w_500m + w[ m500, j, i ] * (dx * dy)
#                     w_1km = w_1km + w[ km1, j, i ] * (dx * dy)
#                     w_3km = w_3km + w[ km3, j, i ] * (dx * dy)
#                     w_5km = w_5km + w[ km5, j, i ] * (dx * dy)
#                     zvort_surface = zvort_surface + zvort[ 0, j, i ] * (dx * dy)
#                     zvort_500m = zvort_500m + zvort[ m500, j, i ] * (dx * dy)
#                     zvort_1km = zvort_1km + zvort[ km1, j, i ] * (dx * dy)
#                     zvort_3km = zvort_3km + zvort[ km3, j, i ] * (dx * dy)
#                     zvort_5km = zvort_5km +  zvort[ km5, j, i ] * (dx * dy)
                
#                 # End the loop if the updraft moves outside of the domain (Not likely but needs to be explicitly declared to prevent code from breaking)    
#                 else:
#                     break
                
#             # End: Loop through pre-defined updraft diameter in the y-direction
#             #-----------------------------------------------------------------------------
        
#         # End: Loop through pre-defined updraft diameter in the x-direction
#         #-----------------------------------------------------------------------------
        
#         # Compute the areal average of updraft variables 
#         uh_area = uh_area / float( num_gridpoints )
#         w_500m = w_500m / float( num_gridpoints )
#         w_1km = w_1km / float( num_gridpoints )
#         w_3km = w_3km / float( num_gridpoints )
#         w_5km = w_5km / float( num_gridpoints )
#         zvort_surface = zvort_surface / float( num_gridpoints )
#         zvort_500m = zvort_500m / float( num_gridpoints )
#         zvort_1km = zvort_1km / float( num_gridpoints )
#         zvort_3km = zvort_3km / float( num_gridpoints )
#         zvort_5km = zvort_5km / float( num_gridpoints )    
        
#         #-----------------------------------------------------------------------------
#         # End: Compute areal averaging about the UH maxima
#         #-----------------------------------------------------------------------------

            
#         #-----------------------------------------------------------------------------
#         # Begin: Supercell Classification
#         #-----------------------------------------------------------------------------
        
#         # Check to see if areal UH exceeds classification threshold
#         if( uh_area >= uh_thresh ):
            
#             # Boolean classification (Supercell)
#             storm_class = 1
            
#             # Report to terminal
#             print( '\n Storm is Supercellular...' )
        
#         else:
            
#             # Boolean classification (Not Supercell)
#             storm_class = 0
            
#             # Report to terminal
#             print( '\n Storm is not Supercellular...' )
        
#         # Store classification into array for later usage
#         classification[k] = storm_class
        
#         #-----------------------------------------------------------------------------
#         # End: Supercell Classification
#         #-----------------------------------------------------------------------------



#     if( k >= track_on ):
        
#         # Get storm-motion vectors from intital contions
#         mean_uv, rm, lm = cm1_initial_conditions( ds, sxi, syj )
        
#         # Representative of the individual plots that compose the figure (Each letter is a separate subplot and . skips that position)
#         figure_mosaic = """
#                           AABB
#                           CCDD
#                         """ 
        
#         # Create figure and axes objects
#         fig, axes = plt.subplot_mosaic( 
#                                        mosaic = figure_mosaic, 
#                                        figsize = ( 20, 20 ), 
#                                        tight_layout = True
#                                       )
        
#         # Construct axis objects for the radar plot
#         ax_A = fig.add_subplot( axes['A'] )
        
#         # Run the Radar Plotter function
#         cm1_radar_plotter( 
#                           ds = ds, ax = ax_A, fig = fig, nx = nx, ny = ny,
#                           x1 = x1, x2 = x2 , y1 = y1, y2 = y2, 
#                           z = km1, sx = sx, sy = sy, uh_x = uhx, uh_y = uhy,
#                           skip_val = 20, shortname = shortname, current_model_time = current_model_time,
#                           terrain = True, peak_elv = 750.0, 
#                           cint = 75.0, id_terrain = True
#                          )
        
#         # Get the original specs for subplot D
#         ss = axes['B'].get_subplotspec()
        
#         # Remove the original instance of subplot D
#         axes['B'].remove()
        
#         # Re-construct the subplot using the original specs with the skewed x-axis projection 
#         axes['B'] = fig.add_subplot( ss, projection = 'skewx' )
        
#         # Run the Sounding Plotter function
#         cm1_sounding_plotter( 
#                               ds = ds, fig = fig, ax = axes['B'],
#                               sx = sxi, sy = syj, 
#                               shortname = shortname, current_model_time = current_model_time,
#                               knots = True, terrain = True, peak_rel = True
#                             )
        
#         # Construct axis objects for the radar plot
#         ax_C = fig.add_subplot( axes['C'] )

#         # (UH_X1 - 50 UH_X2 + 200 in original script)
#         cm1_cross_section_plotter( 
#                                   ds = ds, ds_interp = ds_interp, fig = fig, ax = ax_C,
#                                   x1 = x1 - 50, x2 = x2 + 100, z1 = 0, z2 = 43, uh_y = uh_y,
#                                   sx = sx, sy = sy, rm_x = rm[0], rm_y  = rm[1],
#                                   shortname = shortname, current_model_time = current_model_time, 
#                                   skip_val = 20, terrain = True, knots = True,
#                                   id_terrain = True, peak_pos = 350.0
#                                   )
        
#         # Construct axis objects for the radar plot
#         ax_D = fig.add_subplot( axes['D'] )
        
#         cm1_mesoanalysis_plotter( 
#                                   ds = ds, ax = ax_D, x1 = x1, x2 = x2, y1 = y1, y2 = y2, 
#                                   z1 = 0, z2 = 48, shortname = shortname , current_model_time =  current_model_time,
#                                   sx = sx, sy = sy, storm_motion = rm,
#                                   uh_x = uh_x, uh_y = uh_y, bb_box_dist = bb_box_dist, dx = dx, dy = dy,
#                                   terrain = True, id_terrain = True, knots = True,
#                                   km1 = 6, km3 = 13, km5 = 18, skip_val = 100,
#                                   nx = nx, ny = ny, nz = nz, smooth = 10,
#                                   peak_pos = 350.0, peak_elv = 750, cint = 75.0 
#                                 )

#         plt.show()

#         # Save the current RH plot
#         fig.savefig(
#                     fname = directory + '/plots/{}_4panel_{}_min.jpeg'.format( shortname, current_model_time ),
#                     dpi = 300,
#                     bbox_inches = "tight"
#                    )

# # Report the time required to run the function
# print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )

