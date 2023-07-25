#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:56:14 2022

@author: roger
"""

# Import external modules
#-------------------------------------------------------------------
import sys
import os
import numpy as np
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as scipynd
import matplotlib.lines as mlines
import pint_xarray
import cf_xarray
import pint_xarray
import xarray as xr
from metpy.units import units


# Begin Read-In Data Function
#-------------------------------------------------------------------

def read_stats_data( filename ):
    
    # Report program status to the terminal
    print( '\n\t Now opening input file: {} '.format( filename ) )
    
    # Store data columns
    ( 
     mode, time, x1, y1, x2, y2, zs, 
     w500m, w1km, w3km, w5km,
     zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
     uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
     mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
     shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
     srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
     lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
     zA500m, zA1km, zA3km, zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth,                                
     cp_intensity, cp_area                    
    ) = np.genfromtxt(
                      fname = filename,
                      delimiter = ',',
                      skip_header = 3,
                      usecols = np.arange( 0, 70, 1 ),
                      unpack = True
                    )
        
    # Report program status to terminal
    print( '\n\t Raw data successfully stored into dataframe structure...')

    # Return read in data arrays
    return ( 
            mode, time, x1, y1, x2, y2, zs, 
            w500m, w1km, w3km, w5km,
            zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
            uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
            mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
            shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
            srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
            lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
            zA500m, zA1km, zA3km, zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth,
            cp_intensity, cp_area                                                        
           )

# End Read-In Data Function
#-------------------------------------------------------------------


# Begin Moving Average Function
#-------------------------------------------------------------------

def gauss_smoother( arr, smooth = 2 ) :
    ret_arr = scipynd.gaussian_filter1d( arr, smooth )
    return( ret_arr )

# End Moving Average Function
#-------------------------------------------------------------------


#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Begin Main Script
#-------------------------------------------------------------------
#-------------------------------------------------------------------

# Store program name
program = 'CM1_Analysis_RLTRN_Storm_Tracks_w_Swath.py'

nx = 2400
ny = 2000
uh_bot = 500
uh_top = 5500
uh_int = 1000
uh_alpha = 0.35
uh_sigma = 3

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( '\nBegin Program: {}'.format ( program ) )

# Arrays to loop through simulations and assign color
sim = [ 'cs_rltrn_out', 'cs_rltrn_no_out', 'cs_rltrn_all', 'nc_rltrn_all', 'nc_rltrn_no_out' ]
color = [ 'tab:red', 'tab:green' , 'tab:blue', 'tab:purple', 'tab:cyan' ]

# Arrays for demarcating storm mode and track
sup_start = [ 44, 43, 43, 46, 47 ]
sup_end = [ 65, 53, 66, 65, 82 ]
dis_end = [ 82, 64, 77, 97, 97 ]
linear = [ 0, 0, 0, 0, 96 ]

# Adjust all timestamps by 2 hrs to acccount for delayed CI
for x in range (0, len( sup_start ) ):
    sup_start[x] = sup_start[x] - 24
    sup_end[x] = sup_end[x] - 24
    dis_end[x] = dis_end[x] - 24
    if( linear[x] > 0 ):
        linear[x] = linear[x] - 24
        
# Open terrain file for background data (!!! REQUIRES output file with interpolated terrain !!!)
terr_file = '/Users/roger/Desktop/Realistic_Terrain_Track_Swath/cs_rltrn_no_out_cm1out_000097.nc' 
ds = xr.open_dataset( terr_file, engine = "netcdf4", decode_cf = True )
xh = ds.coords[ 'xh' ].values * units( 'km' )
yh = ds.coords[ 'yh' ].values * units( 'km' )
zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0 )

# Create a plotting grid based on the x,y coordinates
X,Y = np.meshgrid( xh, yh )

# Set tick size (Must be before calling plot object)
plt.rcdefaults()
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

# Create figure templates to add data to during loop
figure_mosaic = """
                .AAA.
                .AAA.
                .AAA.
                .CCC.
                .BBB.
                .BBB.
                .BBB.
                .DDD.
                """
                
fig, axes =  plt.subplot_mosaic(
                                mosaic = figure_mosaic,
                                figsize = ( 15, 30 ),
                                tight_layout = True
                               )

# Create callable axes objects
ax_A = fig.add_subplot( axes['A'] )
ax_A.set_title( 'a: Realistic Terrain: Steady-State Crossing Environment', fontsize = 20, fontweight = 'bold' )
ax_A.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_A.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_A.set_xlim( 115, 495 )
ax_A.set_ylim( 125, 345 )

ax_A.xaxis.set_ticklabels([])

# Terrain Contours
terrain_contour = ax_A.contourf( 
                                X, Y, zs,
                                levels = np.arange( 0, 2250, 250), 
                                cmap = plt.cm.get_cmap( "copper" ).reversed(),
                                alpha = 0.35
                               )

ax_C = fig.add_subplot( axes['C'] )
ax_C.set_xlim( 115, 495 )
ax_C.set_ylim( -50, 1250 )

ax_C.set_title( 'b: Realistic Terrain Cross-Sections', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_C.set_xlabel( 'Zonal Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_C.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )

# Create callable axes objects
ax_B = fig.add_subplot( axes['B'] )
ax_B.set_title( 'a: Realistic Terrain Steady-State Non-Crossing Environment', fontsize = 20, fontweight = 'bold' )
ax_B.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_B.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_B.xaxis.set_ticklabels([])

ax_B.set_xlim( 115, 495 )
ax_B.set_ylim( 125, 345 )

ax_D = fig.add_subplot( axes['D'] )
ax_D.set_xlim( 115, 495 )
ax_D.set_ylim( -50, 1250 )

ax_D.set_title( 'b: Realistic Terrain Cross-Sections', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_D.set_xlabel( 'Zonal Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_D.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )

# Terrain Contours
terrain_contour = ax_B.contourf( 
                                X, Y, zs,
                                levels = np.arange( 0, 2250, 250), 
                                cmap = plt.cm.get_cmap( "copper" ).reversed(),
                                alpha = 0.35
                               )

# Begin: Loop through each simulation
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Construct filename string (!!! REQUIRES Appropriate link to directory containing RLTRN CSVs !!!)
    filename = '/Users/roger/Desktop/Realistic_Terrain_Track_Swath/' + sim[i] + '_output_stats.csv'
    
    # Read in current dataset
    ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km,
      zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
      uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
      mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
      shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
      srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
      lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
      zA500m, zA1km, zA3km, zA5km, zA8km,
      wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth,
      cp_intensity, cp_area ) = read_stats_data( filename )
    
    uhx_smooth = gauss_smoother( uhx, 1. )
    uhy_smooth = gauss_smoother( uhy, 1. )
    
    
    if( i <= 2 ):
        storm_track = ax_A.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], uhy_smooth[ sup_start[i]:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
        # Demarcate where storm is supercellular
        ax_A.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
        ax_A.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    
        # Demarcate where storm is Linear
        if(linear[i]> 0):
            ax_A.scatter( uhx_smooth[ sup_end[i]+1  ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
            ax_A.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
    
    
        # Open .NC with swath data (!!! REQUIRES Appropriate link to directory containing RLTRN Final Output files !!!)
        filename2 = '/Users/roger/Desktop/Realistic_Terrain_Track_Swath/' + sim[i] + '_cm1out_000097.nc'
        
        # Open the netCDF file with xarray 
        ds = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
      
        # Get the x-axis coordinates (km)
        xh = ds.coords[ 'xh' ].values * units( 'km' )
       
        # Get the y-axis coordinates (km)      
        yh = ds.coords[ 'yh' ].values * units( 'km' )
        
        # Create a plotting grid based on the x,y coordinates
        X,Y = np.meshgrid( xh, yh )
        
        # Get the 2-5 km Updraft Helicity Swath (m^2/s^2)
        shs = ds.metpy.parse_cf( 'shs' ).isel( time = 0, yh = slice( 0, ny), xh = slice( 0, nx )  )
    
        # UH Swath Contours         
        cs1 = ax_A.contourf( 
                          X, Y,
                          scipynd.gaussian_filter( shs[ :, : ], uh_sigma ),
                          levels = np.arange( uh_bot, uh_top+uh_int, uh_int ), 
                          alpha = uh_alpha,
                          extend = 'max',
                          cmap = mpl.cm.gist_gray,
                          linewidths = 2.0
                         )
        
        # Run terrain profile through guassian smoother    
        zs_alt = gauss_smoother( zs * 1000, 1. )
        
        # Same plotting methods for panel G
        ax_C.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color ='k', linewidth = 1.5, linestyle = '-', label = sim[i].upper() + ' Elevation' )
        
        # Fill in the area under the terrain curve
        ax_C.fill_between( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color = color[i], alpha = 0.75 )
        
        # Demarcate start/end supercell mode
        ax_C.scatter( uhx_smooth[ sup_start[i] ], zs_alt[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        ax_C.scatter( uhx_smooth[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

        # Demarcate start/end linear mode
        if( linear[i] > 0 ):
            ax_C.scatter( uhx_smooth[ sup_end[i] + 1 ], zs_alt[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
            ax_C.scatter( uhx_smooth[ linear[i] ], zs_alt[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

        
        
    else:
        storm_track = ax_B.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], uhy_smooth[ sup_start[i]:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
        # Demarcate where storm is supercellular
        ax_B.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
        ax_B.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    
        # Demarcate where storm is Linear
        if(linear[i]> 0):
            ax_B.scatter( uhx_smooth[ sup_end[i]+1 ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s = 80, edgecolor = 'k', zorder = 10  )
            ax_B.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s = 80, edgecolor = 'k', zorder = 10  )
    
        # Open .NC with swath data (!!! REQUIRES Appropriate link to directory containing RLTRN Final Output files !!!)
        filename2 = '/Users/roger/Desktop/Realistic_Terrain_Track_Swath/' + sim[i] + '_cm1out_000097.nc'
        
        # Open the netCDF file with xarray 
        ds = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
      
        # Get the x-axis coordinates (km)
        xh = ds.coords[ 'xh' ].values * units( 'km' )
        
        # Get the y-axis coordinates (km)      
        yh = ds.coords[ 'yh' ].values * units( 'km' )
        
        # Create a plotting grid based on the x,y coordinates
        X,Y = np.meshgrid( xh, yh )
        
        # Get the 2-5 km Updraft Helicity Swath (m^2/s^2)
        shs = ds.metpy.parse_cf( 'shs' ).isel( time = 0, yh = slice( 0, ny), xh = slice( 0, nx )  )
    
        # UH Swath Contours         
        cs1 = ax_B.contourf( 
                          X, Y,
                          scipynd.gaussian_filter( shs[ :, : ], uh_sigma ),
                          levels = np.arange( uh_bot, uh_top+uh_int, uh_int ), 
                          alpha = uh_alpha,
                          cmap = mpl.cm.gist_gray,
                          linewidths = 2.0
                         )
        
        # Run terrain profile through guassian smoother    
        zs_alt = gauss_smoother( zs * 1000, 1. )
        
        # Same plotting methods for panel G
        ax_D.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color = 'k', linewidth = 1.5, linestyle = '-', label = sim[i].upper() + ' Elevation' )
        
        # Fill in the area under the terrain curve
        ax_D.fill_between( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color = color[i], alpha = 0.75 )
       
        # Demarcate start/end supercell mode
        ax_D.scatter( uhx_smooth[ sup_start[i] ], zs_alt[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        ax_D.scatter( uhx_smooth[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

        # Demarcate start/end linear mode
        if( linear[i] > 0 ):
            ax_D.scatter( uhx_smooth[ sup_end[i] + 1 ], zs_alt[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
            ax_D.scatter( uhx_smooth[ linear[i] ], zs_alt[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

        
sup_lab = mlines.Line2D( [], [], color = 'k', marker = 'v', label = 'Supercell' )
lin_lab = mlines.Line2D( [], [], color = 'k', marker = 's', label = 'Linear' )

legend_1 = ax_A.legend( title = 'Simulations', loc = 'upper right' )
ax_A.add_artist( legend_1 )

legend_2 = ax_B.legend( title = 'Simulations', loc = 'upper right' )
ax_B.add_artist( legend_2 )

ax_A.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right' )
ax_B.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right' )


ax_A.grid()
ax_B.grid()
ax_C.grid()
ax_D.grid()

# Include color bar legend
cbar  = plt.colorbar( 
                     terrain_contour,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.001,
                     ticks = np.arange( 0, 2250, 250 ),
                     drawedges = False,
                     ax = ax_A
                    ) 

# Include color bar legend
cbar1  = plt.colorbar( 
                     terrain_contour,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.001,
                     ticks = np.arange( 0, 2250, 250 ),
                     drawedges = False,
                     ax = ax_B
                    ) 

# Include color bar legend
cbar2  = plt.colorbar( 
                     cs1,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.05,
                     ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                     drawedges = True,
                     ax = ax_A
                    ) 

cbar3  = plt.colorbar( 
                     cs1,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.05,
                     ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                     drawedges = True,
                     ax = ax_B
                    ) 

# Set colorbar label 
cbar.set_label( label = 'Z (m)', fontsize = 16, fontweight = 'bold' )
cbar2.set_label( label = 'UH (m$^{2}$s$^{-2}$)', fontsize = 16, fontweight = 'bold' )
cbar1.set_label( label = 'Z (m)', fontsize = 16, fontweight = 'bold' )
cbar3.set_label( label = 'UH (m$^{2}$s$^{-2}$)', fontsize = 16, fontweight = 'bold' )


# Save figure
fig.savefig(
            fname = 'RLTRN_storm_track_w_swath.jpeg',
            dpi = 300,
            bbox_inches = "tight"
           )
