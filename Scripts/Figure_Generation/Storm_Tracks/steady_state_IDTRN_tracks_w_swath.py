#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 13:39:37 2022

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
import xarray as xr
from metpy.units import units


from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

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
     zA500m, zA1km, zA3km, zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km                                                    
    ) = np.genfromtxt(
                      fname = filename,
                      delimiter = ',',
                      skip_header = 3,
                      usecols = np.arange( 0, 67, 1 ),
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
            zA500m, zA1km, zA3km, zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km                                                        
           )

# End Read-In Data Function
#-------------------------------------------------------------------


# Begin Moving Average Function
#-------------------------------------------------------------------

def guass_smoother( arr, smooth = 6 ) :
    ret_arr = scipynd.gaussian_filter1d( arr, smooth )
    return( ret_arr )

# End Moving Average Function
#-------------------------------------------------------------------


# Begin Domain Height Array Function
#-------------------------------------------------------------------

# Default values reflect CM1 model domain for this study
def witch_of_agnesi( nx = 2401, ny = 1601, dx = 250.0, 
                     h0 = 750.0, x0 = 350000.0, a = 50000.0 ):
    
    # Initialize arrays
    x = np.arange( 0, nx )                        # gridpoint array units: 0
    y = np.arange( 0, ny )                        # gridpoint array units: 0
    hx = np.zeros( shape = (nx) )                 # Elevation array units: m 
    
    # Compute x over the model domain
    for i in range( 0, nx ):
        x[i] = x[i] * dx
        
    # Compute y over the model domain    
    for i in range( 0, ny ):
        y[i] = y[i] * dx
        
    # Compute height over the model domain (Witch of Agnesi)
    for i in range( 0, len(hx) ):
        hx[i] = (h0) / ( 1 + ( (x[i]-x0) / a ) ** 2 )
            
    HX,HY = np.meshgrid( hx, y )
    
    # Return the height array back to user
    return( x, y, hx, HX, HY )

# End Domain Height Array Function
#-------------------------------------------------------------------



#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Begin Main Script
#-------------------------------------------------------------------
#-------------------------------------------------------------------

# Store program name
program = 'CM1_Analysis_IDTRN_Storm_Tracks_w_Swath.py'

nx = 2400
ny = 1600
uh_bot = 500
uh_top = 5500
uh_int = 1000
uh_alpha = 0.35
uh_sigma = 3

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( '\nBegin Program: {}'.format ( program ) )

x, y, hx, HX, HY = witch_of_agnesi()

# Convert to peak-relative
x = x - 350000.0

# Arrays to loop through simulations and assign color
sim = [ 'cs_ctl' , 'cs_trn', 'cs_mod', 'nc_ctl', 'nc_trn', 'nc_mod' ]
color = [ 'tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:cyan', 'tab:orange' ]

sup_start = [ 23, 25, 21, 28, 23, 26 ]
sup_end = [ 51, 52, 51, 67, 54, 48 ]
dis_end = [ 61, 68, 61, 72, 72, 71 ]
linear = [ 0, 0, 0, 71, 71, 71 ]

# Set tick size (Must be before calling plot object)
plt.rcdefaults()
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

# Create figure templates to add data to during loop
figure_mosaic = """
                A
                A
                B
                B
                C
                """
                
fig, axes =  plt.subplot_mosaic(
                                mosaic = figure_mosaic,
                                figsize = ( 12, 16 ),
                                tight_layout = True
                               )

# Create callable axes objects
ax_A = fig.add_subplot( axes['A'] )
ax_A.set_title( 'a: Idealized Terrain Steady-State Crossing Environment', fontsize = 20, fontweight = 'bold' )
ax_A.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_A.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_A.set_xlim( -250, 50 )
ax_A.set_ylim( 140, 325 )

ax_A.xaxis.set_ticklabels([])

terrain_fill = ax_A.contourf( 
                             x/1000.0, y/1000.0, HX,
                             levels = np.arange(0, 750.0 + 75.0, 75.0), 
                             cmap = mpl.cm.copper.reversed(),
                             alpha = 0.35
                            )


ax_B = fig.add_subplot( axes['B'] )
ax_B.set_title( 'b: Idealized Terrain Steady-State Non-Crossing Environment', fontsize = 20, fontweight = 'bold' )
ax_B.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_B.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_B.set_xlim( -250, 50 )
ax_B.set_ylim( 50, 235 )

ax_B.xaxis.set_ticklabels([])

ax_C = fig.add_subplot( axes['C'] )
ax_C.set_title( 'c: Idealized Terrain Cross-Section', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_C.set_xlabel( 'Distance from Peak (km)', fontsize = 18, fontweight = 'bold' )
ax_C.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )

ax_C.set_xlim( -250, 50 )
ax_C.set_ylim( 0, 1000 )

# Make background color sky blue
ax_C.set_facecolor( 'skyblue' )

# Plot the terrain profile
terrain_profile = ax_C.plot( x/1000.0, hx, linewidth = 2.5, fillstyle = 'full', color = 'k' )

# Fill in the area under the terrain curve
ax_C.fill_between( x/1000.0, hx, color = 'saddlebrown' )

# Add grid
ax_C.grid()

terrain_fill = ax_B.contourf( 
                             x/1000.0, y/1000.0, HX,
                             levels = np.arange(0, 750.0 + 75.0, 75.0), 
                             cmap = mpl.cm.copper.reversed(),
                             alpha = 0.35
                            )


# Begin: Loop through each simulation
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Construct filename string (!!!REQUIRES cleaned model output CSVs!!!)
    filename = '/Users/roger/Desktop/Steady_State_Track_Swath/' + sim[i] + '_model_output_stats.csv'
    
    # Read in current dataset
    ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km,
      zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
      uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
      mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
      shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
      srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
      lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
      zA500m, zA1km, zA3km, zA5km, zA8km,
      wzA500m, wzA1km, wzA3km, wzA5km, wzA8km  ) = read_stats_data( filename )
    
    uhx = uhx - 350.0 
    
    uhx_smooth = guass_smoother( uhx, 2 )
    uhy_smooth = guass_smoother( uhy, 2 )
    
    
    if( i <= 2 ):
        storm_track = ax_A.plot( uhx_smooth[ 20:dis_end[i] ], uhy_smooth[ 20:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
        # Demarcate where storm is supercellular
        ax_A.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
        ax_A.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    
        # Demarcate where storm is Linear
        if(linear[i]> 0):
            ax_A.scatter( uhx_smooth[ sup_end[i]+1  ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
            ax_A.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
    
    
        # Open .NC with swath data (!!!REQUIRES Final output file for UH Swath!!!)
        filename2 = '/Users/roger/Desktop/Steady_State_Track_Swath/' + sim[i] + '_cm1out_000073.nc'
        
        # Open the netCDF file with xarray 
        ds = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
      
        # Get the x-axis coordinates (km)
        xh = ds.coords[ 'xh' ].values * units( 'km' )
        xh = xh - 350.0 * units.km
        
        
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
                          levels = np.arange( uh_bot, uh_top, uh_int ), 
                          cmap = mpl.cm.gist_gray,
                          alpha = uh_alpha,
                          linewidths = 2.0
                         )
        
    else:
        storm_track = ax_B.plot( uhx_smooth[ 20:dis_end[i] ], uhy_smooth[ 20:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
        # Demarcate where storm is supercellular
        ax_B.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
        ax_B.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    
        # Demarcate where storm is Linear
        if(linear[i]> 0):
            ax_B.scatter( uhx_smooth[ sup_end[i]+1 ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s = 80, edgecolor = 'k', zorder = 10  )
            ax_B.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s = 80, edgecolor = 'k', zorder = 10  )
    
        # Open .NC with swath data (!!!REQUIRES Final output file for UH Swath!!!)
        filename2 = '/Users/roger/Desktop/Steady_State_Track_Swath/' + sim[i] + '_cm1out_000073.nc'
        
        # Open the netCDF file with xarray 
        ds = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
      
        # Get the x-axis coordinates (km)
        xh = ds.coords[ 'xh' ].values * units( 'km' )
        xh = xh - 350.0 * units.km
        
        
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
                          levels = np.arange( uh_bot, uh_top, uh_int ), 
                          cmap = mpl.cm.gist_gray,
                          alpha = uh_alpha,
                          linewidths = 2.0
                         )



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


# Include color bar legend
cbar  = plt.colorbar( 
                     terrain_fill,
                     location = 'right',
                     shrink = 0.9,
                     ticks = np.arange( 0, 750.0 + 75.0, 75.0 ),
                     drawedges = False,
                     ax = ax_A
                    ) 

# Include color bar legend
cbar2  = plt.colorbar( 
                     cs1,
                     location = 'right',
                     shrink = 0.9,
                     ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                     drawedges = True,
                     ax = ax_B
                    ) 

# Include color bar legend
cbar2  = plt.colorbar( 
                     cs1,
                     location = 'right',
                     shrink = 0.9,
                     ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                     drawedges = True,
                     ax = ax_C
                    ) 

# Set colorbar label 
cbar.set_label( label = 'Z (m)', fontsize = 16, fontweight = 'bold' )
cbar2.set_label( label = 'UH (m$^{2}$s$^{-2}$)', fontsize = 16, fontweight = 'bold' )

# Save figure
fig.savefig(
            fname = 'IDTRN_storm_tracks_w_swath.jpeg',
            dpi = 200,
            bbox_inches = "tight"
           )
