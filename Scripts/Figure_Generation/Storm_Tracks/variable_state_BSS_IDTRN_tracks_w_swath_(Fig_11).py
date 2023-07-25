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

from metpy.units import units
import metpy.calc as mpcalc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.plots as plots
import pint
import pint_xarray
import cf_xarray
from datetime import datetime


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

## Default values reflect CM1 model domain for this study
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
program = 'CM1_BSS_IDTRN_Analysis_Storm_Tracks_w_Swath.py'

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
sim = [ 'bss_cs_ctl', 'bss_cs_trn', 'bss_nc_ctl', 'bss_nc_trn' ]
color = [ 'tab:blue', 'tab:green', 'tab:purple', 'tab:cyan' ]

sup_start = [ 23, 28, 21, 23]
sup_end = [ 45, 45, 60, 46 ]
dis_end = [ 53, 53, 72, 61 ]
linear = [ 0, 0, 0, 0 ]

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
                .BBB.
                .CDE.
                .FFF.
                .FFF.
                .FFF.
                .GHI.
                .JJJ.
                """
                
fig, axes =  plt.subplot_mosaic(
                                mosaic = figure_mosaic,
                                figsize = ( 20, 40 ),
                                tight_layout = True
                               )

# Create callable axes objects
ax_A = fig.add_subplot( axes['A'] )
ax_A.set_title( 'a: Idealized Terrain Variable-State (BSS) Crossing Environment', fontsize = 20, fontweight = 'bold' )
ax_A.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_A.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_A.set_xlim( -225, 25 )
ax_A.set_ylim( 180, 300 )

terrain_fill = ax_A.contourf( 
                             x/1000.0, y/1000.0, HX,
                             levels = np.arange(0, 750.0 + 75.0, 75.0), 
                             cmap = mpl.cm.copper.reversed(),
                             alpha = 0.35
                            )


ax_F = fig.add_subplot( axes['F'] )
ax_F.set_title( 'a: Idealized Terrain Variable-State (BSS) Non-Crossing Environment', fontsize = 20, fontweight = 'bold' )
ax_F.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_F.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_F.set_xlim( -250, 50 )
ax_F.set_ylim( 110, 230 )

terrain_fill = ax_F.contourf( 
                             x/1000.0, y/1000.0, HX,
                             levels = np.arange(0, 750.0 + 75.0, 75.0), 
                             cmap = mpl.cm.copper.reversed(),
                             alpha = 0.35
                            )


ax_B = fig.add_subplot( axes['B'] )
ax_B.set_title( 'b: Idealized Terrain Cross-Section', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_B.set_xlabel( 'Distance from Peak (km)', fontsize = 18, fontweight = 'bold' )
ax_B.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )

ax_B.set_xlim( -225, 25 )
ax_B.set_ylim( 0, 1000 )

# Make background color sky blue
ax_B.set_facecolor( 'skyblue' )

# Plot the terrain profile
terrain_profile = ax_B.plot( x/1000.0, hx, linewidth = 2.5, fillstyle = 'full', color = 'k' )

# Fill in the area under the terrain curve
ax_B.fill_between( x/1000.0, hx, color = 'saddlebrown' )

# Add grid
ax_B.grid()


ax_J = fig.add_subplot( axes['J'] )
ax_J.set_title( 'b: Idealized Terrain Cross-Section', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_J.set_xlabel( 'Distance from Peak (km)', fontsize = 18, fontweight = 'bold' )
ax_J.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )

ax_J.set_xlim( -250, 50 )
ax_J.set_ylim( 0, 1000 )

# Make background color sky blue
ax_J.set_facecolor( 'skyblue' )

# Plot the terrain profile
terrain_profile = ax_J.plot( x/1000.0, hx, linewidth = 2.5, fillstyle = 'full', color = 'k' )

# Fill in the area under the terrain curve
ax_J.fill_between( x/1000.0, hx, color = 'saddlebrown' )

# Add grid
ax_J.grid()


# Begin: Loop through each simulation
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Construct filename string
    filename = '/scratch/rriggin/bss_cm1r20.3_fixed/Stat_Spreadsheets/' + sim[i] + '_model_output_stats.csv'
               
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
    # mode = guass_smoother( mode, 2.5 )
    
    
    if( i <= 1 ):
        storm_track = ax_A.plot( uhx_smooth[ 21:dis_end[i] ], uhy_smooth[ 21:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
        # Demarcate where storm is supercellular
        ax_A.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
        ax_A.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    
        # Demarcate where storm is Linear
        if(linear[i]> 0):
            ax_A.scatter( uhx_smooth[ sup_end[i]+1  ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
            ax_A.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
        
        ax_A.text( uhx_smooth[21] -10, uhy_smooth[21] +2, s = 'BSS0', fontsize = 18, fontweight = 'bold', color = 'k' )
        ax_A.text( uhx_smooth[36], uhy_smooth[36]+5, s = 'BSS1', fontsize = 18, fontweight = 'bold' , color = 'k' )
        ax_A.text( uhx_smooth[ dis_end[i] ]-5, uhy_smooth[ dis_end[i] ], s = 'BSS1->BSS2', fontsize = 18, fontweight = 'bold', color = 'k' )
        
    
        # Open .NC with swath data (!!! Requires appropriate directory containing BSS IDTRN final output files !!!)
        filename2 = '/scratch/rriggin/bss_cm1r20.3_fixed/' + sim[i].upper() + '/cm1out_000073.nc'
        
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
        storm_track = ax_F.plot( uhx_smooth[ 21:dis_end[i] ], uhy_smooth[ 21:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
        # Demarcate where storm is supercellular
        ax_F.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
        ax_F.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    
        # Demarcate where storm is Linear
        if(linear[i]> 0):
            ax_F.scatter( uhx_smooth[ sup_end[i]+1 ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s = 80, edgecolor = 'k', zorder = 10  )
            ax_F.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s = 80, edgecolor = 'k', zorder = 10  )
    
        if( i == 2 ):
            ax_F.text( uhx_smooth[ 20 ]-10, uhy_smooth[ 20 ]-5, s = 'BSS0', fontsize = 18, fontweight = 'bold', color = 'k' )
            ax_F.text( uhx_smooth[36], uhy_smooth[36]-5, s = 'BSS1', fontsize = 18, fontweight = 'bold', color = 'k' )
            ax_F.text( uhx_smooth[ dis_end[i] ], uhy_smooth[ dis_end[i] ] -5 , s = 'BSS1->BSS2', fontsize = 18, fontweight = 'bold', color = 'k' )
            
        if( i == 3 ):
            ax_F.text( uhx_smooth[ 20 ]-10, uhy_smooth[ 20 ], s = 'BSS0', fontsize = 18, fontweight = 'bold', color = 'k' )
            ax_F.text( uhx_smooth[36], uhy_smooth[36]+2, s = 'BSS0->BSS1', fontsize = 18, fontweight = 'bold', color = 'k' )
            ax_F.text( uhx_smooth[ dis_end[i] ] -2, uhy_smooth[ dis_end[i] ] +2, s = 'BSS1->BSS2', fontsize = 18, fontweight = 'bold', color = 'k' )
       
        
    
        # Open .NC with swath data (!!! Requires appropriate directory containing BSS IDTRN final output files !!!)
        filename2 = '/scratch/rriggin/bss_cm1r20.3_fixed/' + sim[i].upper() + '/cm1out_000073.nc'
        
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
        cs1 = ax_F.contourf( 
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

legend_2 = ax_F.legend( title = 'Simulations', loc = 'upper right' )
ax_F.add_artist( legend_2 )

ax_A.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right' )
ax_F.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right' )


ax_A.grid()
ax_F.grid()


# Include color bar legend
cbar  = plt.colorbar( 
                     terrain_fill,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.001,
                     ticks = np.arange( 0, 900, 150 ),
                     drawedges = False,
                     ax = ax_A
                    ) 

# Include color bar legend
cbar2  = plt.colorbar( 
                     cs1,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.075,
                     ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                     drawedges = True,
                     ax = ax_A
                    ) 

cbar3  = plt.colorbar( 
                     terrain_fill,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.001,
                     ticks = np.arange( 0, 900, 150 ),
                     drawedges = False,
                     ax = ax_F
                    ) 

# Include color bar legend
cbar4  = plt.colorbar( 
                     cs1,
                     location = 'bottom',
                     shrink = 0.8,
                     pad = 0.075,
                     ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                     drawedges = True,
                     ax = ax_F
                    ) 

# Set colorbar label 
cbar.set_label( label = 'Z (m)', fontsize = 16, fontweight = 'bold', x = 1.1, y = 2 )
cbar2.set_label( label = 'UH (m$^{2}$s$^{-2}$)', fontsize = 16, fontweight = 'bold', x = 1.1 , y = 2)
cbar3.set_label( label = 'Z (m)', fontsize = 16, fontweight = 'bold', x = 1.1, y = 2 )
cbar4.set_label( label = 'UH (m$^{2}$s$^{-2}$)', fontsize = 16, fontweight = 'bold', x = 1.1, y = 2 )



#-----------------------------------------------------------------------------
# BSS Sounding Section
#-----------------------------------------------------------------------------

# (!!! Requires appropriate directory containing Initial Conditions Model Output Files !!!)
output_dir = '/scratch/rriggin/bss_cm1r20.3_fixed/IC_Output_Files/'
sim_type = [ 'Crosser', 'Non-Crosser' ]
sound_type = [ 'Upstream', 'Peak', 'Downstream', 'Upstream', 'Peak', 'Downstream' ]
bss_time = [ 't = BSS0', 't = BSS1', 't = BSS2', 't = BSS0', 't = BSS1', 't = BSS2' ]

# List for looping through each axis
current_ax = [ 'C', 'D', 'E', 'G', 'H', 'I' ]
current_ax_label = [ 'C', 'D', 'E', 'C', 'D', 'E' ]

#------------------------------------------------------------------------------
# Begin: Loop to plot on each axis
#------------------------------------------------------------------------------
for i in range( 0, len( current_ax ) ):
    
    # Get the original specs for subplot D
    ss = axes[ current_ax[i] ].get_subplotspec()
    
    # Remove the original instance of subplot D
    axes[ current_ax[i] ].remove()
    
    # Re-construct the subplot using the original specs with the skewed x-axis projection 
    axes[ current_ax[i] ] = fig.add_subplot( ss, projection = 'skewx' )
    
    # Call the current axes object
    ax = fig.add_subplot( axes[ current_ax[i] ] )
    
    # Create a SkewT plotting object
    skew = plots.SkewT( fig, rotation = 45, subplot = ax, aspect = 100 )
    
    # Custom axes limits
    skew.ax.set_ylim( 1000, 100 )
    skew.ax.set_xlim( -65, 35 ) 
    
    # Add fiducial lines
    skew.plot_dry_adiabats( alpha = 0.35 )
    skew.plot_moist_adiabats( alpha = 0.35 )
    skew.plot_mixing_lines( alpha = 0.35 )
    skew.ax.axvline( 0, color = 'c', linestyle = '--', linewidth = 2.5, alpha = 0.35 )
    
    # Set the x-axis label
    ax.set_xlabel( 
                  xlabel = 'Temperature (\u00B0C)',
                  fontsize = 14, fontweight = 'bold'
                 )
    
    # Set the y-axis label
    if( i == 0 or i == 3 ):
        ax.set_ylabel( 
                      ylabel = 'Pressure (hPa)',
                      fontsize = 14, fontweight = 'bold'
                     )
        
    # Label every other tickmark    
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    
    
    # Place a hodograph into the current axis
    hax = inset_axes( skew.ax, '35%', '35%', loc= "upper left" )
    
    # # Label the hodograph
    # hax.set_xlabel( 
    #                'Hodograph (ms$^{-1}$)',
    #                 loc = 'center',
    #                 fontsize = 12, fontweight = 'bold'
    #                )

    # Logic for Axis Title
    if( i < 3 ):
        ax.set_title( 
                     current_ax_label[i] + ": " + sound_type[i], loc = 'left',
                     fontsize = 18, fontweight = 'bold'
                    )
        ax.set_title( 
                     bss_time[i], loc = 'right',
                     fontsize = 16, fontweight = 'bold'
                    )
    else:
        ax.set_title( 
                     current_ax_label[i] + ": " + sound_type[i], loc = 'left',
                     fontsize = 18, fontweight = 'bold'
                    )
        ax.set_title( 
                     bss_time[i], loc = 'right',
                     fontsize = 16, fontweight = 'bold'
                    )

    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Pull sounding from center of domain
    sx = 300
    sy = 200
    
    # Open the current file
    if( i < 3 ):
        fname = output_dir + 'cm1out_' + sim_type[0] + '_' + sound_type[i] + '.nc'
    else:
        fname = output_dir + 'cm1out_' + sim_type[1] + '_' + sound_type[i] + '.nc'
        
    # Report to terminal
    print( '\tNow opening {}'.format( fname ) )
        
    # Open the current netCDF file with xarray 
    ds = xr.open_dataset( fname, engine = "netcdf4", decode_cf = True ) 
   
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the u-component of the wind (m/s)
    u = ds.metpy.parse_cf( 'uinterp' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get the v-component of the wind (m/s)
    v = ds.metpy.parse_cf( 'vinterp' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get potential temperature (K)
    th = ds.metpy.parse_cf( 'th' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get pressure (Pa)
    prs = ds.metpy.parse_cf( 'prs' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get water-vapor mixing ratio (kg/kg)
    qv = ds.metpy.parse_cf( 'qv' ).isel( time = 0, yh = sy, xh = sx ) 
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get heights (km)
    zh = ds.coords[ 'zh' ].values * units( 'km' )
    
    # Report to terminal
    print( '\tNow closing {}'.format( fname ) )
    
    # Close the netCDF file
    ds.close()
    
    #-----------------------------------------------------------------------------
    # End: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Unit conversions required for MetPy functions
    #   Notes: metpy.quantify() adds native units from Dataset
    #   Notes: pint_xarray provides wrapper to convert units using pint.to()
    #-----------------------------------------------------------------------------
       
    # Get native units for required variables
    u = u.metpy.quantify()
    v = v.metpy.quantify()
    prs = prs.metpy.quantify()
    qv = qv.metpy.quantify()
    
    # Convert to required units for Air Temp & Td calculations
    prs = prs.pint.to( 'hPa' )  
    qv = qv.pint.to( 'g/kg' )
    
    #-----------------------------------------------------------------------------
    # End: Unit conversions required for MetPy functions
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Air Temperature & Dew Point Calculations
    #-----------------------------------------------------------------------------
    
    # Convert potential temperature to air temperature (K)
    temp = mpcalc.temperature_from_potential_temperature( prs, th ) 
    
    # Get native units
    temp = temp.metpy.quantify()
    
    # Convert from K to degC
    temp = temp.pint.to( 'degC' )
    
    # Compute the water vapor pressure using pressure and mixing-ratio
    e  = mpcalc.vapor_pressure( prs, qv )
    
    # Compute dew point from vapor pressure (degC)
    td = mpcalc.dewpoint( e ) 
    
    #-----------------------------------------------------------------------------
    # End: Air Temperature & Dew Point Calculations
    #-----------------------------------------------------------------------------

    # Plot temp and dewpoint
    skew.plot( prs, temp, 'tab:red', linewidth = 2.5, label = 'Temp ' )
    skew.plot( prs, td, 'tab:green', linewidth = 2.5, label = 'Dewpoint' )
    
    # Calculate full parcel profile
    prof = mpcalc.parcel_profile( prs, temp[0], td[0] ).pint.to('degC')

    # Calculate LCL height
    lcl_pressure, lcl_temperature = mpcalc.lcl( prs[0], temp[0], td[0] )
    
    # Calculate LFC height
    lfc_pressure, lfc_temperature = mpcalc.lfc( prs, temp, td )
    
    # Calculate Equilibrium Level
    el_pressure, el_temperature = mpcalc.el( prs, temp, td )
    
    # Plot parcel profile
    skew.plot( prs, prof, linewidth = 2.5, linestyle = '--', color = 'black', alpha = 0.65, label = 'Parcel Path' )
    
    # Plot LCL, LFC, and EL levels
    # skew.plot( lcl_pressure, lcl_temperature, marker = '_', markersize = 10, color = 'b' )
    # skew.plot( lfc_pressure, lfc_temperature, marker = '_', markersize = 10, color = 'r' )
    # skew.plot( el_pressure, el_temperature, marker = '_', markersize = 10, color = 'k' )
    
    # Label LCL, LFC, and EL levels
    # plt.text( lcl_temperature.to( 'degC').m + 5, lcl_pressure.m + 25 , 'LCL', color = 'b', size = 18 )
    # plt.text( lfc_temperature.to( 'degC').m + 5, lfc_pressure.m + 20 , 'LFC', color = 'r', size = 18 )
    # plt.text( el_temperature.to( 'degC').m + 5, el_pressure.m + 5, 'EL', color = 'k', size = 18 )
    
    # Shade areas of CAPE and CIN
    skew.shade_cin( prs, temp, prof, alpha = 0.2, label = 'SBCIN' )
    skew.shade_cape( prs, temp, prof, alpha = 0.2, label = 'SBCAPE' )

    # Create the hodograph object
    h = plots.Hodograph( hax, component_range = 45.0 )
    hax.yaxis.tick_right()
   
    # Add a polar grid
    h.add_grid( increment = 10.0 )

    # Height (km AGL) intervals used for hodograph plotting
    intervals = np.array( [ 0.0, 3.0, 6.0, 9.0, 12.0, 15.0 ] ) * units( 'km' )
    
    # Height (km AGL) Mask
    zh_mask = zh <= 15.0 * units( 'km' )

    # Colormap scheme to match interval bins for hodograph
    cmap = [ 'tab:red', 'tab:green', 'tab:olive', 'tab:cyan', 'tab:purple' ]

    # Plot the hodograph
    hodograph = h.plot_colormapped( u[zh_mask], v[zh_mask], c = zh[ zh_mask ], intervals = intervals, colors = cmap )
    
    # Modify hodograph if requested
    #----------------------------------------------
    # if( mod_hodo == True ):
        
    #     mod_angle = units.Quantity( 25, 'degrees' )   
        
    #     # Get magnitude and direction from original components
    #     wspd = mpcalc.wind_speed( u, v )
    #     wdir = mpcalc.wind_direction( u , v )
        
    #     # Modify the direction of the hodograph evenly throughout all layers
    #     for j in range( 0, len( wdir) ):
    #         wdir[j] = wdir[j] + mod_angle
        
    #     # Compute the new wind components
    #     u2,v2 = mpcalc.wind_components( wspd, wdir )
        
    # # Plot the hodograph (m/s) with a color map based on height
    # hodo2 = h.plot( 
    #               u = u2,
    #               v = v2,
    #               # color = hodo_col[i],
    #               # linestyle = line[0],
    #               # alpha = 0.45,
    #              )
    
    # Set hodograph plot limits
    hax.set_xlim( -5, 45 )
    hax.set_ylim( -5, 35 )    
    
    # End Modify Hodograph
    #---------------------------------------------

    # # Compute MLCAPE/CIN
    # mlcape, mlcin = mpcalc.mixed_layer_cape_cin( prs, temp, td )
    
    # Compute Bunkers Storm Motions and Mean Wind
    rm, lm, mean_uv = mpcalc.bunkers_storm_motion(
                                                  pressure = prs,
                                                  u = u,
                                                  v = v,
                                                  height = zh
                                                 )

    # Plot the RM-motion as a vector
    h.wind_vectors( u = rm[0], v = rm[1], color = 'dimgrey', scale = 1, alpha = 0.75, width = 1, headwidth = 3 )
    
    # # Plot the Bunker's Right & Left Motions as labeled pointd
    # txt5 = plt.text( rm[0].m, rm[1].m, 'RM', size = 8 )
    # # txt6 = plt.text( lm[0].m, lm[1].m, 'LM', size = 8 )
    
    # # Place circles around RM & LM labels
    # plt.scatter( rm[0].m+2.5, rm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
    # # plt.scatter( lm[0].m+2.5, lm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
    

    # # Compute 0-3 km SRH
    # srh_3km = mpcalc.storm_relative_helicity(
    #                                          u = u,
    #                                          v = v,
    #                                          height = zh,
    #                                          depth = units.Quantity(3000, 'm'),
    #                                          storm_u = rm[0],
    #                                          storm_v = rm[1]
    #                                         )
    
    # srh_1km = mpcalc.storm_relative_helicity(
    #                                          u = u,
    #                                          v = v,
    #                                          height = zh,
    #                                          depth = units.Quantity(1000, 'm'),
    #                                          storm_u = rm[0],
    #                                          storm_v = rm[1]
    #                                         )

    # p2 = np.arange( -45, 20, 6)
    
    # txt1 = plt.text( 0.0, p2[3], 'MLCAPE: ' + str( round( float(mlcape.m), 2 ) ) + ' Jkg$^{-1}$', fontsize = 14 )
    # txt2 = plt.text( 0.0, p2[2], 'MLCIN: '  + str( round( float(mlcin.m),  2 ) ) + ' Jkg$^{-1}$', fontsize = 14 )   
    # txt3 = plt.text( 0.0, p2[1], '0-3 km SRH: ' + str( round( float(srh_3km[0].m), 2 ) ) + ' m$^2$s$^{-2}$', fontsize = 14 )
    # txt4 = plt.text( 0.0, p2[0], '0-1 km SRH: ' + str( round( float(srh_1km[0].m), 2 ) ) + ' m$^2$s$^{-2}$', fontsize = 14 )
    
    # Create a mask for wind barb plotting
    mask = prs >= 100.0 * units.hPa
    
    # Determine which units of wind should be plotted
    skew.plot_barbs( prs[ mask ][::2], u[ mask ][::2], v[ mask ][::2] )     

    # Add legend 
    if( i == 0 or i == 3 ):
        skew.ax.legend( loc = 3, facecolor = 'grey',  framealpha = 0.2, fontsize = 14 )
    
#------------------------------------------------------------------------------
# End: Loop to plot on each axis
#------------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# End: BSS Sounding Section
#-----------------------------------------------------------------------------

# Save figure
fig.savefig(
            fname = 'BSS_IDTRN_storm_tracks_w_swath.jpeg',
            dpi = 300,
            bbox_inches = "tight"
           )
