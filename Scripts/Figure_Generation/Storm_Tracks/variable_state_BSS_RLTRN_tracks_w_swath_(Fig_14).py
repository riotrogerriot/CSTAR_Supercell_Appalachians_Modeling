#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
import matplotlib.patheffects as PathEffects


#------------------------------------------------------------------------------
# Begin: User-defined Functions
#------------------------------------------------------------------------------

# CSV Read-In Function Specific to Tracking-Algorithm Output
#------------------------------------------------------------------------------
def read_stats_data( filename ):
    
    # Report program status to the terminal
    print( '\n\t Now opening input file: {} '.format( filename ) )
    
    # Store data columns
    ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km, zvort_surface, zvort_500m, 
      zvort_1km, zvort_3km, zvort_5km, uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
      mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir, shear1km, shear1km_dir, shear3km, shear3km_dir, 
      shear6km, shear6km_dir, srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, 
      cape3, lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km, zA500m, zA1km, zA3km, 
      zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth, cp_intensity, cp_area 
     ) = np.genfromtxt( fname = filename, delimiter = ',', skip_header = 3, usecols = np.arange( 0, 70, 1 ), unpack = True )
        
    # Report program status to terminal
    print( '\n\t Raw data successfully stored into dataframe structure...')

    # Return read in data arrays
    return ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km, zvort_surface, zvort_500m, 
      zvort_1km, zvort_3km, zvort_5km, uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
      mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir, shear1km, shear1km_dir, shear3km, shear3km_dir, 
      shear6km, shear6km_dir, srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, 
      cape3, lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km, zA500m, zA1km, zA3km, 
      zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth, cp_intensity, cp_area )


# Moving Average Function
#------------------------------------------------------------------------------
def gauss_smoother( arr, smooth = 2 ) :
    ret_arr = scipynd.gaussian_filter1d( arr, smooth )
    return( ret_arr )

#------------------------------------------------------------------------------
# End: User-defined Functions
#------------------------------------------------------------------------------


# Program name, main directory, and sub-directories required to run this script
program = 'CM1_Analysis_RLTRN_Storm_Tracks_w_Swath.py'
os_type = 0
if( os_type == 0 ):
    wdir = '/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR_Modeling_Project/Simulations/Variable_State_Sims/Stats_Spreadsheets/'
else:
    wdir = 'C:/Users/rriggin/OneDrive - University of North Carolina at Charlotte/CSTAR_Modeling_Project/Simulations/Variable_State_Sims/Stats_Spreadsheets/'
terr_file = wdir + 'BSS_RLTRN_Final_Output_Files/bss_cs_rltrn_cm1out_000093.nc' 

# User-defined parameters
#------------------------------------------------------------------------------
nx = 2400           # Num of x-gridpoints
ny = 2000           # Num of y-gridpoints
uh_bot = 500        # Min UH plot value (m2/s2)
uh_top = 4500       # Max UH plot value (m2/s2)
uh_int = 1000       # UH plotting interval (m2/s2)
uh_sigma = 5        # UH smoothing parameter
t_max = 1400        # Max terrain height (m)
t_int = 200         # Terrain plotting interval (m)
x1 = 125            # X-axis min (km)
x2 = 330            # X-axis max (km)
alpha = 0.45        # Contourf transparency parameter
#------------------------------------------------------------------------------


# Simulation-Specific Parameters
#------------------------------------------------------------------------------
sim = [ 'bss_nc_rltrn', 'bss_cs_rltrn' ]        # Sim names
color = [ 'tab:blue', 'tab:red'  ]              # Sim colors
sup_start = [ 43, 40 ]                          # Supercell Start Index
sup_end = [ 62, 55 ]                            # Supercell End Index
dis_end = [ 79, 58 ]                            # Dissipation Time Index
linear = [ 0, 0 ]                               # Upscale Growth Index
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Begin Main Script
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( '\nBegin Program: {}'.format ( program ) )

# Adjust all timestamps by 2 hrs to acccount for delayed CI w/ RLTRN
for x in range (0, len( sup_start ) ):
    sup_start[x] = sup_start[x] - 24
    sup_end[x] = sup_end[x] - 24
    dis_end[x] = dis_end[x] - 24
    if( linear[x] > 0 ):
        linear[x] = linear[x] - 24
        
# Open terrain file for background data
ds = xr.open_dataset( terr_file, engine = "netcdf4", decode_cf = True )
xh = ds.coords[ 'xh' ].values * units( 'km' )
yh = ds.coords[ 'yh' ].values * units( 'km' )
zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0 )

# Create a plotting grid based on the x,y coordinates
X,Y = np.meshgrid( xh, yh )

# Set tick size (Must be before calling plot object)
plt.rcdefaults()
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)

# Plot layout layout
figure_mosaic = """
                .AAA.
                .AAA.
                .BBB.
                """
                
# Create figure with axes based on the layout
fig, axes =  plt.subplot_mosaic( mosaic = figure_mosaic, figsize = ( 16, 10.5 ), tight_layout = True )

# Create callable Panel A objects, define label, and set plotting limits
ax_A = fig.add_subplot( axes['A'] )
ax_A.set_title( 'a: RLTRN Variable-State (BSS) Environments', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_A.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_A.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_A.set_xlim( x1, x2 )
ax_A.set_ylim( 175, 275 )

# Do not label x-axis ticks on Panel A
ax_A.xaxis.set_ticklabels([])

# Fill in Panel A with a filled countour of the RLTRN field
terrain_contour = ax_A.contourf( X, Y, zs, levels = np.arange( 0, t_max, t_int ), cmap = mpl.colormaps['copper'].reversed(),
                                 extend = 'max', alpha = alpha )

# Create callable Panel A objects, define label, and set plotting limits
ax_B = fig.add_subplot( axes['B'] )
ax_B.set_title( 'b: RLTRN Cross-Sections', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_B.set_xlabel( 'Zonal Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_B.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )
ax_B.set_xlim( x1, x2)
ax_B.set_ylim( -50, 1000 )
# ax_B.set_yticks( np.arange( 0, 1500, 250) )

#------------------------------------------------------------------------------
# Begin: Loop through each simulation
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Construct filename string (!!! REQUIRES Appropriate link to directory containing RLTRN CSVs !!!)
    filename = wdir + sim[i] + '_model_output_stats.csv'
    
    # Read in current simulation spreadsheet summary dataset
    ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km, zvort_surface, zvort_500m, 
      zvort_1km, zvort_3km, zvort_5km, uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
      mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir, shear1km, shear1km_dir, shear3km, shear3km_dir, 
      shear6km, shear6km_dir, srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, 
      cape3, lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km, zA500m, zA1km, zA3km, 
      zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth, cp_intensity, cp_area ) = read_stats_data( filename )
    
    
    # 1-sigma smoothing applied to points to create clean storm tracks
    uhx_smooth = gauss_smoother( uhx, 1.25 )
    uhy_smooth = gauss_smoother( uhy, 1.25 )

    # Plot the storm tracks as a solid line
    if( i == 1 ):
        storm_track = ax_A.plot( uhx_smooth[ 14:dis_end[i] ], uhy_smooth[ 14:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
    else:
        storm_track = ax_A.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], uhy_smooth[ sup_start[i]:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
        
    # Demarcate where storm begins and ends supercellular mode
    ax_A.scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    ax_A.scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )

    # Demarcate where storm begins and ends linear mode (if applicable)
    #------------------------------------------------------------------------------
    if(linear[i]> 0):
        ax_A.scatter( uhx_smooth[ sup_end[i]+1  ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
        ax_A.scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10  )
    #------------------------------------------------------------------------------

    # Annotate BSS Transitions and add colored halos around them
    if( i == 1 ):
        txt1 = ax_A.text( uhx_smooth[ sup_start[i] ]-35, uhy_smooth[ sup_start[i] ]-10, s = 'BSS0', fontsize = 15, fontweight = 'bold', color = 'k' )
        txt2 = ax_A.text( uhx_smooth[31], uhy_smooth[31]+12.5, s = 'BSS0', fontsize = 15, fontweight = 'bold' , color = 'k' )
        txt3 = txt2
    else:
        txt1 = ax_A.text( uhx_smooth[ sup_start[i] ]-30, uhy_smooth[ sup_start[i] ]-2.5, s = 'BSS0', fontsize = 15, fontweight = 'bold', color = 'k' )
        txt2 = ax_A.text( uhx_smooth[36]-15, uhy_smooth[36]-2, s = 'BSS1', fontsize = 15, fontweight = 'bold' , color = 'k' )
        txt3 = ax_A.text( uhx_smooth[ dis_end[i] ]-25, uhy_smooth[ dis_end[i] ]-5, s = 'BSS1->BSS2', fontsize = 15, fontweight = 'bold', color = 'k' )
    plt.setp( [txt1, txt2, txt3], path_effects = [ PathEffects.withStroke( linewidth = 3, foreground = color[i], alpha = 0.65 ) ] )  
    

    # Open respective netCDF with UH swath data w/ Xarray
    filename2 = wdir + 'BSS_RLTRN_Final_Output_Files/' + sim[i] + '_cm1out_000093.nc' 
    ds = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
    
    # Get the 2-5 km Updraft Helicity Swath (m^2/s^2)
    shs = ds.metpy.parse_cf( 'shs' ).isel( time = 0, yh = slice( 0, ny), xh = slice( 0, nx )  )

    # UH Swath Contours (uses same XY grid from terrain file [line 122] )
    cs1 = ax_A.contourf( X, Y, scipynd.gaussian_filter( shs[ :, : ], uh_sigma ), levels = np.arange( uh_bot, uh_top+uh_int, uh_int ), 
                          alpha = alpha, extend = 'max', cmap = mpl.cm.gist_gray )
    
    
    # Run terrain profile through guassian smoother    
    zs_alt = gauss_smoother( zs * 1000, 0.75 )
    
    # Plot terrain profiles underneath xy locations of storm tracks
    ax_B.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color ='k', linewidth = 1.5, linestyle = '-' )
    
    # Fill in the area under the terrain curve
    ax_B.fill_between( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color = color[i], alpha = 0.75, label = sim[i].upper() )
    
    # Demarcate start/end supercell mode on terrain curve
    ax_B.scatter( uhx_smooth[ sup_start[i] ], zs_alt[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
    ax_B.scatter( uhx_smooth[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

    # Demarcate start/end linear mode on terrain curve (if applicable)
    #------------------------------------------------------------------------------
    if( linear[i] > 0 ):
        ax_B.scatter( uhx_smooth[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        ax_B.scatter( uhx_smooth[ linear[i] ], zs_alt[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
    #------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# End: Loop through each simulation
#------------------------------------------------------------------------------

# Make labels for storm-mode and create legends
sup_lab = mlines.Line2D( [], [], color = 'k', marker = 'v', label = 'Supercell' )
lin_lab = mlines.Line2D( [], [], color = 'k', marker = 's', label = 'Linear' )
legend_1 = ax_A.legend( title = 'Simulations', loc = 'upper left', facecolor = 'lightgrey', fontsize = 14 )
legend_2 = ax_B.legend( title = 'Cross-Sections', loc = 'upper left', facecolor = 'lightgrey', fontsize = 14 )
ax_A.add_artist( legend_1 )
ax_A.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right', facecolor = 'lightgrey', fontsize = 14 )

# Add grids
ax_A.grid()
ax_B.grid()

# Include color bar legend (Terrain) and set label
cbar  = plt.colorbar( terrain_contour, location = 'right', shrink = 0.85,
                      ticks = np.arange( 0, t_max, t_int ), drawedges = True, ax = ax_A ) 
cbar.ax.set_title( 'Z (m)', fontsize = 16, fontweight = 'bold', pad = 20.0 )


# Include color bar legend (UH Swath ) and set label
cbar2  = plt.colorbar( cs1, location = 'right', shrink = 0.85, ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                       drawedges = True, ax = ax_B ) 
cbar2.ax.set_title( '\tUH (m$^{2}$ s$^{-2}$)', fontsize = 16, fontweight = 'bold', pad = 20.0 )


# Save figure
fig.savefig( fname = 'BSS_RLTRN_storm_track_w_swath_V2.jpeg', dpi = 300, bbox_inches = "tight" )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# End Main Script
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Confirm that script successfully ran
print( "\n{} successfully completed!".format( program ) )

# Report the time required to run the function
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )
