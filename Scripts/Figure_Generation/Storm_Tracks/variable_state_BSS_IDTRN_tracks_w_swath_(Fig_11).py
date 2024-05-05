#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### NEED BSS ANNOTATIONS


# Import external modules
#------------------------------------------------------------------------------
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
def gauss_smoother( arr, smooth = 6 ) :
    ret_arr = scipynd.gaussian_filter1d( arr, smooth )
    return( ret_arr )


# Begin Domain Height Array Function (Default values reflect CM1 model domain for this study)
#-------------------------------------------------------------------
def witch_of_agnesi( nx = 2401, ny = 1601, dx = 250.0, h0 = 750.0, x0 = 350000.0, a = 50000.0 ):

    x = np.arange( 0, nx )                        # gridpoint array units: 0
    y = np.arange( 0, ny )                        # gridpoint array units: 0
    hx = np.zeros( shape = (nx) )                 # Elevation array units: m 
    
    # Compute x & y over the model domain
    for i in range( 0, nx ):
        x[i] = x[i] * dx
    for i in range( 0, ny ):
        y[i] = y[i] * dx
        
    # Compute height over the model domain (Witch of Agnesi)
    for i in range( 0, len(hx) ):
        hx[i] = (h0) / ( 1 + ( (x[i]-x0) / a ) ** 2 )
            
    HX,HY = np.meshgrid( hx, y )
    
    # Return the height array back to user
    return( x, y, hx, HX, HY )


#------------------------------------------------------------------------------
# End: User-defined Functions
#------------------------------------------------------------------------------


# Program name, main directory, and sub-directories required to run this script
program = 'CM1_BSS_Analysis_Storm_Tracks_w_Swath.py'

# OS boolean (0 = MacOS, 1 = Windows)
os_type = 0
if( os_type == 0 ):
    wdir = '/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR_Modeling_Project/Simulations/Variable_State_Sims/'
else:
    wdir = 'C:/Users/rriggin/OneDrive - University of North Carolina at Charlotte/CSTAR_Modeling_Project/Simulations/Variable_State_Sims/'


# User-defined parameters
#------------------------------------------------------------------------------
nx = 2400             # Num of x-gridpoints
ny = 1600             # Num of y-gridpoints
uh_bot = 500          # Min UH plot value (m2/s2)
uh_top = 4500         # Max UH plot value (m2/s2)
uh_int = 1000         # UH plotting interval (m2/s2)
uh_sigma = 5          # UH smoothing parameter
x1 = -250             # X-axis min (km)
x2 = 50               # X-axis max (km)
alpha = 0.45          # Contourf transparency parameter
ax = ['A', 'B', 'C' ] # Array for varrying axis plotting within simulation loop
ax_i = 0              # Index for calling appropriate axis
#------------------------------------------------------------------------------


# Simulation-Specific Parameters
#------------------------------------------------------------------------------
sim = [ 'bss_cs_ctl' , 'bss_cs_trn', 'bss_nc_ctl', 'bss_nc_trn' ] # Sim names
color = [ 'darkorange', 'tab:red', 'lightskyblue', 'tab:blue' ]   # Sim colors
sup_start = [ 23, 28, 21, 23 ]                                    # Supercell Start Index
sup_end = [ 45, 45, 60, 46 ]                                      # Supercell End Index
dis_end = [ 53, 53, 72, 61 ]                                      # Dissipation Time Index
linear = [ 0, 0, 0, 0 ]                                           # Upscale Growth Index
#------------------------------------------------------------------------------


#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Begin Main Script
#-------------------------------------------------------------------
#-------------------------------------------------------------------

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( '\nBegin Program: {}'.format ( program ) )

x, y, hx, HX, HY = witch_of_agnesi()

# Convert to peak-relative
x = x - 350000.0

# Set tick size (Must be before calling plot object)
plt.rcdefaults()
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

# Plot layout layout
figure_mosaic = """
                .AAA.
                .AAA.
                .BBB.
                .BBB.
                .CCC.
                """
                
# Create figure with axes based on the layout
fig, axes =  plt.subplot_mosaic( mosaic = figure_mosaic, figsize = ( 16, 12.5 ), tight_layout = True )

# Create callable Panel A objects, define label, and set plotting limits
ax_A = fig.add_subplot( axes['A'] )
ax_A.set_title( 'a: IDTRN Variable-State (BSS) Crossing Environment', fontsize = 16, fontweight = 'bold', loc = 'left' )
ax_A.set_xlabel( '', fontsize = 18, fontweight = 'bold', loc = 'left' )
ax_A.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_A.set_xlim( x1, x2 )
ax_A.set_ylim( 190, 290 )

# Fill in Panel A with a filled countour of the IDTRN field
terrain_fill = ax_A.contourf( x/1000.0, y/1000.0, HX, levels = np.arange(0, 750.0 + 75.0, 75.0), 
                              cmap = mpl.colormaps['copper'].reversed(), alpha = alpha )

# Create callable Panel B objects, define label, and set plotting limits
ax_B = fig.add_subplot( axes['B'] )
ax_B.set_title( 'b: IDTRN Variable-State (BSS) Non-Crossing Environment', fontsize = 16, fontweight = 'bold', loc = 'left' )
ax_B.set_xlabel( '', fontsize = 18, fontweight = 'bold',  loc = 'left' )
ax_B.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_B.set_xlim( x1, x2 )
ax_B.set_ylim( 110, 210 )

# Do not label x-axis ticks on Panel A or B
ax_A.xaxis.set_ticklabels([])
ax_B.xaxis.set_ticklabels([])

# Fill in Panel A with a filled countour of the IDTRN field
terrain_fill = ax_B.contourf( x/1000.0, y/1000.0, HX, levels = np.arange(0, 750.0 + 75.0, 75.0), 
                              cmap = mpl.colormaps['copper'].reversed(), alpha = alpha )


# Create callable Panel C objects, define label, and set plotting limits
ax_C = fig.add_subplot( axes['C'], aspect = 'auto' )
ax_C.set_title( 'c: Idealized Terrain Cross-Section', fontsize = 16, fontweight = 'bold', loc = 'left' )
ax_C.set_xlabel( 'Distance from Peak (km)', fontsize = 18, fontweight = 'bold' )
ax_C.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )
ax_C.set_xlim( x1, x2 )
ax_C.set_ylim( 0, 1000 )

# Make background color sky blue
ax_C.set_facecolor( 'skyblue' )

# Plot the terrain profile and fill in the area under the terrain curve
terrain_profile = ax_C.plot( x/1000.0, hx, linewidth = 2.5, fillstyle = 'full', color = 'k' )
ax_C.fill_between( x/1000.0, hx, color = 'saddlebrown' )

# Add grids
ax_A.grid()
ax_B.grid()
ax_C.grid()

#------------------------------------------------------------------------------
# Begin: Loop through each simulation ### REMOVE UH Swath logic once BSS_NC_TRN re-run is complete
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Variable axis logic( i < 1 = CS in Panel A, i > 1 = NC in Panel B )
    if( i > 1 ):
        ax_i = 1
    
    # Construct filename string (!!!REQUIRES cleaned model output CSVs!!!)
    filename = wdir + 'Stats_Spreadsheets/' + sim[i] + '_model_output_stats.csv'
    
    # Read in current simulation spreadsheet summary dataset
    ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km, zvort_surface, zvort_500m, 
     zvort_1km, zvort_3km, zvort_5km, uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
     mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir, shear1km, shear1km_dir, shear3km, shear3km_dir, 
     shear6km, shear6km_dir, srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, 
     cape3, lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km, zA500m, zA1km, zA3km, 
     zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth, cp_intensity, cp_area ) = read_stats_data( filename )
   
    # 1-sigma smoothing applied to points to create clean storm tracks
    uhx_smooth = gauss_smoother( uhx - 350.0, 2. )
    uhy_smooth = gauss_smoother( uhy, 2. )
    
    # Plot the storm tracks as a solid line
    storm_track = axes[ ax[ax_i] ].plot( uhx_smooth[ sup_start[i]:dis_end[i] ], uhy_smooth[ sup_start[i]:dis_end[i] ], color = color[i], linewidth = 2.5, linestyle = '-', label = sim[i].upper() )
    
    # Demarcate where storm begins and ends supercellular mode
    axes[ ax[ax_i] ].scatter( uhx_smooth[ sup_start[i] ], uhy_smooth[ sup_start[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )
    axes[ ax[ax_i] ].scatter( uhx_smooth[ sup_end[i] ], uhy_smooth[ sup_end[i] ], color = color[i], marker ='v', s = 100, edgecolor = 'k', zorder = 10 )

    # Demarcate where storm begins and ends linear mode (if applicable)
    #------------------------------------------------------------------------------
    if(linear[i]> 0):
        axes[ ax[ax_i] ].scatter( uhx_smooth[ sup_end[i]+1  ], uhy_smooth[ sup_end[i]+1 ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10 )
        axes[ ax[ax_i] ].scatter( uhx_smooth[ linear[i] ], uhy_smooth[ linear[i] ], color = color[i], marker = 's', s= 80, edgecolor = 'k', zorder = 10 )
    #------------------------------------------------------------------------------

    # Annotate BSS Transitions and add colored halos around them
    if( i == 3 ):
        txt1 = axes[ ax[ax_i] ].text( uhx_smooth[ sup_start[i] ]-20, uhy_smooth[ sup_start[i] ]+6, s = 'BSS0', fontsize = 15, fontweight = 'bold', color = 'k' )
        txt2 = axes[ ax[ax_i] ].text( uhx_smooth[36]-30, uhy_smooth[36]+6, s = 'BSS1', fontsize = 15, fontweight = 'bold' , color = 'k' )
        txt3 = axes[ ax[ax_i] ].text( uhx_smooth[ dis_end[i] ]-30, uhy_smooth[ dis_end[i] ]+2.5, s = 'BSS1->BSS2', fontsize = 15, fontweight = 'bold', color = 'k' )
    else:
        txt1 = axes[ ax[ax_i] ].text( uhx_smooth[ sup_start[i] ]-32.5, uhy_smooth[ sup_start[i] ]-5, s = 'BSS0', fontsize = 15, fontweight = 'bold', color = 'k' )
        txt2 = axes[ ax[ax_i] ].text( uhx_smooth[36]-35, uhy_smooth[36]-10, s = 'BSS1', fontsize = 15, fontweight = 'bold' , color = 'k' )
        txt3 = axes[ ax[ax_i] ].text( uhx_smooth[ dis_end[i] ]-30, uhy_smooth[ dis_end[i] ]-5, s = 'BSS1->BSS2', fontsize = 15, fontweight = 'bold', color = 'k' )
    plt.setp( [txt1, txt2, txt3], path_effects = [ PathEffects.withStroke( linewidth = 3, foreground = color[i], alpha = 0.65 ) ] )  
    
    # Open .NC with swath data (!!!REQUIRES Final output file for UH Swath!!!)
    filename2 = wdir + '/Stats_Spreadsheets/BSS_IDTRN_Final_Output_Files/' + sim[i] + '_cm1out_000073.nc'
    
    # if( i < len(sim)-1 ):
    ds = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
      
    # Get the x and y-axis coordinates (km)
    xh = ds.coords[ 'xh' ].values * units( 'km' )
    xh = xh - 350.0 * units.km     
    yh = ds.coords[ 'yh' ].values * units( 'km' )
     
    # Create a plotting grid based on the x,y coordinates
    X,Y = np.meshgrid( xh, yh )
    
    # Get the 2-5 km Updraft Helicity Swath (m^2/s^2)
    shs = ds.metpy.parse_cf( 'shs' ).isel( time = 0, yh = slice( 0, ny), xh = slice( 0, nx )  )
    
    # UH Swath Contours (uses same XY grid from terrain file [line 122] )
    cs1 = axes[ ax[ax_i] ].contourf( X, Y, scipynd.gaussian_filter( shs[ :, : ], uh_sigma ), levels = np.arange( uh_bot, uh_top+uh_int, uh_int ), 
                      alpha = alpha, extend = 'max', cmap = mpl.cm.gist_gray )
    
#------------------------------------------------------------------------------
# End: Loop through each simulation
#------------------------------------------------------------------------------

# Make labels for storm-mode and create legends
sup_lab = mlines.Line2D( [], [], color = 'k', marker = 'v', label = 'Supercell' )
lin_lab = mlines.Line2D( [], [], color = 'k', marker = 's', label = 'Linear' )
legend_1 = ax_A.legend( title = 'Simulations', loc = 'upper left', facecolor = 'lightgrey', fontsize = 14 )
legend_2 = ax_B.legend( title = 'Simulations', loc = 'lower left', fontsize = 14, facecolor = 'lightgrey' )
ax_A.add_artist( legend_1 )
ax_B.add_artist( legend_2 )
ax_A.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right', facecolor = 'lightgrey', fontsize = 14 )
ax_B.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right', facecolor = 'lightgrey', fontsize = 14 )

# Include color bars (1: Terrain 2: UH 3: Placeholder to be removed via Photoshop) and set label
cbar  = plt.colorbar( terrain_fill, location = 'right', shrink = 0.85, ticks = np.arange( 0, 825.0, 75.0 ), drawedges = True, ax = ax_A ) 
cbar.ax.set_title( 'Z (m)', fontsize = 16, fontweight = 'bold', pad = 20.0 )
cbar2  = plt.colorbar( cs1, location = 'right', shrink = 0.85, ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ), drawedges = True, ax = ax_B ) 
cbar2.ax.set_title( '\tUH (m$^{2}$ s$^{-2}$)', fontsize = 16, fontweight = 'bold', pad = 20.0 )
cbar3  = plt.colorbar( cs1, location = 'right', shrink = 0.85, ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ), drawedges = True, ax = ax_C ) 

# Save figure
fig.savefig( fname = 'BSS_IDTRN_storm_tracks_w_swath.jpeg', dpi = 300, bbox_inches = "tight" )

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# End Main Script
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Confirm that script successfully ran
print( "\n{} successfully completed!".format( program ) )

# Report the time required to run the function
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )