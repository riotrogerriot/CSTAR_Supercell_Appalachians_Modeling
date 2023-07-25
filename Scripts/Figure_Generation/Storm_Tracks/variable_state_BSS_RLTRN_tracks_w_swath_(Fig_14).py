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
import xarray as xr
from metpy.units import units
import matplotlib.patheffects as PathEffects
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
program = 'CM1_BSS_RLTRN_Analysis_Storm_Tracks_w_Swaths.py'

nx = 2400
ny = 2000
uh_bot = 500
uh_top = 4500
uh_int = 1000
uh_alpha = 0.35
uh_sigma = 5.0

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( '\nBegin Program: {}'.format ( program ) )

# Arrays to loop through simulations and assign color
sim = [ 'bss_nc_rltrn_all', 'bss_cs_rltrn_all' ]
color = [ 'tab:blue', 'tab:red' ]


# Arrays for demarcating storm mode and track
sup_start = [ 43, 39 ]
sup_end = [  69, 56 ]
dis_end = [ 80, 57 ]
linear = [ 0, 0 ]


# Adjust all timestamps by 2 hrs to acccount for delayed CI
for x in range (0, len( sup_start ) ):
    sup_start[x] = sup_start[x] - 25
    sup_end[x] = sup_end[x] - 25
    dis_end[x] = dis_end[x] - 25
    if( linear[x] > 0 ):
        linear[x] = linear[x] - 25
        
# Open terrain file for background data (!!! REQUIRES OUTPUT FILE !!!)
terr_file = '/Users/roger/Desktop/BSS_CS_RLTRN_ALL_plots/cm1out_000093.nc' 
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
                .BBB.
                .CDE.
                .FGH.
                """
                
fig, axes =  plt.subplot_mosaic(
                                mosaic = figure_mosaic,
                                figsize = ( 22, 22.5 ),
                                tight_layout = True
                               )

# Create callable axes objects
ax_A = fig.add_subplot( axes['A'] )
ax_A.set_title( 'a: Realistic Terrain: Variable-State (BSS) Simulations', fontsize = 20, fontweight = 'bold' )
ax_A.set_xlabel( '', fontsize = 18, fontweight = 'bold' )
ax_A.set_ylabel( 'Meridional Distance (km)', fontsize = 18, fontweight = 'bold' )

ax_A.set_xlim( 125, 350 )
ax_A.set_ylim( 185, 275 )

ax_A.xaxis.set_ticklabels([])

# Terrain Contours
terrain_contour = ax_A.contourf( 
                                X, Y, zs,
                                levels = np.arange( 0, 2250, 250), 
                                cmap = plt.cm.get_cmap( "copper" ).reversed(),
                                alpha = 0.35
                               )

ax_B = fig.add_subplot( axes['B'] )
ax_B.set_xlim( 125, 350 )
ax_B.set_ylim( -50, 1250 )
ax_B.set_yticks( np.arange( 0, 1500, 250) )

ax_B = fig.add_subplot( axes['B'] )
ax_B.set_title( 'b: Realistic Terrain Cross-Sections', fontsize = 20, fontweight = 'bold', loc = 'left' )
ax_B.set_xlabel( 'Zonal Distance (km)', fontsize = 18, fontweight = 'bold' )
ax_B.set_ylabel( 'Elevation (m)', fontsize = 18, fontweight = 'bold' )

# # Add white space for manually adding in Soundings
# ax_B = fig.add_subplot( axes['B'] )
# ax_B.axis('off')

# Begin: Loop through each simulation
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Construct filename string (!!! REQUIRES Directory containing BSS_RLTRN CSVs !!!)
    filename = '/Users/roger/Desktop/' + sim[i].upper() + '_plots/model_output_stats.csv'
    
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
    
        #########################################
        ### Denote RLTRN BSS times!!! ###
        #########################################
        
        txt1 = ax_A.text( uhx_smooth[ sup_start[i] ] -25, uhy_smooth[ sup_start[i] ] -1, s = 'BSS0', fontsize = 18, fontweight = 'bold', color = 'k' )
        if( sim[i] == 'bss_cs_rltrn_all' ):
            txt2 = ax_A.text( uhx_smooth[33]-15, uhy_smooth[33]+2.5, s = 'BSS0->BSS1', fontsize = 18, fontweight = 'bold' , color = 'k' )
        else:
            txt2 = ax_A.text( uhx_smooth[37]-20, uhy_smooth[37]+5, s = 'BSS0->BSS1', fontsize = 18, fontweight = 'bold' , color = 'k' )
        # txt3 = ax_A.text( uhx_smooth[61], uhy_smooth[61], s = 'BSS1->BSS2', fontsize = 18, fontweight = 'bold', color = 'k' )
        # txt4 = ax_A.text( uhx_smooth[73], uhy_smooth[73], s = 'BSS2', fontsize = 18, fontweight = 'bold', color = 'k' )

        # Open .NC with swath data (!!! REQUIRES OUTPUT FILE !!!)
        filename2 = '/Users/roger/Desktop/' + sim[i].upper() + '_plots/cm1out_000093.nc'
        
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
        ax_B.plot( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color ='k', linewidth = 1.5, linestyle = '-', label = sim[i].upper() + ' Elevation' )
        
        # Fill in the area under the terrain curve
        ax_B.fill_between( uhx_smooth[ sup_start[i]:dis_end[i] ], zs_alt[ sup_start[i]:dis_end[i] ], color = color[i], alpha = 0.75 )
        
        # Demarcate start/end supercell mode
        ax_B.scatter( uhx_smooth[ sup_start[i] ], zs_alt[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        ax_B.scatter( uhx_smooth[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

        # Demarcate start/end linear mode
        if( linear[i] > 0 ):
            ax_B.scatter( uhx_smooth[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
            ax_B.scatter( uhx_smooth[ linear[i] ], zs_alt[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

    
        plt.setp( 
                 [txt1, txt2], path_effects = [
                 PathEffects.withStroke( linewidth = 3, foreground = color[i], alpha = 0.75 ) ]
                )  
     

sup_lab = mlines.Line2D( [], [], color = 'k', marker = 'v', label = 'Supercell' )
lin_lab = mlines.Line2D( [], [], color = 'k', marker = 's', label = 'Linear' )

legend_1 = ax_A.legend( title = 'Simulations', loc = 'upper right' )
ax_A.add_artist( legend_1 )

ax_A.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, lin_lab], loc = 'lower right' )

ax_A.grid()
ax_B.grid()

# Include color bar legend
cbar  = plt.colorbar( 
                     terrain_contour,
                     location = 'right',
                     shrink = 0.8,
                     fraction = 0.075,
                     # pad = 0.05,
                     # anchor = (0.5, 1.0),
                     ticks = np.arange( 0, 2250, 250 ),
                     drawedges = False,
                     ax = ax_A
                    ) 

# Include color bar legend
cbar2  = plt.colorbar( 
                      cs1,
                      location = 'right',
                      shrink = 0.8,
                      fraction = 0.075,
                      # pad = 0.015,
                      # anchor = (0.5, 1.0),
                      ticks = np.arange( uh_bot, uh_top+uh_int, uh_int ),
                      drawedges = True,
                      ax = ax_B
                    ) 


# Set colorbar label 
cbar.ax.set_title( 'Z (m)', fontsize = 16, fontweight = 'bold', pad = 20.0 )
cbar2.ax.set_title( '\tUH (m$^{2}$s$^{-2}$)', fontsize = 16, fontweight = 'bold', pad = 20.0 )

#-----------------------------------------------------------------------------
# BSS Sounding Section
#-----------------------------------------------------------------------------

# (!!! Requires appropriate directory containing Initial Conditions Model Output Files !!!)
output_dir = '/scratch/rriggin/bss_cm1r20.3_fixed/IC_Output_Files/'
sim_type = [ 'Crosser', 'Non-Crosser' ]
sound_type = [ 'Upstream', 'Peak', 'Downstream', 'Upstream', 'Peak', 'Downstream' ]
bss_time = [ 't = BSS0', 't = BSS1', 't = BSS2', 't = BSS0', 't = BSS1', 't = BSS2' ]

# List for looping through each axis
current_ax = [ 'C', 'D', 'E', 'F', 'G', 'H' ]
current_ax_label = [ 'c', 'd', 'e', 'f', 'g', 'h' ]

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
            fname = 'BSS_RLTRN_storm_tracks_w_swath.jpeg',
            dpi = 300,
            bbox_inches = "tight"
           )
