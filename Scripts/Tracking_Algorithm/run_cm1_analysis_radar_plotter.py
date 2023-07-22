#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 23:04:23 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Function: cm1_radar_plotter
#
#   Update Records:
#       (3/21/22) Script Created
#
# Notes:
#   Need to add cross-section reference line once  plotting dimensions are finalized!!!
#       Use discrete color choice for this
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
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import metpy.plots as plots


#-----------------------------------------------------------------------------
# Begin: Define cm1_radar_plotter
#-----------------------------------------------------------------------------
def cm1_radar_plotter( 
                      ds, fig, ax, nx, ny,
                      x1, x2, y1, y2, z,
                      sx, sy, uh_x, uh_y, skip_val,
                      shortname, current_model_time,
                      terrain = False, peak_elv = 0.0, cint = 0.0,
                      id_terrain = False, peak_pos = 350.0, 
                      knots = False, meso_box = False, mx = 0, my = 0,
                      mdx = 0.0, mdy = 0.0
                     ):

    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the u-component of the wind (m/s)    
    u = ds.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
     # Get the v-component of the wind (m/s)   
    v = ds.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    # Get the radar reflectivity (dBZ)
    dbz = ds.metpy.parse_cf( 'dbz' ).isel( time = 0, zh = z, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get potential temperature (K)
    thpert = ds.metpy.parse_cf( 'thpert' ).isel( time = 0, zh = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the terrain field if simulation includes terrain (m)
    if( terrain == True ):
        zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the 2-5 km Updraft Helicity (m^2/s^2)
    uh = ds.metpy.parse_cf( 'uh' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
     # Convert winds to knots if desired
    #----------------------------------
    if( knots == True ):
        u = u.pint.to( 'knots' )
        v = v.pint.to( 'knots' )
        wind_units = 'kts'
    else:
        wind_units = 'ms$^{-1}$'
    #----------------------------------
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get the x-axis coordinates (km)
    xh = ds.coords[ 'xh' ].isel( xh = slice( 0, nx ) ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)      
    yh = ds.coords[ 'yh' ].isel( yh = slice( 0, ny) ).values * units( 'km' )
    
    
    # Logic to determine if x-axis should be peak-relative or not
    #------------------------------------------------------------
    if( id_terrain == True ):
        
        # Convert x-coordinates from grid-relative to peak-relative coordinates
        xh = xh - peak_pos * units( 'km' )
        
        # Change x-axis label to reflect coordinate shift
        ax.set_xlabel( "Distance from Peak (km)" , fontsize = 18, fontweight = 'bold' )
        
        sx = int( sx - peak_pos )
        
        uh_x = int( uh_x - peak_pos )
        
    # Label the axis as grid-relative    
    else:
        ax.set_xlabel( "Zonal Distance (km)" , fontsize = 18, fontweight = 'bold' )
    #------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # End: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Plotting Data
    #-----------------------------------------------------------------------------
    
    # Create a plotting grid based on the x,y coordinates
    X,Y = np.meshgrid( xh, yh )
    
    # Begin: Terrain Logic
    #---------------------
    
    # If simulation includes terrain
    if( terrain == True ):
        
        # Mask out lowest level of terrain (Needed for Realistic Terrain, Does not affect Idealized Terrain)
        trn_mask = np.ma.array( zs, mask = zs < 0.01 )
        
        # Create filled contours of the terrain field
        cf_trn = ax.contourf(
                             X,Y,
                             trn_mask,
                             levels = np.arange( 0, peak_elv + cint, cint ),
                             cmap = plt.cm.get_cmap( 'copper' ).reversed(),
                             alpha = 0.35
                            )
        
        # Begin: Logic for terrain contour lines
        #---------------------------------------
        if( id_terrain == False ):
            
            # Add contour lines for realistic terrain
            c_trn = ax.contour(
                                X,Y,
                                trn_mask,
                                levels = np.arange( 0, peak_elv + (cint*5), cint*5 ),
                                linewidths = 0.5,
                                alpha = 0.5
                               )           
        else:
            
            # Add contour lines for idealized terrain
            c_trn = ax.contour(
                               X,Y,
                               trn_mask,
                               levels = np.arange( 0, peak_elv + cint, cint ),
                               colors = 'saddlebrown',
                               linewidths = 1.0,
                               alpha = 1.0
                              )
    

        # End: Logic for terrain contour lines
        #---------------------------------------
    
    else:
        # Set background color to match 0m elevation
        ax.set_facecolor( '#fff0dd' )
            
    # End: Terrain Logic
    #---------------------
            
    
    # Mask out reflectivity values below 0 dBZ
    dbz_mask = np.ma.array( dbz , mask = dbz <= 10.0 )
    
    # Add reflectivity as a filled contour
    cf = ax.contourf(
                     X,Y,
                     dbz_mask, 
                     levels = np.arange( start = -10, stop = 80, step = 5.0 ),
                     cmap = plots.ctables.registry.get_colortable( "NWSReflectivity" ),
                     alpha = 0.60,
                     norm = Normalize(-10, 80)
                    )
    
    # Generate a reflectivity colorbar
    cbar = fig.colorbar ( cf, pad = 0.05 )
    
    # Modify the colorbar label
    cbar.set_label( label = 'Radar Reflectivity (dBZ)', fontsize = 18, fontweight = 'bold' )
    
    
    # Add UH Contours    
    cs1 = ax.contourf(
                      X,Y,
                      uh, 
                      levels = np.arange( 200, 6000, 500 ), 
                      colors = "black",
                      alpha = 0.40
                     )
    
    # Create legend label for UH contour
    cs1_lab = mpatches.Patch( color = 'k', alpha = 0.4, label = 'UH > 200 m$^2$s$^{-2}$' )


    # Slicing index used for clean wind vector plotting
    skip = ( slice( None, None, skip_val ), slice( None, None, skip_val ) )
    
    # Add wind vectors
    q = ax.quiver( 
                  X[skip],Y[skip],
                  u[skip], v[skip],
                  color = 'black', units = 'xy', angles = 'uv',
                  scale = 2.5, scale_units = 'xy', pivot = 'tail', alpha = 0.65
                 )
  
    # Add reference vector
    ax.quiverkey( 
                 Q = q,
                 X = 1.09, 
                 Y = 1.025, 
                 U = 10.0,
                 label = '10 '+ wind_units,
                 labelpos = 'N',
                 color = 'black',
                )
    
    
    # Add gust front contour as inferred by the surface -3 and -1K thpert contour
    cs2 = ax.contour(
                     X,Y,
                     thpert, levels = np.arange(-3, 0, 2), 
                     linestyles = "--" , linewidths = 2.0, 
                     colors = [ 'royalblue', 'blue' ], alpha = 0.75
                    )
    
    # Create legend label for gust front contour
    cs2a_lab = mlines.Line2D( 
                            [], [], color = 'royalblue',
                            linestyle = '--', linewidth = 2.0,
                            label = '-3 K \u03F4 Pert.'
                           )
    
    # Create legend label for gust front contour
    cs2b_lab = mlines.Line2D( 
                            [], [], color = 'blue',
                            linestyle = '--', linewidth = 2.0,
                            label = '-1 K \u03F4 Pert.'
                           )
    
    # Label inflow sounding point
    snd_xy = ax.plot(
                     sx, sy,
                     marker = '*', markersize = 20,
                     color = 'yellow', markeredgecolor = 'k',
                     label = 'Sounding Loc.'
                    )
    
    # Label UH maxima point
    uh_xy = ax.plot(
                    uh_x, uh_y,
                    marker = 'X', markersize = 15, color = 'aqua', 
                    markeredgecolor = 'k', label = 'UH Maxima'
                   )
    
    # Cross-Section reference line
    cross_x = ax.axhline( 
                         y = uh_y, color = 'k',
                         linestyle = '--', linewidth = 2.5,
                         label = 'Cross-Section'
                        )
    
    # Set y-axis label
    ax.set_ylabel( "Meridional Distance (km)" , fontsize = 18, fontweight = 'bold' )
    
    # Right Hand Plot Title (Model Desc)
    plt.title( 
              "\nCM1-{}: t = {} min. \na: Near-Surface Analysis".format( shortname, current_model_time),
              loc = "left", fontsize = 14, fontweight = 'bold' 
             )       
    
    # Right Hand Plot Title (Model Desc)
    plt.title( 
              "1 km Radar Reflectivity > 10 dBZ \n0.05 km Wind Vectors ({})".format( wind_units ),
              loc = "right", fontsize = 12, fontweight = 'bold' 
             )  
    
    # Add a grid to the axis
    # ax.grid()
    
    if( meso_box == True ):
        p = mpatches.Rectangle( ( xh[mx].m, yh[my].m ), mdx, mdy, fill = False, label = 'Area Window' )
        ax.add_patch( p )
    
    # Create Legend
    ax.legend( 
              handles = [ cs1_lab, cs2a_lab, cs2b_lab, snd_xy[0], uh_xy[0], cross_x, p ],
              loc = 1, facecolor = 'white', 
              framealpha = 0.75, fontsize = 14
             )
    
    # if( x2 >=  nx ):
    #     x2 = nx -1
    # if( y2 >= ny ):
    #     y2 = ny - 1
        
    # Set plot windows
    ax.set_xlim( xh[x1], xh[x2] )
    ax.set_ylim( yh[y1], yh[y2] )
    
    
    
    #-----------------------------------------------------------------------------
    # End: Plotting Data
    #-----------------------------------------------------------------------------
    
    # Return the plot
    return( ax )

#-----------------------------------------------------------------------------
# End: Define cm1_radar_plotter
#-----------------------------------------------------------------------------

