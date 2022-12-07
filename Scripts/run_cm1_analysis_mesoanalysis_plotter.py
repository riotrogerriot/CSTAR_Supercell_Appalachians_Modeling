#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 21:38:01 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Function: cm1_mesoanalysis_plotter
#
#   Update Records:
#       (3/29/22) Script Created
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

# Notes about srh_cy
# This function utilizes Cython to compute SRH over the model domain for a given
# depth. Cython uses a combo of python and C to speed up computations. This function
# requires the C code to be compiled using a terminal command: python setup_srh.py build_ext --inplace
# before the code can be ran properly!
from srh_cy import srh




import pint
import pint_xarray
import cf_xarray

#-----------------------------------------------------------------------------
# Begin: Define cm1_cross_section_plotter
#-----------------------------------------------------------------------------
def cm1_mesoanalysis_plotter( ds, ax, x1, x2, y1, y2, z1, z2, shortname, current_model_time,
                              sx, sy, storm_motion, uh_x, uh_y, bb_box_dist, dx, dy, 
                              terrain = False, id_terrain = True, knots = True,
                              km1 = 6, km3 = 13, km5 = 18, skip_val = 100,
                              nx = 2400, ny = 1600, nz = 48, smooth = 8,
                              peak_pos = 350.0, peak_elv = 750.0, cint = 75.0 
                             ):

    #-----------------------------------------------------------------------------
    # Begin: Scale bb_box for plot
    #-----------------------------------------------------------------------------
    
    bb_box_dist = 2 * bb_box_dist
    
    x1 = int( uh_x - ( ( 0.50 * bb_box_dist * 1000.0 ) / dx ) )
    x2 = int( uh_x + ( ( 1.50 * bb_box_dist * 1000.0 ) / dx ) )
    y1 = int( uh_y - ( ( 0.75 * bb_box_dist * 1000.0 ) / dy ) )
    y2 = int( uh_y + ( ( 1.25 * bb_box_dist * 1000.0 ) / dy ) )
    
    if( x2 < nx ):
        x2 = x2
    else:
        x2 = nx-1
        
    if( y2 < ny ):
        y2 = y2
    else:
        y2 = ny-1
    
    #-----------------------------------------------------------------------------
    # End: Scale bb_box for plot
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the u-component of the wind (m/s)    
    u = ds.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
     # Get the v-component of the wind (m/s)   
    v = ds.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the radar reflectivity (dBZ)
    dbz = ds.metpy.parse_cf( 'dbz' ).isel( time = 0, zh = km1, yh = slice( 0, ny ), xh = slice( 0, nx ) )

    
    # Get the 2-5 km Updraft Helicity Swath (m^2/s^2)
    shs = ds.metpy.parse_cf( 'shs' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get CAPE (J/kg)
    cape = ds.metpy.parse_cf( 'cape' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get CIN (J/kg)
    cin = ds.metpy.parse_cf( 'cin' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # # Get pressure (Pa)
    # pressure = ds.metpy.parse_cf( 'prs' ).isel( time = 0,  zh = km1, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the terrain field if simulation includes terrain (m)
    if( terrain == True ):
        zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
        
        zs = zs.values
        
    else:
        zs = np.zeros( shape = ( ny, nx ), dtype = 'float32' )
    
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get the x-axis coordinates (km)
    xh = ds.coords[ 'xh' ].isel( xh = slice( 0, nx ) ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)      
    yh = ds.coords[ 'yh' ].isel( yh = slice( 0, ny ) ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)      
    zh = ds.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' )
    
    
    # Logic to determine if x-axis should be peak-relative or not
    #------------------------------------------------------------
    if( id_terrain == True ):
        
        # Convert x-coordinates from grid-relative to peak-relative coordinates
        xh = xh - peak_pos * units( 'km' )
        
        # Change x-axis label to reflect coordinate shift
        ax.set_xlabel( "Distance from Peak (km)" , fontsize = 18, fontweight = 'bold' )
        
        # sx = int( sx - peak_pos )
        
        # uh_x = int( uh_x - peak_pos )
        
    # Label the axis as grid-relative    
    else:
        ax.set_xlabel( "Zonal Distance (km)" , fontsize = 18, fontweight = 'bold' )
    #------------------------------------------------------------
    
    
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
    
    # pressure = pressure.metpy.quantify()
    # pressure = pressure.pint.to( 'hPa' )  
    
    #-----------------------------------------------------------------------------
    # End: Unit conversions required for MetPy functions
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Calculations
    #-----------------------------------------------------------------------------
    
    # # Compute Bunkers Storm Motions and Mean Wind for given inflow sounding 
    # rm, lm, mean_uv = mpcalc.bunkers_storm_motion( 
    #                                               pressure = pressure[:, sy, sx],
    #                                               u = u[:, sy ,sx],
    #                                               v = v[:, sy, sx],
    #                                               height = zh
    #                                              )
    
    # Convert zh into 3D array required for SRH function
    zh_grid = np.zeros( shape = (nz, ny, nx) )
    
    # Interpolate heights to a 3D domain required for SRH function
    #-------------------------------------------------------------
    for k in range( 0, nz ):
        zh_grid[k,:,:] = zh[k].m
        
        # Convert to meters
        zh_grid[k,:,:] = zh_grid[k,:,:] * 1000.0
    #-------------------------------------------------------------
    
    # Compute 0-1 km SRH using Cython function (Must convert from pint to np.array using .values method)
    srh1km = srh(
                 u = u.values,                        # units: m/s
                 v = v.values,                        # units: m/s
                 zh = zh_grid,                        # units: m
                 zs = zs,                             # units: m
                 u_storm = storm_motion[0].m,         # units: m/s
                 v_storm = storm_motion[1].m,         # units: m/s
                 layer = 1000.0                       # units: m
                )
    
    #-----------------------------------------------------------------------------
    # End: Calculations
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Plotting Data
    #-----------------------------------------------------------------------------
    
    # Convert winds to knots if desired (Must be after SRH call for proper calculation!)
    #----------------------------------
    if( knots == True ):
        u = u.pint.to( 'knots' )
        v = v.pint.to( 'knots' )
        storm_motion = storm_motion.to( 'knots' )
        wind_units = 'kts'
    else:
        wind_units = 'ms$^{-1}$'
    #----------------------------------
    
    # Create a plotting grid based on the x,y coordinates
    X,Y = np.meshgrid( xh, yh )
    
    # Add a grid to the axis
    ax.grid()
    
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
        
        # # Generate a reflectivity colorbar
        # cbar = fig.colorbar ( cf_trn, pad = 0.05 )
        
        # # Modify the colorbar label
        # cbar.set_label( label = 'Elevation (m)', fontsize = 18, fontweight = 'bold' )
        
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
            
    
    #-----------------------------------------------------------------------------
    # Begin: Plot data
    #-----------------------------------------------------------------------------
   
    # Add UH swath contours    
    cs4 = ax.contourf(
                      X,Y,
                      scipy.ndimage.gaussian_filter( input = shs, sigma = smooth ), 
                      levels = np.arange( 200, 8200, 4000 ), 
                      colors = "black",
                      alpha = 0.35
                     )
    
    
    # Slicing index used for clean wind vector plotting
    skip = ( slice( None, None, skip_val ), slice( None, None, skip_val ) )
    
    # Add SR wind vectors at surface
    q1 = ax.barbs( 
                 X[skip], Y[skip],
                 u[0,:,:][skip], v[0,:,:][skip],
                 length = 6.0, color = 'black', alpha = 0.75, label = 'Surface Wind {}'.format( wind_units )
                )
    
    # Add SR wind vectors at 3km
    q2 = ax.barbs( 
                 X[skip], Y[skip],
                 u[km3,:,:][skip], v[km3,:,:][skip],
                 length = 6.0, color = 'blue', alpha = 0.75, label = '3 km Wind {}'.format( wind_units )
                )
    
    # Add SR wind vectors at 5km
    q3 = ax.barbs( 
                 X[skip], Y[skip],
                 u[km5,:,:][skip], v[km5,:,:][skip],
                 length = 6.0, color = 'tab:red', alpha = 0.75, label = '5 km Wind {}'.format( wind_units )
                )
    
    # Add 0-1 km SRH contours
    cs1 = ax.contour(
                      X,Y,
                      scipy.ndimage.gaussian_filter( input = srh1km, sigma = smooth ),
                      levels = np.arange( 100, 600, 100 ), 
                      linestyles = "-" , linewidths = 1.5, 
                      colors = [ 'indigo' ], alpha = 0.75
                    )
    
    # Add CIN as hatched filled contour
    cf3 = ax.contourf(
                      X,Y, 
                      scipy.ndimage.gaussian_filter( input = cin, sigma = smooth ),
                      levels = np.arange( 50, 550, 250 ), 
                      hatches = [ '/' ], colors = 'tab:cyan',
                      alpha = 0.40
                     )
    
    # Add CIN contour line                
    cs3 = ax.contour(
                      X,Y, 
                      scipy.ndimage.gaussian_filter( input = cin, sigma = smooth ),
                      levels = np.arange( 50, 250, 200 ), 
                      linestyles = "--" , linewidths = 2.0, 
                      colors = [ 'b' ], alpha = 0.75
                    )
    
    # Add CAPE contours
    cs2 = ax.contour(
                      X,Y, 
                      scipy.ndimage.gaussian_filter( input = cape, sigma = smooth ),
                      levels = np.arange( 500, 3000, 500 ), 
                      linestyles = "--" , linewidths = 2.0, 
                      colors = [ 'crimson' ], alpha = 0.75
                    )
    
    # Mask out reflectivity values below 10 dBZ
    dbz_mask = np.ma.array( dbz , mask = dbz <= 10.0 )
    
    # Add reflectivity as a filled contour
    cf = ax.contourf(
                      X,Y,
                      dbz_mask, 
                      levels = np.arange( start = -10, stop = 80, step = 5.0 ),
                      cmap = plots.ctables.registry.get_colortable( "NWSReflectivity" ),
                      alpha = 0.45,
                      norm = Normalize(-10, 80)
                    )
    
    # # Label inflow sounding point
    # snd_xy = ax.plot(
    #                   sx, sy,
    #                   marker = '*', markersize = 25,
    #                   color = 'yellow', markeredgecolor = 'k',
    #                   label = 'Sounding Loc.'
    #                 )
    
    # # Label UH maxima point
    # uh_xy = ax.plot(
    #                 uh_x, uh_y,
    #                 marker = 'X', markersize = 10, color = 'aqua', 
    #                 markeredgecolor = 'k', label = 'UH Maxima'
    #                 )
    
    #-----------------------------------------------------------------------------
    # End: Plot Data
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Label contoured data
    #-----------------------------------------------------------------------------
    
    label_flag = True
    
    # Label the 0-1 km SRH contours
    try:
        cs1_lab = ax.clabel(
                            cs1,
                            levels = np.arange( 100, 600, 100 ),
                            fmt = '%3.0f',
                            fontsize = 10,
                            inline = True,
                            inline_spacing = 0.25,
                            colors = [ 'indigo' ],
                            use_clabeltext = True
                           )
    except ValueError:
        print( 'SRH too low to contour' )
        label_flag = False
        
    # Label CAPE contours
    try:
        cs2_lab = ax.clabel(
                            cs2,
                            levels = np.arange( 500, 3000, 500 ),
                            fmt = '%3.0f',
                            fontsize = 10,
                            inline = True,
                            inline_spacing = 0.25,
                            colors = [ 'crimson' ],
                            use_clabeltext = True
                           )
    except ValueError:
        print( 'CAPE too low to contour' )
        label_flag = False
        
    #-----------------------------------------------------------------------------
    # End: Label contoured data
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Create legend labels
    #-----------------------------------------------------------------------------
    
    if( label_flag == True ):
        # Create an SRH legend label
        cs1_leg = mlines.Line2D( 
                                [], [], color = 'indigo',
                                linestyle = '-', linewidth = 2.0,
                                label = '0-1 km SRH > 100 m$^2$s$^{-2}$'
                               )
        # Create an CAPE legend label
        cs2_leg = mlines.Line2D( 
                                [], [], color = 'crimson',
                                linestyle = '--', linewidth = 2.0,
                                label = 'CAPE > 500 J$kg^{-1}$'
                               )
        
    # Create an CIN legend label
    cs3_leg = mpatches.Patch( facecolor = 'tab:cyan', edgecolor = 'k', hatch = '/', 
                              alpha = 0.75, label = 'CIN > 50 J$kg^{-1}$' )
    
    # Create legend label for UH contour
    cs4_leg = mpatches.Patch( color = 'black', alpha = 0.35, 
                              label = 'UH Swath > 200 m$^2$s$^{-2}$' )
    
    #-----------------------------------------------------------------------------
    # End: Create legend labels
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Optimize plot
    #-----------------------------------------------------------------------------
    
    if( label_flag == True ):
        
        # Create a halo around contour labels
        plt.setp(
                 [cs1_lab, cs2_lab], path_effects =[ 
                 PathEffects.withStroke( linewidth = 3, foreground = 'w', alpha = 0.75 ) ] 
                )
    
    # Set y-axis label
    ax.set_ylabel( "Meridional Distance (km)" , fontsize = 18, fontweight = 'bold' )
    
    # Right Hand Plot Title (Model Desc)
    plt.title( 
              "CM1-{}: t = {} min.\nd: Environment Mesoanalysis".format(shortname, current_model_time),
              loc = "left", fontsize = 14, fontweight = 'bold' 
              )       
    
    # Right Hand Plot Title (Model Desc)
    plt.title( 
              "1.0 km Radar Reflectivity > 10 dBZ \n {}-Sigma Gaussian Smoothed Parameters".format( smooth ),
              loc = "right", fontsize = 12, fontweight = 'bold' 
             )  
    
    # Create Legend
    if( label_flag == True ):
        ax.legend( 
                  handles = [ cs1_leg, cs2_leg, cs3_leg, cs4_leg,  q1, q2, q3 ], #, uh_xy[0] ],
                  loc = 1, facecolor = 'white', 
                  framealpha = 0.75, fontsize = 14
                 )
    else:
        ax.legend( 
                  handles = [ cs3_leg, cs4_leg,  q1, q2, q3 ], #, uh_xy[0] ],
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
    # End: Optimize plot
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # End: Plotting Data
    #-----------------------------------------------------------------------------
    
    # Return the plot
    return( ax )

#-----------------------------------------------------------------------------
# End: Define cm1_radar_plotter
#-----------------------------------------------------------------------------





# nx = 2400
# ny = 1200
# nz = 48

# x1 = 750 - 150
# x2 = 1150 + 450
# y1 = 800 - 175
# y2 = 1100 + 225

# z1 = 0
# z2 = 48
# sx = -160
# sy = 195

# km1 = 6
# km3 = 13
# km5 = 18

# peak_pos = 350.0
# peak_elv = 750.0
# cint = 75.0
# id_terrain = True
# knots = True
# terrain = True
# skip_val = 100

# smooth = 8 

# filename = '/scratch/rriggin/cm1r20.3/original_sims_corrected_snd/crosser/cs_trn/cm1out_000037.nc'

# # Open the current terrain-interpolated netCDF file with xarray 
# ds = xr.open_dataset( filename, engine = "netcdf4", decode_cf = True )

# # Representative of the individual plots that compose the figure (Each letter is a separate subplot and . skips that position)
# figure_mosaic = """
#                   AA
#                   AA
#                 """

# # Create figure and axes objects
# fig, axes = plt.subplot_mosaic( 
#                                 mosaic = figure_mosaic, 
#                                 figsize = ( 11, 8.5 ), 
#                                 constrained_layout = True
#                               )

#   # Construct axis objects for the radar plot
# ax= fig.add_subplot( axes['A'] )


# #-----------------------------------------------------------------------------
# # Begin: Get variables from xarray dataset
# #-----------------------------------------------------------------------------

# # Data Variables
# #-----------------------------------------------------------------------------

# # Get the u-component of the wind (m/s)    
# u = ds.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh =  slice( z1, z2 ), 
#                                           yh = slice( y1, y2 ), xh = slice( x1, x2 ) )

#   # Get the v-component of the wind (m/s)   
# v = ds.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh =  slice( z1, z2 ), 
#                                           yh = slice( y1, y2 ), xh = slice( x1, x2 ) )
# # Get the radar reflectivity (dBZ)
# dbz = ds.metpy.parse_cf( 'dbz' ).isel( time = 0, zh = km1, 
#                                         yh = slice( y1, y2), xh = slice( x1, x2) )

# # Get the 2-5 km Updraft Helicity Swath (m^2/s^2)
# shs = ds.metpy.parse_cf( 'shs' ).isel( time = 0, yh = slice( y1, y2 ), xh = slice( x1, x2 ) )

# # Get CAPE (J/kg)
# cape = ds.metpy.parse_cf( 'cape' ).isel( time = 0, yh = slice( y1, y2 ), xh = slice( x1, x2 ) )

# # Get CIN (J/kg)
# cin = ds.metpy.parse_cf( 'cin' ).isel( time = 0, yh = slice( y1, y2 ), xh = slice( x1, x2 ) )

# # Get pressure (Pa)
# pressure = ds.metpy.parse_cf( 'prs' ).isel( time = 0, zh =  slice( z1, z2 ), 
#                                           yh = slice( y1, y2 ), xh = slice( x1, x2 ) )

# # Get the terrain field if simulation includes terrain (m)
# if( terrain == True ):
#     zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0,
#                                           yh = slice( y1, y2 ), xh = slice( x1, x2 ) )


# # Coordinate Variables
# #-----------------------------------------------------------------------------

# # Get the x-axis coordinates (km)
# xh = ds.coords[ 'xh' ].isel( xh = slice( x1, x2) ).values * units( 'km' )

# # Get the y-axis coordinates (km)      
# yh = ds.coords[ 'yh' ].isel( yh = slice( y1, y2) ).values * units( 'km' )

# # Get the y-axis coordinates (km)      
# zh = ds.coords[ 'zh' ].isel( zh = slice( z1, z2) ).values * units( 'km' )


# # Logic to determine if x-axis should be peak-relative or not
# #------------------------------------------------------------
# if( id_terrain == True ):
    
#     # Convert x-coordinates from grid-relative to peak-relative coordinates
#     xh = xh - peak_pos * units( 'km' )
    
#     # Change x-axis label to reflect coordinate shift
#     ax.set_xlabel( "Distance from Peak (km)" , fontsize = 18, fontweight = 'bold' )
    
#     # sx = int( sx - peak_pos )
    
#     # uh_x = int( uh_x - peak_pos )
    
# # Label the axis as grid-relative    
# else:
#     ax.set_xlabel( "Zonal Distance (km)" , fontsize = 18, fontweight = 'bold' )
# #------------------------------------------------------------


# #-----------------------------------------------------------------------------
# # End: Get variables from xarray dataset
# #-----------------------------------------------------------------------------


# #-----------------------------------------------------------------------------
# # Begin: Unit conversions required for MetPy functions
# #   Notes: metpy.quantify() adds native units from Dataset
# #   Notes: pint_xarray provides wrapper to convert units using pint.to()
# #-----------------------------------------------------------------------------
   
# # Get native units for required variables
# u = u.metpy.quantify()
# v = v.metpy.quantify()

# pressure = pressure.metpy.quantify()
# pressure = pressure.pint.to( 'hPa' )  

# #-----------------------------------------------------------------------------
# # End: Unit conversions required for MetPy functions
# #-----------------------------------------------------------------------------


# #-----------------------------------------------------------------------------
# # Begin: Calculations
# #-----------------------------------------------------------------------------

# # Compute Bunkers Storm Motions and Mean Wind for given inflow sounding 
# rm, lm, mean_uv = mpcalc.bunkers_storm_motion( 
#                                               pressure = pressure[:, sy, sx],
#                                               u = u[:, sy ,sx],
#                                               v = v[:, sy, sx],
#                                               height = zh
#                                               )

# # Convert zh into 3D array required for SRH function
# zh_grid = np.zeros( shape = (nz, ny, nx) )

# # Interpolate heights to a 3D domain required for SRH function
# #-------------------------------------------------------------
# for k in range( 0, nz ):
#     zh_grid[k,:,:] = zh[k].m
    
#     # Convert to meters
#     zh_grid[k,:,:] = zh_grid[k,:,:] * 1000.0
# #-------------------------------------------------------------

# # Compute 0-1 km SRH using Cython function (Must convert from pint to np.array using .values method)
# srh1km = srh(
#               u = u.values,              # units: m/s
#               v = v.values,              # units: m/s
#               zh = zh_grid,              # units: m
#               zs = zs.values,            # units: m
#               u_storm = rm[0].m,         # units: m/s
#               v_storm = rm[1].m,         # units: m/s
#               layer = 1000.0             # units: m
#             )

# #-----------------------------------------------------------------------------
# # End: Calculations
# #-----------------------------------------------------------------------------


# #-----------------------------------------------------------------------------
# # Begin: Plotting Data
# #-----------------------------------------------------------------------------

# # Convert winds to knots if desired (Must be after SRH call for proper calculation!)
# #----------------------------------
# if( knots == True ):
#     u = u.pint.to( 'knots' )
#     v = v.pint.to( 'knots' )
#     rm = rm.to( 'knots' )
#     wind_units = 'kts'
# else:
#     wind_units = 'ms$^{-1}$'
# #----------------------------------

# # Create a plotting grid based on the x,y coordinates
# X,Y = np.meshgrid( xh, yh )

# # Add a grid to the axis
# # ax.grid()

# # Begin: Terrain Logic
# #---------------------

# # If simulation includes terrain
# if( terrain == True ):
    
#     # Mask out lowest level of terrain (Needed for Realistic Terrain, Does not affect Idealized Terrain)
#     trn_mask = np.ma.array( zs, mask = zs < 0.01 )
    
#     # Create filled contours of the terrain field
#     cf_trn = ax.contourf(
#                           X,Y,
#                           trn_mask,
#                           levels = np.arange( 0, peak_elv + cint, cint ),
#                           cmap = plt.cm.get_cmap( 'copper' ),#.reversed(),
#                           alpha = 0.35
#                         )
    
#     # # Generate a reflectivity colorbar
#     # cbar = fig.colorbar ( cf_trn, pad = 0.05 )
    
#     # # Modify the colorbar label
#     # cbar.set_label( label = 'Elevation (m)', fontsize = 18, fontweight = 'bold' )
    
#     # Begin: Logic for terrain contour lines
#     #---------------------------------------
#     if( id_terrain == False ):
        
#         # Add contour lines for realistic terrain
#         c_trn = ax.contour(
#                             X,Y,
#                             trn_mask,
#                             levels = np.arange( 0, peak_elv + (cint*5), cint*5 ),
#                             linewidths = 0.5,
#                             alpha = 0.5
#                             )           
#     else:
        
#         # Add contour lines for idealized terrain
#         c_trn = ax.contour(
#                             X,Y,
#                             trn_mask,
#                             levels = np.arange( 0, peak_elv + cint, cint ),
#                             colors = 'saddlebrown',
#                             linewidths = 1.0,
#                             alpha = 1.0
#                           )

#     # End: Logic for terrain contour lines
#     #---------------------------------------
        
# # End: Terrain Logic
# #---------------------
        

# # Add UH swath contours    
# cs4 = ax.contourf(
#                   X,Y,
#                   scipy.ndimage.gaussian_filter( input = shs, sigma = smooth ), 
#                   levels = np.arange( 200, 8200, 4000 ), 
#                   colors = "black",
#                   alpha = 0.35
#                  )

# # # Create legend label for UH contour
# cs4_leg = mpatches.Patch( color = 'k', alpha = 0.35, label = 'UH Swath > 200 m$^2$s$^{-2}$' )


# # Slicing index used for clean wind vector plotting
# skip = ( slice( None, None, skip_val ), slice( None, None, skip_val ) )

# # Add SR wind vectors at surface
# q1 = ax.barbs( 
#               X[skip], Y[skip],
#               u[0,:,:][skip] - rm[0], v[0,:,:][skip] - rm[1],
#               length = 6.0, color = 'black', alpha = 0.75, label = 'Surface SR Wind {}'.format( wind_units )
#             )

# # Add SR wind vectors at 3km
# q2 = ax.barbs( 
#               X[skip], Y[skip],
#               u[km3,:,:][skip] - rm[0], v[km3,:,:][skip] - rm[1],
#               length = 6.0, color = 'blue', alpha = 0.75, label = '3 km SR Wind {}'.format( wind_units )
#             )

# # Add SR wind vectors at 5km
# q3 = ax.barbs( 
#               X[skip], Y[skip],
#               u[km5,:,:][skip] - rm[0], v[km5,:,:][skip] - rm[1],
#               length = 6.0, color = 'tab:red', alpha = 0.75, label = '5 km SR Wind {}'.format( wind_units )
#             )

# # Add 0-1 km SRH contours
# cs1 = ax.contour(
#                   X,Y,
#                   scipy.ndimage.gaussian_filter( input = srh1km, sigma = smooth ),
#                   levels = np.arange( 100, 600, 100 ), 
#                   linestyles = "-" , linewidths = 1.5, 
#                   colors = [ 'indigo' ], alpha = 0.75
#                 )

# # Label the 0-1 km SRH contours
# cs1_lab = ax.clabel(
#                     cs1,
#                     levels = np.arange( 100, 600, 100 ),
#                     fmt = '%3.0f',
#                     fontsize = 8,
#                     inline = True,
#                     colors = [ 'indigo' ],
#                     use_clabeltext = True
#                     )

# # Add CIN as hatched filled contour
# cf3 = ax.contourf(
#                   X,Y, 
#                   scipy.ndimage.gaussian_filter( input = cin, sigma = smooth ),
#                   levels = np.arange( 50, 550, 250 ), 
#                   hatches = [ '/' ], colors = 'tab:cyan',
#                   alpha = 0.40
#                  )

# # Add CIN contour line                
# cs3 = ax.contour(
#                   X,Y, 
#                   scipy.ndimage.gaussian_filter( input = cin, sigma = smooth ),
#                   levels = np.arange( 50, 250, 200 ), 
#                   linestyles = "--" , linewidths = 2.0, 
#                   colors = [ 'b' ], alpha = 0.75
#                 )

# # Add CAPE contours
# cs2 = ax.contour(
#                   X,Y, 
#                   scipy.ndimage.gaussian_filter( input = cape, sigma = smooth ),
#                   levels = np.arange( 500, 3000, 500 ), 
#                   linestyles = "--" , linewidths = 2.0, 
#                   colors = [ 'crimson' ], alpha = 0.75
#                 )

# # Label CAPE contours
# cs2_lab = ax.clabel(
#                     cs2,
#                     levels = np.arange( 500, 3000, 500 ),
#                     fmt = '%3.0f',
#                     fontsize = 8,
#                     inline = True,
#                     colors = [ 'crimson' ],
#                     use_clabeltext = True
#                    )

# # Create a halo around contour labels
# plt.setp(
#           [cs1_lab, cs2_lab], path_effects =[ 
#           PathEffects.withStroke( linewidth = 3, foreground = 'w', alpha = 0.75 ) ] 
#         )

# # Mask out reflectivity values below 10 dBZ
# dbz_mask = np.ma.array( dbz , mask = dbz <= 10.0 )

# # Add reflectivity as a filled contour
# cf = ax.contourf(
#                   X,Y,
#                   dbz_mask, 
#                   levels = np.arange( start = -10, stop = 80, step = 5.0 ),
#                   cmap = plots.ctables.registry.get_colortable( "NWSReflectivity" ),
#                   alpha = 0.45,
#                   norm = Normalize(-10, 80)
#                 )

# # Create an SRH legend label
# cs1_leg = mlines.Line2D( 
#                         [], [], color = 'indigo',
#                         linestyle = '-', linewidth = 2.0,
#                         label = '0-1 km SRH > 100 m$^2$s$^{-2}$'
#                         )
# # Create an CAPE legend label
# cs2_leg = mlines.Line2D( 
#                         [], [], color = 'crimson',
#                         linestyle = '--', linewidth = 2.0,
#                         label = 'CAPE > 500 J$kg^{-1}$'
#                         )

# # Create an CIN legend label
# cs3_leg = mpatches.Patch( facecolor = 'tab:cyan', edgecolor = 'k', hatch = '/', 
#                           alpha = 0.75, label = 'CIN > 50 J$kg^{-1}$' )

# # Label inflow sounding point
# snd_xy = ax.plot(
#                   sx, sy,
#                   marker = '*', markersize = 20,
#                   color = 'yellow', markeredgecolor = 'k',
#                   label = 'Inflow Sounding'
#                 )

# # # Label UH maxima point
# # uh_xy = ax.plot(
# #                 uh_x, uh_y,
# #                 marker = 'X', markersize = 10, color = 'aqua', 
# #                 markeredgecolor = 'k', label = 'UH Maxima'
# #                 )

# # Set y-axis label
# ax.set_ylabel( "Meridional Distance (km)" , fontsize = 18, fontweight = 'bold' )

# # Right Hand Plot Title (Model Desc)
# plt.title( 
#           "\nEnvironment Mesoanalysis",
#           loc = "left", fontsize = 14, fontweight = 'bold' 
#           )       

# # Right Hand Plot Title (Model Desc)
# plt.title( 
#           "1.0 km Radar Reflectivity > 10 dBZ \n {}-Sigma Gaussian Smoothed Parameters".format( smooth ),
#           loc = "right", fontsize = 12, fontweight = 'bold' 
#           )  


# # Create Legend
# ax.legend( 
#           handles = [ cs1_leg, cs2_leg, cs3_leg, cs4_leg,  q1, q2, q3, snd_xy[0] ], #, uh_xy[0] ],
#           loc = 1, facecolor = 'white', 
#           framealpha = 0.75, fontsize = 14
#          )

# #  # Save the current RH plot
# # fig.savefig(
# #             fname = 'test_meso.jpeg',
# #             dpi = 300,
# #             bbox_inches = "tight"
# #            )

# #-----------------------------------------------------------------------------
# # End: Plotting Data
# #-----------------------------------------------------------------------------

