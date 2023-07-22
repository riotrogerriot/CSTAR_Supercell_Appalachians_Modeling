#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 23:04:23 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Function: cm1_cross_section_plotter
#
#   Update Records:
#       (3/23/22) Script Created
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


import pint
import pint_xarray
import cf_xarray


#-----------------------------------------------------------------------------
# Begin: Define cm1_cross_section_plotter
#-----------------------------------------------------------------------------
def cm1_cross_section_plotter( 
                              ds_interp, fig, ax, x1, x2, z1, z2, nx, nz, uh_y,
                              sx, sy, rm_x, rm_y, shortname, current_model_time,
                              skip_val = 25, terrain = False, knots = False,
                              id_terrain = False, peak_pos = 350.0
                             ):

    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the u-component of the wind (m/s)    
    u = ds_interp.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get the v-component of the wind (m/s)   
    v = ds_interp.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get the w-component of the wind (m/s)   
    w = ds_interp.metpy.parse_cf( 'winterp' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get pressure (Pa)
    pressure = ds_interp.metpy.parse_cf( 'prs' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get the radar reflectivity (dBZ)
    dbz = ds_interp.metpy.parse_cf( 'dbz' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get Vertical Vorticity (1/s)
    zvort = ds_interp.metpy.parse_cf( 'zvort' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get potential temperature perturbation (K)
    thpert = ds_interp.metpy.parse_cf( 'thpert' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get potential temperature perturbation (K)
    th = ds_interp.metpy.parse_cf( 'th' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get water-vapor mixing ratio (kg/kg)
    qv = ds_interp.metpy.parse_cf( 'qv' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get liquid cloud mixing ratio (kg/kg)
    qc = ds_interp.metpy.parse_cf( 'qc' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    # Get liquid water-vapor mixing ratio (kg/kg)
    qi = ds_interp.metpy.parse_cf( 'qi' ).isel( time = 0, zh = slice( 0, nz ), yh = uh_y, xh = slice( 0, nx ) )
    
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get the x-axis coordinates (km)
    xh = ds_interp.coords[ 'xh' ].isel( xh = slice( 0, nx) ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)
    yh = ds_interp.coords[ 'yh' ].isel( yh = uh_y ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)      
    zh = ds_interp.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' ) 
    
    
    # Logic to determine if x-axis should be peak-relative or not
    #------------------------------------------------------------
    if( id_terrain == True ):
        
        # Convert x-coordinates from grid-relative to peak-relative coordinates
        xh = xh - peak_pos * units( 'km' )
        
        # Change x-axis label to reflect coordinate shift
        ax.set_xlabel( "Distance from Peak (km)" , fontsize = 18, fontweight = 'bold' )
        
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
    pressure = pressure.metpy.quantify()
    qv = qv.metpy.quantify()
    
    # Convert to required units for Air Temp & Td calculations
    pressure = pressure.pint.to( 'hPa' )  
    qv = qv.pint.to( 'g/kg' )
    
    # Convert to required units for Air Temp & Td calculations
    pressure = pressure.pint.to( 'hPa' )  
    
    #-----------------------------------------------------------------------------
    # End: Unit conversions required for MetPy functions
    #-----------------------------------------------------------------------------
    
    
    
    #-----------------------------------------------------------------------------
    # Begin: Calculations
    #-----------------------------------------------------------------------------
    
    # Total Cloud Mixing-Ratio (kg/kg)
    qci = qc + qi
    
    # Convert potential temperature to air temperature (K)
    temp = mpcalc.temperature_from_potential_temperature( pressure, th ) 
    
    # Get native units
    temp = temp.metpy.quantify()
    
    # Convert from K to degC
    temp = temp.pint.to( 'degC' )
    
    # Compute the water vapor pressure using pressure and mixing-ratio
    e  = mpcalc.vapor_pressure( pressure, qv )
    
    # Compute dew point from vapor pressure (degC)
    td = mpcalc.dewpoint( e ) 
    
    # Compute Relative Humidity from Temp and Dew Point
    rh = mpcalc.relative_humidity_from_dewpoint( temp, td )
    
    # Mask out missing values below the terrain using missing data code (ncdump)
    # u_mask = np.ma.array( u.values, mask = u == -999999.875 )
    # v_mask = np.ma.array( v.values, mask = v == -999999.875 )
    w_mask = np.ma.array( w, mask = u == -999999.875 )
    qci_mask = np.ma.array( qci, mask = qci == -999999.875 )
    thpert_mask = np.ma.array( thpert, mask = thpert == -999999.875 )
    rh_mask = np.ma.array( rh, mask = rh == -999999.875 )
    zvort_mask = np.ma.array( zvort, mask = zvort == -999999.875 )
    
    # Mask out values of reflectivity below given threshold
    dbz_mask = np.ma.array( dbz , mask = dbz <= 30.0 )
    
    # Convert winds to knots if desired
    #----------------------------------
    if( knots == True ):
        u = u.pint.to( 'knots' )
        v = v.pint.to( 'knots' )
        rm_x = rm_x.to( 'knots' )
        rm_y = rm_y.to( 'knots' )
        wind_units = 'kts'
    else:
        wind_units = 'ms$^{-1}$'
    #----------------------------------
    
    # Calculate Helicity
    wzvort = w_mask * zvort_mask
    
    #-----------------------------------------------------------------------------
    # End: Calculations
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Plotting Data
    #-----------------------------------------------------------------------------
    
    # Create a plotting grid based on the x,y coordinates
    X,Z = np.meshgrid( xh, zh )
    
    
    # Plot data
    #-----------------------------------------------------------------------------
       
    # Contour relative humidity
    cf1 = ax.contourf(
                      X, Z, rh_mask, 
                      levels = np.arange( 0, 1.3, 0.1 ),
                      linestyles = 'solid', cmap = 'binary',
                      alpha = 0.8
                      )
    
    # Generate a relative humidity colorbar
    cbar = fig.colorbar ( cf1, pad = 0.05 )
    
    # Modify the colorbar label
    cbar.set_label( label = 'Relative Humidity', fontsize = 18, fontweight = 'bold' )
    
    
    # Contour helicity
    cf2 = ax.contour(
                      X,Z,wzvort, 
                      levels = np.arange( 0.1, 0.2, 0.1 ),
                      linestyles = 'dashed', linewidths = 2.5, 
                      colors = 'limegreen', alpha = 1.0
                     )
    
    # Contour cloud mixing-ratio
    cf3 = ax.contour(
                      X,Z,qci,
                      levels = np.arange( 0.0001, 0.0006, 0.0005 ),
                      linestyles = 'solid', linewidths = 3.0, 
                      colors = 'black', alpha = 1.0
                     )
    
    # Contour potential temperature perturbation (0-3 km only)
    cf4 = ax.contour(
                      X[0:13,:],Z[0:13,:],thpert_mask[0:13,:], 
                      levels = np.arange( -1, 0, 1 ),
                      linestyles = 'dotted', linewidths = 2.5,
                      colors = 'blue', alpha = 1.0
                     )
    
    # Fill Contour radar reflectivity
    cf = ax.contourf(
                      X,Z,
                      dbz_mask, 
                      levels = np.arange( start = -10, stop = 80, step = 5.0 ),
                      cmap = plots.ctables.registry.get_colortable( "NWSReflectivity" ),
                      alpha = 0.40, norm = Normalize(-10, 80),
                    )
    
    # Slicing index used for clean wind vector plotting
    skip = ( slice( None, None, 2 ), slice( None, None, skip_val ) )
    
    u = u.values
    v = v.values 
    
    u_mask = np.ma.array( u, mask = u == -1943844.25 )
    v_mask = np.ma.array( v, mask = v == -1943844.25 )
    
    # Add storm-relative wind vectors
    q = ax.barbs( 
                 X[skip], Z[skip],
                 u_mask[skip] - rm_x , v_mask[skip] - rm_y,
                 length = 5.0, alpha = 0.65
                )
    
    
    # Create legend
    #-----------------------------------------------------------------------------
    
    # Create legend label for helicity contour
    cs2_lab = mlines.Line2D( 
                            [], [], color = 'limegreen',
                            linestyle = 'dashed', linewidth = 2.5,
                            label = 'w\u03b6 > 0.1 ms$^{-2}$'
                           )
    
    # Create legend label for cloud mixing-ratio contour
    cs3_lab = mlines.Line2D( 
                            [], [], color = 'black',
                            linestyle = 'solid', linewidth = 3.0,
                            label = '0.5 gkg$^{-1}$ $r_{Cloud}$'
                           )
    
    # Create legend label for cloud mixing-ratio contour
    cs4_lab = mlines.Line2D( 
                            [], [], color = 'blue',
                            linestyle = 'dotted', linewidth = 2.5,
                            label = '0-3 km -1K \u03F4 Pert.'
                           )
    
    # Create Legend
    ax.legend( 
              handles = [ cs2_lab, cs3_lab, cs4_lab ],
              loc = 1, facecolor = 'white', 
              framealpha = 0.55, fontsize = 14
             )
    
    
    # Finalize Formatting
    #-----------------------------------------------------------------------------
    
    # Set y-axis label
    ax.set_ylabel( "Altitude (km)" , fontsize = 18, fontweight = 'bold' )
    
    # Right Hand Plot Title (Model Desc)
    plt.title( 
              "CM1-{}: t = {} min.\nc: Zonal Cross-Section: y = {} km ".format( shortname, current_model_time, round(yh.m,0) ),
              loc = "left", fontsize = 14, fontweight = 'bold' 
             )
    
    # Right Hand Plot Title (Model Desc)
    plt.title( 
              "Radar Reflectivity > 30 dBZ \nStorm-Relative Winds ({})".format( wind_units ),
              loc = "right", fontsize = 12, fontweight = 'bold' 
             )               
    
    # Set plot windows
    ax.set_xlim( xh[x1], xh[x2] )
    ax.set_ylim( zh[z1], zh[z2] )
    
    # Add a grid to the axis
    ax.grid()
    
    #-----------------------------------------------------------------------------
    # End: Plotting Data
    #-----------------------------------------------------------------------------
    
    # Return the plot
    return( ax )

#-----------------------------------------------------------------------------
# End: Define cm1_radar_plotter
#-----------------------------------------------------------------------------







# x1 = 850
# x2 = 1050
# uh_y = 880
# z1 = 0
# z2 = 38
# peak_pos = 400.0
# id_terrain = True
# knots = False
# terrain = True
# skip_val = 20

# filename = '/Users/roger/Desktop/test_cm1_analysis/cm1out_000036.nc'
# filename2 = '/Users/roger/Desktop/test_cm1_analysis/cm1out_000036_i.nc'

# # Open the current terrain-interpolated netCDF file with xarray 
# ds = xr.open_dataset( filename, engine = "netcdf4", decode_cf = True )

# # Open the current terrain-interpolated netCDF file with xarray 
# ds_interp = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )

# # Representative of the individual plots that compose the figure (Each letter is a separate subplot and . skips that position)
# figure_mosaic = """
#                  AA
#                  AA
#                 """

# # Create figure and axes objects
# fig, axes = plt.subplot_mosaic( 
#                                mosaic = figure_mosaic, 
#                                figsize = ( 11, 8.5 ), 
#                                constrained_layout = True
#                               )

#  # Construct axis objects for the radar plot
# ax= fig.add_subplot( axes['A'] )


# # Data Variables
# #-----------------------------------------------------------------------------

# # Get the u-component of the wind (m/s)    
# u = ds_interp.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh = slice( z1, z2), 
#                                          yh = uh_y, xh = slice( x1, x2) )

# # Get the v-component of the wind (m/s)   
# v = ds_interp.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh = slice( z1, z2), 
#                                          yh = uh_y, xh = slice( x1, x2) )

# # Get the w-component of the wind (m/s)   
# w = ds_interp.metpy.parse_cf( 'winterp' ).isel( time = 0, zh = slice( z1, z2), 
#                                          yh = uh_y, xh = slice( x1, x2) )

# # Get pressure (Pa)
# pressure = ds_interp.metpy.parse_cf( 'prs' ).isel( time = 0, zh = slice( z1, z2), 
#                                             yh = uh_y, xh = slice( x1, x2) )

# # Get the radar reflectivity (dBZ)
# dbz = ds_interp.metpy.parse_cf( 'dbz' ).isel( time = 0, zh = slice( z1, z2), 
#                                        yh = uh_y, xh = slice( x1, x2) )

# # Get Vertical Vorticity (1/s)
# zvort = ds_interp.metpy.parse_cf( 'zvort' ).isel( time = 0, zh = slice( z1, z2), 
#                                            yh = uh_y, xh = slice( x1, x2) )

# # Get potential temperature perturbation (K)
# thpert = ds_interp.metpy.parse_cf( 'thpert' ).isel( time = 0, zh = slice( z1, z2), 
#                                              yh = uh_y, xh = slice( x1, x2) )

# # Get potential temperature perturbation (K)
# th = ds_interp.metpy.parse_cf( 'th' ).isel( time = 0, zh = slice( z1, z2), 
#                                      yh = uh_y, xh = slice( x1, x2) )

# # Get water-vapor mixing ratio (kg/kg)
# qv = ds_interp.metpy.parse_cf( 'qv' ).isel( time = 0, zh = slice( z1, z2), 
#                                      yh = uh_y, xh = slice( x1, x2) )

# # Get liquid cloud mixing ratio (kg/kg)
# qc = ds.metpy.parse_cf( 'qc' ).isel( time = 0, zh = slice( z1, z2), 
#                                      yh = uh_y, xh = slice( x1, x2) )

# # Get liquid water-vapor mixing ratio (kg/kg)
# qi = ds.metpy.parse_cf( 'qi' ).isel( time = 0, zh = slice( z1, z2), 
#                                      yh = uh_y, xh = slice( x1, x2) )

# # # Get the terrain field if simulation includes terrain (m)
# # if( terrain == True ):
# #     zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0, zh = slice( z1, z2), 
# #                                          yh = uh_y, xh = slice( x1, x2) )


# # Coordinate Variables
# #-----------------------------------------------------------------------------

# # Get the x-axis coordinates (km)
# xh = ds_interp.coords[ 'xh' ].isel( xh = slice( x1, x2) ).values * units( 'km' )

# # Get the y-axis coordinates (km)
# yh = ds_interp.coords[ 'yh' ].isel( yh = uh_y ).values * units( 'km' )

# # Get the y-axis coordinates (km)      
# zh = ds_interp.coords[ 'zh' ].isel( zh = slice( z1, z2) ).values * units( 'km' )


# # Logic to determine if x-axis should be peak-relative or not
# #------------------------------------------------------------
# if( id_terrain == True ):
    
#     # Convert x-coordinates from grid-relative to peak-relative coordinates
#     xh = xh - peak_pos * units( 'km' )
    
#     # Change x-axis label to reflect coordinate shift
#     ax.set_xlabel( "Distance from Peak (km)" , fontsize = 18, fontweight = 'bold' )
    
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
# qv = qv.metpy.quantify()

# # Convert to required units for Air Temp & Td calculations
# pressure = pressure.pint.to( 'hPa' )  
# qv = qv.pint.to( 'g/kg' )

# # Convert to required units for Air Temp & Td calculations
# pressure = pressure.pint.to( 'hPa' )  

# #-----------------------------------------------------------------------------
# # End: Unit conversions required for MetPy functions
# #-----------------------------------------------------------------------------



# #-----------------------------------------------------------------------------
# # Begin: Calculations
# #-----------------------------------------------------------------------------

# # Total Cloud Mixing-Ratio (kg/kg)
# qci = qc + qi

# # Convert potential temperature to air temperature (K)
# temp = mpcalc.temperature_from_potential_temperature( pressure, th ) 

# # Get native units
# temp = temp.metpy.quantify()

# # Convert from K to degC
# temp = temp.pint.to( 'degC' )

# # Compute the water vapor pressure using pressure and mixing-ratio
# e  = mpcalc.vapor_pressure( pressure, qv )

# # Compute dew point from vapor pressure (degC)
# td = mpcalc.dewpoint( e ) 

# # Compute Relative Humidity from Temp and Dew Point
# rh = mpcalc.relative_humidity_from_dewpoint( temp, td )

# # Mask out missing values below the terrain using missing data code (ncdump)
# u_mask = np.ma.array( u.values, mask = u == -1e+06 )
# v_mask = np.ma.array( v.values, mask = v == -1e+06 )
# w_mask = np.ma.array( w, mask = u == -999999.875 )
# qci_mask = np.ma.array( qci, mask = qci == -999999.875 )
# thpert_mask = np.ma.array( thpert, mask = thpert == -999999.875 )
# rh_mask = np.ma.array( rh, mask = rh == -999999.875 )
# zvort_mask = np.ma.array( zvort, mask = zvort == -999999.875 )

# # Mask out values of reflectivity below given threshold
# dbz_mask = np.ma.array( dbz , mask = dbz <= 40.0 )

# # Convert winds to knots if desired
# #----------------------------------
# if( knots == True ):
#     u = u.pint.to( 'knots' )
#     v = v.pint.to( 'knots' )
#     wind_units = 'kts'
# else:
#     wind_units = 'ms$^{-1}$'
# #----------------------------------

# # Calculate Helicity
# wzvort = w_mask * zvort_mask

# #-----------------------------------------------------------------------------
# # End: Calculations
# #-----------------------------------------------------------------------------


# #-----------------------------------------------------------------------------
# # Begin: Plotting Data
# #-----------------------------------------------------------------------------

# # Create a plotting grid based on the x,y coordinates
# X,Z = np.meshgrid( xh, zh )


# # Plot data
# #-----------------------------------------------------------------------------
       
# # Contour relative humidity
# cf1 = ax.contourf(
#                   X,Z,rh_mask, 
#                   levels = np.arange( 0, 1.3, 0.1 ),
#                   linestyles = 'solid', cmap = 'binary',
#                   alpha = 0.65
#                  )

# # Generate a relative humidity colorbar
# cbar = fig.colorbar ( cf1, pad = 0.05 )

# # Modify the colorbar label
# cbar.set_label( label = 'Relative Humidity', fontsize = 18, fontweight = 'bold' )


# # Contour helicity
# cf2 = ax.contourf(
#                   X,Z,wzvort, 
#                   levels = np.arange( 0.1, 2.1, 1.0 ),
#                   hatches = [ '/' ], colors = None,
#                   alpha = 0.50
#                  )

# cf2 = ax.contour(
#                   X,Z,wzvort, 
#                   levels = np.arange( 0.1, 2.1, 1.0 ),
#                   linestyles = '-', linewidths = 1.0, 
#                   colors = 'k',
#                   alpha = 0.50
#                  )

# # Contour cloud mixing-ratio
# cf3 = ax.contour(
#                   X,Z,qci,
#                   levels = np.arange( 0.0001, 0.0006, 0.0005 ),
#                   linestyles = 'solid', linewidths = 2.0, 
#                   colors = 'black', alpha = 1.0
#                  )

# # Contour potential temperature perturbation (0-3 km only)
# cf4 = ax.contour(
#                   X[0:13,:],Z[0:13,:],thpert_mask[0:13,:], 
#                   levels = np.arange( -1, 0, 1 ),
#                   linestyles = '--', linewidths = 3.0,
#                   colors = [ 'blue' ], alpha = 1.0
#                  )

# # Fill Contour radar reflectivity
# cf = ax.contourf(
#                   X,Z,
#                   dbz_mask, 
#                   levels = np.arange( start = -10, stop = 80, step = 5.0 ),
#                   cmap = plots.ctables.registry.get_colortable( "NWSReflectivity" ),
#                   alpha = 0.40, norm = Normalize(-10, 80),
#                 )

# # Slicing index used for clean wind vector plotting
# skip = ( slice( None, None, 2 ), slice( None, None, skip_val ) )

# # Add wind vectors
# q = ax.barbs( 
#              X[skip], Z[skip],
#              u_mask[skip] , v_mask[skip],
#              length = 5.0, alpha = 0.65
#             )


# # Create legend
# #-----------------------------------------------------------------------------

# # Create legend label for helicity contour
# cs2_lab = mpatches.Patch( facecolor = None , edgecolor = 'k', hatch = '/', 
#                           alpha = 0.5, label = 'w\u03b6 > 0.1 ms$^{-2}$' )

# # Create legend label for cloud mixing-ratio contour
# cs3_lab = mlines.Line2D( 
#                         [], [], color = 'black',
#                         linestyle = '-', linewidth = 2.0,
#                         label = '0.5 gkg$^{-1}$ $r_{Cloud}$'
#                        )

# # # Create legend label for gust front contour
# # cs4a_lab = mlines.Line2D( 
# #                         [], [], color = 'royalblue',
# #                         linestyle = '--', linewidth = 3.0,
# #                         label = '0-3 km -3 K \u03F4 Pert.'
# #                        )

# # Create legend label for gust front contour
# cs4b_lab = mlines.Line2D( 
#                         [], [], color = 'blue',
#                         linestyle = '--', linewidth = 3.0,
#                         label = '0-3 km -1 K \u03F4 Pert.'
#                        )

# # Create Legend
# ax.legend( 
#           handles = [ cs2_lab, cs3_lab, cs4b_lab ],
#           loc = 1, facecolor = 'white', 
#           framealpha = 0.75, fontsize = 14
#          )


# # Finalize Formatting
# #-----------------------------------------------------------------------------

# # Set y-axis label
# ax.set_ylabel( "Altitude (km)" , fontsize = 18, fontweight = 'bold' )

# # Right Hand Plot Title (Model Desc)
# plt.title( 
#           "\nZonal Cross-Section Analysis \ny = {} km ".format( round(yh.m,0) ),
#           loc = "left", fontsize = 14, fontweight = 'bold' 
#          )

# # Right Hand Plot Title (Model Desc)
# plt.title( 
#           "Radar Reflectivity > 40 dBZ \nStorm-Relative Winds ({})".format( wind_units ),
#           loc = "right", fontsize = 12, fontweight = 'bold' 
#          )               

# # Add a grid to the axis
# ax.grid()

# #-----------------------------------------------------------------------------
# # End: Plotting Data
# #-----------------------------------------------------------------------------



