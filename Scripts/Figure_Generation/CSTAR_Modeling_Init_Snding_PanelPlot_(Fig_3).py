#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:35:32 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Load external libraries
#-----------------------------------------------------------------------------
import os
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.plots as plots
import metpy.calc as mpcalc
import pint
import pint_xarray
import cf_xarray
from datetime import datetime

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Begin Main Script
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Record Script start time
startTime = datetime.now()

# Program name
program = 'cm1_init_conditions_figure.py'

# Report program status to terminal
print( "\nBegin {}...".format( program ) )


#-----------------------------------------------------------------------------
# Begin: Plot Configuration
#-----------------------------------------------------------------------------

# (!!! Requires output files that were to large to upload to repo please contact corresponding author for data if needed !!!)
output_dir = '/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR/Initial_Conditions_Output/Panel_Plot/'

sound_type = [ 'Upstream', 'Peak', 'Downstream' ]
sim_type = [ 
            'CS', 'CS', 'CS', 'NC', 'NC', 'NC',
            'CS', 'CS', 'CS', 'NC', 'NC', 'NC',
            'CS', 'CS', 'CS', 'NC', 'NC', 'NC',
           ]

subset_type = [ '_All', '_All', '_All', '_All', '_All', '_All',
                '_Out', '_Out', '_Out', '_Out', '_Out', '_Out',
                '_No_Out', '_No_Out', '_No_Out', '_No_Out', '_No_Out', '_No_Out', ]

bss_time = [ 't = 0', 't = BSS1', 't = BSS2', 't = 0', 't = BSS1', 't = BSS2',
             't = 0', 't = BSS1', 't = BSS2', 't = 0', 't = BSS1', 't = BSS2',
             't = 0', 't = BSS1', 't = BSS2', 't = 0', 't = BSS1', 't = BSS2' ]


filenames = [ 'CS_Upstream_All', 'CS_Peak_All', 'CS_Downstream_All',
              'NC_Upstream_All', 'NC_Peak_All', 'NC_Downstream_All',
              'CS_Upstream_Out', 'CS_Peak_Out', 'CS_Downstream_OUT',
              'NC_Upstream_Out', 'NC_Peak_Out', 'NC_Downstream_OUT',
              'CS_Upstream_No_Out', 'CS_Peak_No_Out', 'CS_Downstream_No_Out',
              'NC_Upstream_No_Out', 'NC_Peak_No_Out', 'NC_Downstream_No_Out'
            ]


# Set tickmark size
plt.rcdefaults()
plt.rc( 'font', weight = 'bold' )
plt.rc( 'xtick', labelsize = 14 )
plt.rc( 'ytick', labelsize = 14 )    


# Representative of the individual plots that compose the figure (Each letter is a separate subplot and . skips that position)
figure_mosaic = """
                ABC
                DEF
                GHI
                JKL
                MNO
                PQR
                """ 

# Create figure and axes objects
fig, axes = plt.subplot_mosaic( 
                               mosaic = figure_mosaic, 
                               figsize = ( 18, 20 ), 
                               tight_layout = True,
                               subplot_kw = { 'projection': 'skewx' }
                              )

# List for looping through each axis
current_ax = [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
               'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R' ]

#------------------------------------------------------------------------------
# Begin: Loop to plot on each axis
#------------------------------------------------------------------------------
for i in range( 0, len(axes) ):
    
    # Call the current axes object
    ax = fig.add_subplot( axes[ current_ax[i] ] )
    
    # Create a SkewT plotting object
    skew = plots.SkewT( fig, rotation = 45, subplot = ax, aspect = 100 )
    
    # Custom axes limits
    skew.ax.set_ylim( 1000, 100 )
    skew.ax.set_xlim( -70, 35 ) 
    
    # Label every other tickmark
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible( False )
    for label in ax.yaxis.get_ticklabels()[::2]:
        label.set_visible( False )

    if( current_ax[i] == 'A' or current_ax[i] == 'B' or current_ax[i] == 'C' ):    
        ax.xaxis.set_ticklabels([])
    if( current_ax[i] == 'G' or current_ax[i] == 'H' or current_ax[i] == 'I' ):
        ax.xaxis.set_ticklabels([])
    if( current_ax[i] == 'M' or current_ax[i] == 'N' or current_ax[i] == 'O' ):    
        ax.xaxis.set_ticklabels([])
        
        
    if( current_ax[i] == 'B' or current_ax[i] == 'C' or current_ax[i] == 'E' or current_ax[i] == 'F' ):
        ax.yaxis.set_ticklabels([])
    if( current_ax[i] == 'H' or current_ax[i] == 'I' or current_ax[i] == 'K' or current_ax[i] == 'L' ):
        ax.yaxis.set_ticklabels([])
    if( current_ax[i] == 'N' or current_ax[i] == 'O' or current_ax[i] == 'Q' or current_ax[i] == 'R' ):
        ax.yaxis.set_ticklabels([])
    
    # # Create a SkewT plotting object
    # skew = plots.SkewT( fig, rotation = 45, subplot = ax, aspect = 100 )
    
    # # Custom axes limits
    # skew.ax.set_ylim( 1000, 100 )
    # skew.ax.set_xlim( -70, 35 ) 
    
    # # Label every other tickmark
    # for label in ax.xaxis.get_ticklabels()[::2]:
    #     label.set_visible( False )
    # for label in ax.yaxis.get_ticklabels()[::2]:
    #     label.set_visible( False )
    
    # Add fiducial lines
    skew.plot_dry_adiabats( alpha = 0.35 )
    skew.plot_moist_adiabats( alpha = 0.35 )
    skew.plot_mixing_lines( alpha = 0.35 )
    skew.ax.axvline( 0, color = 'c', linestyle = '--', linewidth = 2.5, alpha = 0.35 )
    
    # Set the x-axis label
    if( i == 15 or i == 16 or i == 17 ):
        ax.set_xlabel( 
                      xlabel = 'Temperature (\u00B0C)',
                      fontsize = 14, fontweight = 'bold'
                     )
    
    # Set the y-axis label
    if( i == 0 or i == 3 or i == 6 or i == 9 or i == 12 or i == 15 ):
        ax.set_ylabel( 
                      ylabel = 'Pressure (hPa)',
                      fontsize = 14, fontweight = 'bold'
                     )
        
    # Place a hodograph into the current axis
    hax = inset_axes( skew.ax, '35%', '35%', loc= "upper left" )
    
    # Panel title logic
    if( i == 0 ):
        ax.set_title( 
                      sound_type[0] + '  ALL\n', loc = 'center',
                      fontsize = 20, fontweight = 'bold'
                    )
        ax.set_title( 
                      current_ax[i].lower() + ': ' + sim_type[i] + subset_type[i].upper(), loc = 'left',
                      fontsize = 14, fontweight = 'bold'
                    )
    elif( i == 1 ):
        ax.set_title( 
                      sound_type[1] + '  OUT\n', loc = 'center',
                      fontsize = 20, fontweight = 'bold'
                    )
        ax.set_title( 
                      current_ax[i].lower() + ': ' + sim_type[i] + subset_type[i].upper(), loc = 'left',
                      fontsize = 14, fontweight = 'bold'
                    )
    elif( i == 2 ):
        ax.set_title( 
                      sound_type[2] + ' NO_OUT\n', loc = 'center',
                      fontsize = 20, fontweight = 'bold'
                    )
        ax.set_title( 
                      current_ax[i].lower() + ': ' + sim_type[i] + subset_type[i].upper(), loc = 'left',
                      fontsize = 14, fontweight = 'bold'
                    )
    else:    
        ax.set_title( 
                      current_ax[i].lower() + ': ' + sim_type[i] + subset_type[i].upper(), loc = 'left',
                      fontsize = 14, fontweight = 'bold'
                    )

    # Add BSS time to right label
    ax.set_title( 
                  bss_time[i], loc = 'right',
                  fontsize = 14, fontweight = 'bold'
                )


    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Pull sounding from center of domain
    sx = 300
    sy = 200
    
    # Open the current file
    fname = output_dir + 'cm1out_' + filenames[i] + '.nc'
        
    # Report to terminal
    print( '\t{}: tNow opening {}'.format( current_ax[i], fname ) )
        
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
    # print( '\tNow closing {}'.format( fname ) )
    
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
    skew.plot( prs, temp, 'tab:red', linewidth = 2.5, label = 'T ' )
    skew.plot( prs, td, 'tab:green', linewidth = 2.5, label = 'Td' )
    
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
    
    # # Plot LCL, LFC, and EL levels
    # # skew.plot( lcl_pressure, lcl_temperature, marker = '_', markersize = 10, color = 'b' )
    # # skew.plot( lfc_pressure, lfc_temperature, marker = '_', markersize = 10, color = 'r' )
    # # skew.plot( el_pressure, el_temperature, marker = '_', markersize = 10, color = 'k' )
    
    # # Label LCL, LFC, and EL levels
    # # plt.text( lcl_temperature.to( 'degC').m + 5, lcl_pressure.m + 25 , 'LCL', color = 'b', size = 18 )
    # # plt.text( lfc_temperature.to( 'degC').m + 5, lfc_pressure.m + 20 , 'LFC', color = 'r', size = 18 )
    # # plt.text( el_temperature.to( 'degC').m + 5, el_pressure.m + 5, 'EL', color = 'k', size = 18 )
    
    # Shade areas of CAPE and CIN
    skew.shade_cin( prs, temp, prof, alpha = 0.2, label = 'CIN' )
    skew.shade_cape( prs, temp, prof, alpha = 0.2, label = 'CAPE' )

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
    
    # # Modify hodograph if requested
    # #----------------------------------------------
    # # if( mod_hodo == True ):
        
    # #     mod_angle = units.Quantity( 25, 'degrees' )   
        
    # #     # Get magnitude and direction from original components
    # #     wspd = mpcalc.wind_speed( u, v )
    # #     wdir = mpcalc.wind_direction( u , v )
        
    # #     # Modify the direction of the hodograph evenly throughout all layers
    # #     for j in range( 0, len( wdir) ):
    # #         wdir[j] = wdir[j] + mod_angle
        
    # #     # Compute the new wind components
    # #     u2,v2 = mpcalc.wind_components( wspd, wdir )
        
    # # # Plot the hodograph (m/s) with a color map based on height
    # # hodo2 = h.plot( 
    # #               u = u2,
    # #               v = v2,
    # #               # color = hodo_col[i],
    # #               # linestyle = line[0],
    # #               # alpha = 0.45,
    # #              )
    
    # Set hodograph plot limits
    hax.set_xlim( -5, 45 )
    hax.set_ylim( -10, 35 )    
    
    # # End Modify Hodograph
    # #---------------------------------------------

    # Compute MLCAPE/CIN
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin( prs, temp, td )
    
    # Compute Bunkers Storm Motions and Mean Wind
    rm, lm, mean_uv = mpcalc.bunkers_storm_motion(
                                                  pressure = prs,
                                                  u = u,
                                                  v = v,
                                                  height = zh
                                                  )

    # # Plot the mean wind as a vector
    # h.wind_vectors( u = mean_uv[0], v = mean_uv[1], color = 'dimgrey', scale = 1, alpha = 0.75, width = 1, headwidth = 3 )
    
    # # Plot the Bunker's Right & Left Motions as labeled pointd
    # txt5 = plt.text( rm[0].m, rm[1].m, 'RM', size = 8 )
    # txt6 = plt.text( lm[0].m, lm[1].m, 'LM', size = 8 )
    
    # # Place circles around RM & LM labels
    # plt.scatter( rm[0].m+2.5, rm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
    # plt.scatter( lm[0].m+2.5, lm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
    

    # Compute 0-3 km SRH
    srh_3km = mpcalc.storm_relative_helicity(
                                              u = u,
                                              v = v,
                                              height = zh,
                                              depth = units.Quantity(3000, 'm'),
                                              storm_u = rm[0],
                                              storm_v = rm[1]
                                            )
    
    srh_1km = mpcalc.storm_relative_helicity(
                                              u = u,
                                              v = v,
                                              height = zh,
                                              depth = units.Quantity(1000, 'm'),
                                              storm_u = rm[0],
                                              storm_v = rm[1]
                                            )

    p2 = np.arange( -105, 20, 14)
    
    txt1 = plt.text( 0.0, p2[3], 'CAPE: ' + str( round( int(mlcape.m), 0 ) ), fontsize = 14 )
    txt2 = plt.text( 0.0, p2[2], 'CIN: '  + str( round( int(mlcin.m),  0 ) ), fontsize = 14 )   
    txt3 = plt.text( 0.0, p2[1], 'SRH: ' + str( round( int(srh_3km[0].m), 0 ) ), fontsize = 14 )
    # txt4 = plt.text( 0.0, p2[0], 'SRH1KM: ' + str( round( float(srh_1km[0].m), 0 ) ) + ' m$^2$s$^{-2}$', fontsize = 12 )
    
    # Create a mask for wind barb plotting
    mask = prs >= 100.0 * units.hPa
    
    # Determine which units of wind should be plotted
    skew.plot_barbs( prs[ mask ][::2], u[ mask ][::2], v[ mask ][::2] )     

    # # Add legend 
    # if( i == 15 ):
    #     skew.ax.legend( loc = 3, facecolor = 'grey',  framealpha = 0.2, fontsize = 10 )
    
#------------------------------------------------------------------------------
# End: Loop to plot on each axis
#------------------------------------------------------------------------------







# Save the current skew-t
plt.savefig(
            fname = "init_snding_panelplot.jpeg",
            dpi = 300,
            bbox_inches = "tight"
           )
