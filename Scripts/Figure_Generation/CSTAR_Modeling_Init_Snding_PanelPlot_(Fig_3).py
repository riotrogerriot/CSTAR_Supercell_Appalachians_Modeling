#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 22:19:09 2022

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

output_dir = '/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR_Modeling_Project/Simulations/Initial_Conditions/Initial_Conditions_Output_Files/'
sim_type = [ 'Crosser', 'Non-Crosser' ]
sims = [ 'CS', 'NC']
sound_type = [ 'Upstream', 'Peak', 'Downstream', 'Upstream', 'Peak', 'Downstream' ]
bss_time = [ 't = 0', 't = BSS1', 't = BSS2', 't = 0', 't = BSS1', 't = BSS2' ]

# Representative of the individual plots that compose the figure (Each letter is a separate subplot and . skips that position)
figure_mosaic = """
                ABC
                DEF
                """ 
                
# Set tickmark size
plt.rcdefaults()
plt.rc( 'font', weight = 'bold' )
plt.rc( 'xtick', labelsize = 14 )
plt.rc( 'ytick', labelsize = 14 ) 

# Create figure and axes objects
fig, axes = plt.subplot_mosaic( 
                               mosaic = figure_mosaic, 
                               figsize = ( 30, 15 ), 
                               tight_layout = False,
                               subplot_kw = { 'projection': 'skewx' }
                              )

# List for looping through each axis
current_ax = [ 'A', 'B', 'C', 'D', 'E', 'F' ]

mod_hodo = True

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
    skew.ax.set_xlim( -75, 35 ) 
    
    # Label every other tickmark
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible( False )
    for label in ax.yaxis.get_ticklabels()[::2]:
        label.set_visible( False )

    if( current_ax[i] == 'A' or current_ax[i] == 'B' or current_ax[i] == 'C' ):    
        ax.xaxis.set_ticklabels([])
    if( current_ax[i] == 'B' or current_ax[i] == 'C' or current_ax[i] == 'E' or current_ax[i] == 'F' ):
        ax.yaxis.set_ticklabels([])
    
    # Add fiducial lines
    skew.plot_dry_adiabats( alpha = 0.35 )
    skew.plot_moist_adiabats( alpha = 0.35 )
    skew.plot_mixing_lines( alpha = 0.35 )
    skew.ax.axvline( 0, color = 'c', linestyle = '--', linewidth = 2.5, alpha = 0.35 )
    
    # Set the x-axis label
    if( current_ax[i] == 'D' or current_ax[i] == 'E' or current_ax[i] == 'F' ):
        ax.set_xlabel( 
                      xlabel = 'Temperature (\u00B0C)',
                      fontsize = 14, fontweight = 'bold'
                     )
    else:
        ax.set_xlabel( '' )
    
    # Set the y-axis label
    if( current_ax[i] == 'A' or current_ax[i] == 'D' ):
        ax.set_ylabel( 
                      ylabel = 'Pressure (hPa)',
                      fontsize = 14, fontweight = 'bold'
                     )
    else:
        ax.set_ylabel( '' )
    
    # Place a hodograph into the current axis
    hax = inset_axes( skew.ax, '35%', '35%', loc= "upper left" )
    
    # Label the hodograph
    hax.set_xlabel( 
                   'Hodograph (ms$^{-1}$)',
                    loc = 'center',
                    fontsize = 12, fontweight = 'bold'
                   )

    # Logic for Axis Title
    if( i < 3 ):
        ax.set_title( 
                     current_ax[i].lower() + ": " + sims[0] + " " + sound_type[i] + ' Composite', loc = 'left',
                     fontsize = 18, fontweight = 'bold'
                    )
        ax.set_title( 
                     bss_time[i], loc = 'right',
                     fontsize = 16, fontweight = 'bold'
                    )
    else:
        ax.set_title( 
                     current_ax[i].lower() + ": " + sims[1] + " " + sound_type[i] + ' Composite', loc = 'left',
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
    # ----------------------------------------------
    if( current_ax[i] == 'A' or current_ax[i] == 'D' ):
        
        mod_angle = units.Quantity( 25, 'degrees' )   
        
        # Get magnitude and direction from original components
        wspd = mpcalc.wind_speed( u, v )
        wdir = mpcalc.wind_direction( u , v )
        
        # Modify the direction of the hodograph evenly throughout all layers
        for j in range( 0, len( wdir) ):
            wdir[j] = wdir[j] + mod_angle
        
        # Compute the new wind components
        u2,v2 = mpcalc.wind_components( wspd, wdir )
        
        # Plot the hodograph (m/s) with a color map based on height
        hodo2 = h.plot( 
                      u = u2,
                      v = v2,
                      color = 'tab:brown',
                      linestyle = '-',
                      linewidth = 3,
                      alpha = 0.5,
                      )
        
        # Compute Bunkers Storm Motions and Mean Wind
        rm2, lm2, mean_uv2 = mpcalc.bunkers_storm_motion(
                                                      pressure = prs,
                                                      u = u2,
                                                      v = v2,
                                                      height = zh
                                                     )

        # Plot the mean wind as a vector
        h.wind_vectors( u = mean_uv2[0], v = mean_uv2[1], color = 'tab:brown', scale = 1, alpha = 0.5, width = 1, headwidth = 3 )
        
        # Plot the Bunker's Right & Left Motions as labeled pointd
        txt7 = plt.text( rm2[0].m, rm2[1].m, 'RM', size = 8, color = 'tab:brown', alpha  = 0.5 )
        # txt6 = plt.text( lm[0].m, lm[1].m, 'LM', size = 8 )
        
        # Place circles around RM & LM labels
        plt.scatter( rm2[0].m+2.5, rm2[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'tab:brown', alpha= 0.5 )
        # plt.scatter( lm[0].m+2.5, lm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
        
    # Set hodograph plot limits
    hax.set_xlim( -5, 45 )
    hax.set_ylim( -15, 25 )    
    
    # End Modify Hodograph
    #---------------------------------------------

    # Compute MLCAPE/CIN
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin( prs, temp, td )
    
    # Compute Bunkers Storm Motions and Mean Wind
    rm, lm, mean_uv = mpcalc.bunkers_storm_motion(
                                                  pressure = prs,
                                                  u = u,
                                                  v = v,
                                                  height = zh
                                                 )

    # Plot the mean wind as a vector
    h.wind_vectors( u = mean_uv[0], v = mean_uv[1], color = 'k', scale = 1, alpha = 0.9, width = 1, headwidth = 3 )
    
    # Plot the Bunker's Right & Left Motions as labeled pointd
    txt5 = plt.text( rm[0].m, rm[1].m, 'RM', size = 8, alpha = 0.9 )
    # txt6 = plt.text( lm[0].m, lm[1].m, 'LM', size = 8 )
    
    # Place circles around RM & LM labels
    plt.scatter( rm[0].m+2.5, rm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k', alpha = 0.9 )
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

    if( current_ax[i] == 'A' or current_ax[i] == 'D' ):
        p2 = np.arange( -57.5, 20, 6)
    else:
        p2 = np.arange( -75, 20, 6)
        
    txt1 = plt.text( 0.0, p2[3], 'MLCAPE: ' + str( round( float(mlcape.m), 2 ) ) + ' J kg$^{-1}$', fontsize = 14 )
    txt2 = plt.text( 0.0, p2[2], 'MLCIN: '  + str( round( float(mlcin.m),  2 ) ) + ' J kg$^{-1}$', fontsize = 14 )   
    txt3 = plt.text( 0.0, p2[1], '0-3 km SRH: ' + str( round( float(srh_3km[0].m), 2 ) ) + ' m$^2$ s$^{-2}$', fontsize = 14 )
    txt4 = plt.text( 0.0, p2[0], '0-1 km SRH: ' + str( round( float(srh_1km[0].m), 2 ) ) + ' m$^2$ s$^{-2}$', fontsize = 14 )
    
    # Create a mask for wind barb plotting
    mask = prs >= 100.0 * units.hPa
    
    # Determine which units of wind should be plotted
    skew.plot_barbs( prs[ mask ][::2], u[ mask ][::2], v[ mask ][::2] )     

    # Add legend
    if( current_ax[i] == 'A' or current_ax[i] == 'D' ):
        skew.ax.legend( loc = 3, facecolor = 'grey',  framealpha = 0.2, fontsize = 14 )
    
#------------------------------------------------------------------------------
# End: Loop to plot on each axis
#------------------------------------------------------------------------------

# Save the current skew-t
plt.savefig(
            fname = "/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR_Modeling_Project/Figures/Publication_Figs/RAW_Figs/CSTAR_Modeling_IC_PanelPlot_RAW_V2.jpeg",
            dpi = 300,
            bbox_inches = "tight"
            )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# End Main Script
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Confirm that script successfully ran
print( "\n{} successfully completed!".format( program ) )

# Report the time required to run the function
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )

