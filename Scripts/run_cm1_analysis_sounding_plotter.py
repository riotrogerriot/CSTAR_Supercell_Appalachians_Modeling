#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:15:09 2022

@author: roger
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 23:04:23 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Function: cm1_sounding_plotter
#
#   Update Records:
#       (3/22/22) Script Created
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.plots as plots
import metpy.calc as mpcalc
import pint
import pint_xarray
import cf_xarray


#-----------------------------------------------------------------------------
# Begin: Define cm1_sounding_plotter
#-----------------------------------------------------------------------------
def cm1_sounding_plotter( 
                         ds, fig, ax, sx, sy, shortname, current_model_time,
                         knots = False, terrain = False, peak_rel = False, 
                         peak_pos = 350.0
                        ):
    
    
    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the u-component of the wind (m/s)
    u = ds.metpy.parse_cf( 'uinterp' ).isel( time = 0, yh = sy, xh = sx )

    # Get the v-component of the wind (m/s)
    v = ds.metpy.parse_cf( 'vinterp' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get sounding location coordinates (m)
    xh = u.coords[ 'xh' ].values
    yh = u.coords[ 'yh' ].values
    
    # Convert location from m to km
    xh = ( xh / 1000.0 ) * units( 'km' )
    yh = ( yh / 1000.0 ) * units( 'km' )
    
    # Round the location values
    xh = round( xh.m, 0 )
    yh = round( yh.m, 0 )
    
    # Get potential temperature (K)
    th = ds.metpy.parse_cf( 'th' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get pressure (Pa)
    pressure = ds.metpy.parse_cf( 'prs' ).isel( time = 0, yh = sy, xh = sx )
    
    # Get water-vapor mixing ratio (kg/kg)
    qv = ds.metpy.parse_cf( 'qv' ).isel( time = 0, yh = sy, xh = sx ) 
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get heights (km)
    zh = ds.coords[ 'zh' ].values * units( 'km' )
    
    
    # Logic for computing height (AGL) for hodograph 
    #-----------------------------------------------
    if( terrain == True ):
        
        # Get model grid heights (km)
        zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0, yh = sy, xh = sx ) * units( 'm' )
        
        # Add native units
        zs = zs.metpy.quantify()
        
        # Convert to km
        zs = zs.pint.to( 'km' )

        # Compute Height (km AGL)
        zh_agl = zh - ( zs.values * units( 'm' ) )
    
    # If terrain is not included    
    else:
        
        # Height (AGL) is the same as the model grid height
        zh_agl = zh
    #-----------------------------------------------
    
    # Logic to determine if x-axis should be peak-relative or not
    #------------------------------------------------------------
    if( peak_rel == True ):
        
        # Convert x-coordinates from grid-relative to peak-relative coordinates
        xh = xh - peak_pos 
        
    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
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
    
    # Convert winds to knots if desired
    #----------------------------------
    if( knots == True ):
        u = u.pint.to( 'knots' )
        v = v.pint.to( 'knots' )
    #----------------------------------
    
    #-----------------------------------------------------------------------------
    # End: Unit conversions required for MetPy functions
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Air Temperature & Dew Point Calculations
    #-----------------------------------------------------------------------------
    
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
    
    #-----------------------------------------------------------------------------
    # End: Air Temperature & Dew Point Calculations
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: SkewT Plotting
    #-----------------------------------------------------------------------------
    
    # Create SkewT object
    skew = plots.SkewT( fig, rotation = 45, subplot = ax, aspect = 100 )
    
    # Compute the surface-based parcel path (degC)
    prof = mpcalc.parcel_profile(
                                 pressure,
                                 temp[0],
                                 td[0]
                                ).pint.to( 'degC' )
    
    # Calculate LCL height
    lcl_pressure, lcl_temperature = mpcalc.lcl( pressure[0], temp[0], td[0] )
    
    # Calculate LFC height
    lfc_pressure, lfc_temperature = mpcalc.lfc( pressure, temp, td )
    
    # Calculate Equilibrium Level
    el_pressure, el_temperature = mpcalc.el( pressure, temp, td )
    
    
    # Plot temp, td, and parcel path
    skew.plot( pressure, temp, 'red', linewidth = 2.5, label = 'Temp' )
    skew.plot( pressure, td, 'green', linewidth = 2.5, label = 'Dewpoint' )
    skew.plot( pressure, prof, '--', linewidth = 2.5, color = 'k', label = 'Parcel Path' )
    
    # Plot LCL, LFC, and EL levels
    skew.plot( lcl_pressure, lcl_temperature.to('degC'), marker = '_', markersize = 40, color = 'b' )
    skew.plot( lfc_pressure, lfc_temperature.to('degC'), marker = '_', markersize = 40, color = 'r' )
    skew.plot( el_pressure, el_temperature.to('degC'), marker = '_', markersize = 40, color = 'k' )
    
    # Label LCL, LFC, and EL levels
    plt.text( lcl_temperature.to( 'degC').m + 5, lcl_pressure.m + 25 , 'LCL', color = 'b', size = 18 )
    plt.text( lfc_temperature.to( 'degC').m + 5, lfc_pressure.m + 20 , 'LFC', color = 'r', size = 18 )
    plt.text( el_temperature.to( 'degC').m + 5, el_pressure.m + 5, 'EL', color = 'k', size = 18 )
        
    
    # Add fiducial lines
    skew.plot_dry_adiabats( alpha = 0.35 )
    skew.plot_moist_adiabats( alpha = 0.35 )
    skew.plot_mixing_lines( alpha = 0.35 )
    skew.ax.axvline( 0, color = 'c', linestyle = '--', linewidth = 2.5, alpha = 0.35 )
    
    # Shade CAPE/CIN (Needed values to work!)
    skew.shade_cin( pressure.values * units( 'hPa' ), temp.values * units( 'degC' ), prof.values * units( 'degC' ), label = 'SBCIN' )
    skew.shade_cape( pressure.values * units( 'hPa' ), temp.values * units( 'degC' ), prof.values * units( 'degC' ), label = 'SBCAPE' )
        
    # Reformat axes labels
    skew.ax.set_xlabel( 'Temperature (\u00b0C)', fontsize = 18, fontweight = 'bold' )
    skew.ax.set_ylabel( 'Pressure (hPa)', fontsize = 18, fontweight = 'bold' )
    
    # Create Legend
    skew.ax.legend( loc = 3, facecolor = 'grey',  framealpha = 0.75, fontsize = 14 )
        
    #-----------------------------------------------------------------------------
    # End: SkewT Plotting
    #-----------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------
    # Begin: Hodograph Plotting
    #-----------------------------------------------------------------------------
    
    # Create a mask to limit wind barbs to plotting limits
    mask_uv = pressure >= 100.0 * units( 'hPa' )
    
    # Plot Wind Barbs
    skew.plot_barbs( pressure[mask_uv].values , u[mask_uv].values, v[mask_uv].values )
    
    
    # Height (km AGL) intervals used for hodograph plotting
    intervals = np.array( [ 0.0, 3.0, 6.0, 9.0, 12.0, 15.0 ] ) * units( 'km' )
    
    # Height (km AGL) Mask
    zh_mask = zh_agl <= 15.0 * units( 'km' )
    
    # Colormap scheme to match interval bins for hodograph
    cmap = [ 'tab:red', 'tab:green', 'tab:olive', 'tab:cyan', 'tab:purple' ]
    
    # Insert the hodograph axis
    hodo_ax = inset_axes( skew.ax, '40%', '30%' )
    
    # Plot the hodograph on the new inset axes
    hodo = plots.Hodograph( hodo_ax, component_range = 100.0 )
    
    # Plot the hodograph
    hodograph = hodo.plot_colormapped( u[zh_mask], v[zh_mask], c = zh_agl[ zh_mask ], intervals = intervals, colors = cmap )
    
    
    # Compute Bunkers Storm Motions and Mean Wind
    rm, lm, mean_uv = mpcalc.bunkers_storm_motion( pressure, u, v, zh_agl )
    
    # Plot the mean wind as a vector
    hodo.wind_vectors( u = mean_uv[0], v = mean_uv[1], color = 'dimgrey', scale = 1, alpha = 0.75, width = 1, headwidth = 3 )
    
    # Plot the Bunker's Right & Left Motions as labeled pointd
    plt.text( rm[0].m, rm[1].m, 'RM', size = 8 )
    plt.text( lm[0].m, lm[1].m, 'LM', size = 8 )
    
    
    # Logic for hodograph Plot Parameters
    #------------------------------------
    
    # Parameters if using wind units (kts)
    if( knots == True ):
        
        # Hodograph Ring Increments
        hodo.add_grid( increment = 20.0 )
        
        # Hodograph axes limits
        hodo_ax.set_xlim( -20, 80 )
        hodo_ax.set_ylim( -50, 70 )
        
        # Place circles around RM & LM labels
        plt.scatter( rm[0].m+5, rm[1].m+2.5, s= 250, facecolor = 'none', edgecolor = 'k' )
        plt.scatter( lm[0].m+5, lm[1].m+2.5, s= 250, facecolor = 'none', edgecolor = 'k' )
       
        # Add x-axis label 
        hodo_ax.set_xlabel( 
                           'Hodograph (kts)',
                            loc = 'center',
                            fontsize = 10,
                            fontweight = 'bold'
                           )
    
    # Parameters if using wind units (m/s)
    else:
        
        # Hodograph Ring Increments
        hodo.add_grid( increment = 10.0 )
        
        # Hodograph axes limits
        hodo_ax.set_xlim( -10, 40 )
        hodo_ax.set_ylim( -25, 35 )
        
        # Place circles around RM & LM labels
        plt.scatter( rm[0].m+2.5, rm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
        plt.scatter( lm[0].m+2.5, lm[1].m+1.25, s= 250, facecolor = 'none', edgecolor = 'k' )
        
        # Add x-axis label
        hodo_ax.set_xlabel( 
                           'Hodograph (ms$^{-1}$)',
                           loc = 'center',
                           fontsize = 10,
                           fontweight = 'bold'
                          )
    #------------------------------------    
    
    
    # Add plot title
    ax.set_title( 
                 "CM1-{}: t ={} min.\nb: Near-Storm Inflow Sounding".format (shortname, current_model_time),
                 loc = "left", fontsize = 14, fontweight = 'bold' 
                ) 
    ax.set_title( 
                 "x = {} km, y = {} km".format (xh, yh),
                 loc = "right", fontsize = 14, fontweight = 'bold' 
                ) 
    
    
    #-----------------------------------------------------------------------------
    # End: Hodograph Plotting
    #-----------------------------------------------------------------------------
    
    # Return the axes object 
    return( ax )

#-----------------------------------------------------------------------------
# End: Define cm1_sounding_plotter
#-----------------------------------------------------------------------------

