#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 18:33:04 2022

@author: rriggin
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from metpy.units import units
import metpy.calc as mpcalc

import pint
import pint_xarray
import cf_xarray

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.plots as plots

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

# Working directory
# main_dir = '/scratch/rriggin/cm1r20.3/original_sims_req_files/Crosser/Idealized_Terrain_NoBSS/dry_run'
main_dir = '/Users/roger/Desktop/'

# Composite Environment Label
env = 'Crossing'

# Number of gridpoints in x-direction
nx = 2400

# Number of gridpoints in y-direction
ny = 1600

# Number of gridpoints in z-direction
nz = 48

peak_elv = 350000.0

# Moving window to compute mesocyclone area over (km^2)
meso_window = 35.0

# Common Height Indices (Must change if z-stretch setting are modified)
m500 = 4
km1 = 6
km3 = 13
km5 = 18
km8 = 24
km10 = 28

current_ax = [ 'B', 'C', 'D' ]
sounding_loc = [ -200, -100, -50 ]

# Get smooth terrain profile
x, y, hx, HX, HY = witch_of_agnesi( x0 = peak_elv )

hx = hx/1000.0 * units( 'km' )

# Convert to peak-relative
x = x - peak_elv

# Change the working directory to the provided path
os.chdir( main_dir )


# Loop through specified model output files
#------------------------------------------------------------------------------
for i in np.arange( 25, 26 ):
    
    # Get current integration time
    current_model_time = (i-1) * 5.0 
    
    # Logic to open standard cm1out netCDF
    #-------------------------------------
    
    # # Define the current filename as a string
    # if ( i < 10 ):
    #     filename = 'cm1out_00000' + str(i) + '.nc'
    #     print( "\tOpening " + filename + "..." )
        
    # elif( i < 100 and i >= 10 ):
    #     filename = 'cm1out_0000' + str(i) + '.nc'
    #     print( "\nOpening " + filename + "...")  
    
    # else:
    #     filename = 'cm1out_000' + str(i) + '.nc'
    #     print( "\nOpening " + filename + "...")

    # # Open the current netCDF file with xarray 
    # DS = xr.open_dataset( filename, engine = "netcdf4", decode_cf = True )    
    
    
    # Define the terrain-interpolated filename as a string
    if ( i < 10 ):
        filename2 = 'cm1out_00000' + str(i) + '_i.nc'
    elif( i < 100 and i >= 10 ):
        filename2 = 'cm1out_0000' + str(i) + '_i.nc'
    else:
        filename = 'cm1out_000' + str(i) + '_i.nc'
    
    # End: Logic to open standard cm1out netCDF
    #-------------------------------------

    
    # Open the current terrain-interpolated netCDF file with xarray 
    DS_interp = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
    
    
    # Coordinate Variables
    #-----------------------------------------------------------------------------
    
    # Get the x-axis coordinates (km)
    xh = DS_interp.coords[ 'xh' ].isel( xh = slice( 0, nx ) ).values * units( 'km' )
    
    # Convert x-coordinates from grid-relative to peak-relative coordinates
    xh = xh - peak_elv / 1000.0 * units( 'km' )
    
    # Get the y-axis coordinates (km)
    yh = DS_interp.coords[ 'yh' ].isel( yh = slice( 0, ny ) ).values * units( 'km' )
    
    # Get the y-axis coordinates (km)      
    zh = DS_interp.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' ) 
    
    
    # Collect raw and derived variables 
    #-----------------------------------------------------------------------------

    # Get the u-component of the wind (m/s)    
    u = DS_interp.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    u = u.metpy.quantify()
    u_mask = np.ma.array( u, mask = u == -999999.875 )
    u_arr = [ u_mask[:,800,600], u_mask[:,800,1000], u_mask[:,800,1200] ]
    
    # Get the v-component of the wind (m/s)    
    v = DS_interp.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    v = v.metpy.quantify()
    v_mask = np.ma.array( v, mask = u == -999999.875 )
    v_arr = [ v_mask[:,800,600], v_mask[:,800,1000], v_mask[:,800,1200] ]
    
    # Get pressure (Pa), mask out terrain, convert to hPa, and select sounding points
    pressure = DS_interp.metpy.parse_cf( 'prs' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    prs_mask = pressure.where( pressure > -999999.875 ) 
    prs_mask = units.Quantity( prs_mask.values/100.0, 'hPa' )
    prs_arr = [ prs_mask[:,800,600], prs_mask[:,800,1000], prs_mask[:,800,1200] ]
    
    # Get potential temperature perturbation (K) and mask out terrain
    th = DS_interp.metpy.parse_cf( 'th' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    th = th.metpy.quantify()
    th_mask = np.ma.array( th, mask = th == -999999.875 ) 
    th_arr = [ th_mask[:,800,600], th_mask[:,800,1000], th_mask[:,800,1200] ] * units( 'K' )
    
    # Get water-vapor mixing ratio (kg/kg) and convert to g/kg
    qv = DS_interp.metpy.parse_cf( 'qv' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    qv = qv.metpy.quantify()
    qv = qv.pint.to( 'g/kg' )
    
    # Convert potential temperature to air temperature (K), convert to degC, and select sounding points
    temp = mpcalc.temperature_from_potential_temperature( pressure, th ) 
    temp = temp.metpy.quantify()
    temp = temp.pint.to( 'degC' )
    temp_arr = [ temp[:,800,600], temp[:,800,1000], temp[:,800,1200] ]
    
    # Compute the water vapor pressure using pressure and mixing-ratio
    e  = mpcalc.vapor_pressure( pressure, qv )
    
    # Compute dew point from vapor pressure (degC) and select sounding points
    td = mpcalc.dewpoint( e ) 
    td_arr = [ td[:,800,600], td[:,800,1000], td[:,800,1200] ]
    
    # Remove NaNs from second and third sounding points for proper plotting
    prs_arrs = [ prs_arr[0], prs_arr[1][1:], prs_arr[2][3:] ]
    temp_arrs = [ temp_arr[0], temp_arr[1][1:], temp_arr[2][3:] ]
    td_arrs = [ td_arr[0], td_arr[1][1:], td_arr[2][3:] ] 
    u_arrs = [ u_arr[0], u_arr[1][1:], u_arr[2][3:] ]
    v_arrs = [ v_arr[0], v_arr[1][1:], v_arr[2][3:] ]
    th_arrs = [ th_arr[0], th_arr[1][1:], th_arr[2][3:] ]
    zh_arrs = [ zh[:], zh[1:] - hx[1000], zh[3:] - hx[1200] ]
    

    # Create plot
    #-----------------------------------------------------------------------------

    # Set tick size (Must be before calling plot object)
    plt.rcdefaults()
    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    
    # Create figure templates to add data to during loop
    figure_mosaic = """
                    AAA
                    AAA
                    BCD
                    BCD
                    """
                    
    # Create Figure and Axes objects                
    fig, axes =  plt.subplot_mosaic(
                                    mosaic = figure_mosaic,
                                    figsize = ( 14, 12 ),
                                    tight_layout = True
                                   )

    # Create callable axes objects and set main labels
    ax_A = fig.add_subplot( axes['A'] )
    ax_A.set_title( 'a: {} Environment'.format( env ), fontsize = 20, fontweight = 'bold', loc = 'left' )
    ax_A.set_title( 'y = 200 km, t = {} min.'.format( current_model_time ), fontsize = 20, fontweight = 'bold', loc = 'right' )
    ax_A.set_ylabel( 'Elevation (km)', fontsize = 18, fontweight = 'bold' )
    ax_A.set_xlabel( 'Distance from Peak (km)', fontsize = 18, fontweight = 'bold' )
    
    # Set axis_A limits
    ax_A.set_xlim( -350, 250 )
    ax_A.set_ylim( 0, 3.0 )
    
    # Draw lines where soundings are pulled
    ax_A.axvline( x = -200, linestyle = '--', color = 'k' )
    ax_A.axvline( x = -100, linestyle = '--', color = 'k' )
    ax_A.axvline( x = -50, linestyle = '--', color = 'k' )
    
    # Plot the terrain profile
    terrain_profile = ax_A.plot( x/1000.0, hx, linewidth = 2.5, fillstyle = 'full', color = 'k' )
    
    # Fill in the area under the terrain curve
    ax_A.fill_between( x/1000.0, hx, color = 'k' )

    # Plot contours of potential temperature
    c1 = ax_A.contour(
                      xh, zh, th_mask[:,800,:],
                      levels = np.arange( 300, 313, 1 ),
                      linestyles = 'solid', linewidths = 2.0, 
                      colors = 'black', alpha = 1.0
                     )
   
    # Plot filled contours of potential temperature
    cf1 = ax_A.contourf( 
                       xh, zh, th_mask[ :, 800, : ], 
                       levels = np.arange( 300, 313, 1 ), 
                       linestyles = "solid",
                       cmap = "binary",
                       alpha = 0.75
                      )
    
    # Include color bar legend
    cbar = plt.colorbar( cf1, location = 'bottom', fraction = 0.1, shrink = 0.75, pad = 0.15, drawedges = True ) 
    
    # Colorbar label
    cbar.set_label( label = '\u03F4 (K)', fontsize = 16, fontweight = 'bold' )
    
    # Add grid
    ax_A.grid()

    # Loop through each sounding point
    #-----------------------------------------------------------------------------
    for j in range( 0, len( current_ax ) ):
        
        # Axes set-up
        #-----------------------------------------------------------------------------
        
        # Get the original specs for subplot D
        ss = axes[ current_ax[j] ].get_subplotspec()
        
        # Remove the original instance of subplot D
        axes[ current_ax[j] ].remove()
        
        # Re-construct the subplot using the original specs with the skewed x-axis projection 
        axes[ current_ax[j] ] = fig.add_subplot( ss, projection = 'skewx' )
        
        # Call the current axes object
        ax = fig.add_subplot( axes[ current_ax[j] ] )
        
        # Create a SkewT plotting object
        skew = plots.SkewT( fig, rotation = 45, subplot = ax, aspect = 100 )
        
        # Custom axes limits
        skew.ax.set_ylim( 1000, 100 )
        skew.ax.set_xlim( -70, 40 ) 
        
        # Add fiducial lines
        skew.plot_dry_adiabats( alpha = 0.35 )
        skew.plot_moist_adiabats( alpha = 0.35 )
        skew.plot_mixing_lines( alpha = 0.35 )
        skew.ax.axvline( 0, color = 'c', linestyle = '--', linewidth = 2.5, alpha = 0.35 )
        
        # Set the x-axis label
        ax.set_xlabel( xlabel = 'Temperature (\u00B0C)', fontsize = 14, fontweight = 'bold' )
        
        # Set the y-axis labels
        if( j == 0 ):
            ax.set_ylabel( ylabel = 'Pressure (hPa)', fontsize = 14, fontweight = 'bold' )
        else:
            ax.set_ylabel( ylabel = '' )
         
        # Set the sounding plot titles
        ax.set_title( str( current_ax[j].lower() ) + ":  " , loc = 'left', fontsize = 18, fontweight = 'bold' )
        ax.set_title( 'x = ' + str( sounding_loc[j] ) + ' km', loc = 'right', fontsize = 16, fontweight = 'bold' )
            
        # Label every other tickmark    
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in ax.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        
        
        # Sounding Calculations & Plotting
        #-----------------------------------------------------------------------------
        
        # Compute the surface-based parcel path (degC)
        prof = mpcalc.parcel_profile( prs_arrs[j], temp_arrs[j][0], td_arrs[j][0] )
        
        # Calculate LCL height
        lcl_pressure, lcl_temperature = mpcalc.lcl( prs_arrs[j][0], temp_arrs[j][0], td_arrs[j][0] )
        
        # Calculate LFC height
        lfc_pressure, lfc_temperature = mpcalc.lfc( prs_arrs[j], temp_arrs[j], td_arrs[j] )
        
        # Calculate Equilibrium Level
        el_pressure, el_temperature = mpcalc.el( prs_arrs[j], temp_arrs[j], td_arrs[j] )
        
        # Plot temp, td, and parcel path
        skew.plot( prs_arrs[j], temp_arrs[j], 'red', linewidth = 2.5, label = 'Temp' )
        skew.plot( prs_arrs[j], td_arrs[j], 'green', linewidth = 2.5, label = 'Dewpoint' )
        skew.plot( prs_arrs[j], prof, '--', linewidth = 2.5, color = 'k', label = 'Parcel Path' )
        
        # Plot LCL, LFC, and EL levels
        skew.plot( lcl_pressure, lcl_temperature.to('degC'), marker = '_', markersize = 30, color = 'b' )
        skew.plot( lfc_pressure, lfc_temperature.to('degC'), marker = '_', markersize = 30, color = 'r' )
        skew.plot( el_pressure, el_temperature.to('degC'), marker = '_', markersize = 30, color = 'k' )
        
        # Shade CAPE/CIN 
        skew.shade_cin( prs_arrs[j].m * units( 'hPa' ), temp_arrs[j].values * units( 'degC' ), prof.to( 'degC' ), label = 'SBCIN' )
        skew.shade_cape( prs_arrs[j].m * units( 'hPa' ), temp_arrs[j].values * units( 'degC' ), prof.to( 'degC' ), label = 'SBCAPE' )
        
        # Create a mask for wind barbs
        mask_uv = prs_arrs[j] >= 100.0 * units( 'hPa' )
        
        # Plot every other wind barb through 100 hPa
        skew.plot_barbs( prs_arrs[j][mask_uv][::2], u_arrs[j][mask_uv][::2], v_arrs[j][mask_uv][::2] )
       
       
       
        # Height (km AGL) intervals used for hodograph plotting
        intervals = np.array( [ 0.0, 3.0, 6.0, 9.0, 12.0, 15.0 ] ) * units( 'km' )
        
        # Height (km AGL) Mask
        zh_mask = zh_arrs[j] <= 15.0 * units( 'km' )
        
        # Colormap scheme to match interval bins for hodograph
        cmap = [ 'tab:red', 'tab:green', 'tab:olive', 'tab:cyan', 'tab:purple' ]
        
        # Insert the hodograph axis
        hodo_ax = inset_axes( skew.ax, '30%', '30%' )
        
        # Plot the hodograph on the new inset axes
        hodo = plots.Hodograph( hodo_ax, component_range = 100.0 )
        
        # Plot the hodograph
        hodograph = hodo.plot_colormapped( u_arrs[j][zh_mask], v_arrs[j][zh_mask], c = zh_arrs[j][ zh_mask ], intervals = intervals, colors = cmap )       
        
        # Hodograph Ring Increments
        hodo.add_grid( increment = 10.0 )
        
        # Hodograph axes limits
        hodo_ax.set_xlim( -10, 40 )
        hodo_ax.set_ylim( -25, 35 )
        
       
        # Brunt Vasilia (N) and Frode Number (Fr) Calculations
        #-----------------------------------------------------------------------------        
        
        mlcape, mlcin = mpcalc.mixed_layer_cape_cin( prs_arrs[j], temp_arrs[j], td_arrs[j] )
        
        lfc_lcl_delta = lfc_pressure - lcl_pressure
        
        # Compute N (Data Array) through terrain depth (zh[6] ~ 1 km)
        N = mpcalc.brunt_vaisala_frequency( height = zh[0:km1].to('m'), potential_temperature = th_arrs[j][0:km1] )
        
        # Average the results of N for a single value representative of the 0-3 km Layer
        N_avg = np.average( N )
        
        ax.text( x = 0.025, y = 0.7, s = 'CAPE: ' + str( round(mlcape.m,0) ), fontsize = 14, transform = ax.transAxes )
        ax.text( x = 0.025, y = 0.6, s = 'CIN: ' + str( round(mlcin.m,0) ), fontsize = 14, transform = ax.transAxes )
        ax.text( x = 0.025, y = 0.5, s = 'LFC-LCL: ' + str( round(lfc_lcl_delta.m,0) ), fontsize = 14, transform = ax.transAxes )
        ax.text( x = 0.025, y = 0.4, s = 'N1KM: ' + str( round(N_avg.m,3) ), fontsize = 14, transform = ax.transAxes )
        
        # Average the zonal wind through the 0-3 km Layer
        u_avg = np.average( u_arrs[j][0:km3] ) * units( 'm/s' )
        u_avg_1km =np.average( u_arrs[j][0:km1] ) * units( 'm/s' )
        
        Fr_avg_3km = ( u_avg )/ ( N_avg * np.amax( hx ).to( 'm' ) )
        ax.text( x = 0.025, y = 0.3, s ='F3KM: ' + str(round(Fr_avg_3km.m, 2) ), fontsize = 14, transform = ax.transAxes )
        
        Fr_avg_1km = ( u_avg_1km )/ ( N_avg * np.amax( hx ).to( 'm' ) )
        ax.text( x = 0.025, y = 0.2, s = 'F1KM: ' + str( round( Fr_avg_1km.m, 2 ) ), fontsize = 14, transform = ax.transAxes )
        
        Fr_surface = ( u_arrs[j][0] * units( 'm/s' ) )/ ( N_avg * np.amax( hx ).to( 'm' ) )
        ax.text( x = 0.025, y = 0.1, s = 'F0KM: ' + str( round(Fr_surface.m, 2 ) ), fontsize = 14, transform = ax.transAxes )
        
    # End loop through each sounding point
    #-----------------------------------------------------------------------------
            
    # Show current plot    
    plt.show()
    
    # Save figure
    fig.savefig(
                fname = '{}_Dry_Run_Analysis_{}_min{}_.jpeg'.format( env, int(current_model_time), i ),
                dpi = 200,
                bbox_inches = "tight"
               )
    
# End: Loop through specified model output files
#------------------------------------------------------------------------------
