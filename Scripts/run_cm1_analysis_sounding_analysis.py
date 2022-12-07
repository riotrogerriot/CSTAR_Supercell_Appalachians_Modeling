#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:57:44 2022

@author: rriggin
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
def cm1_sounding_analysis( 
                          ds, sx, sy, shortname, current_model_time, 
                          rm, lm, mean_uv, km3, terrain = False
                          # , peak_rel = False, peak_pos = 350.0
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
    
    # # Get sounding location coordinates (m)
    # xh = u.coords[ 'xh' ].values
    # yh = u.coords[ 'yh' ].values
    
    # # Convert location from m to km
    # xh = ( xh / 1000.0 ) * units( 'km' )
    # yh = ( yh / 1000.0 ) * units( 'km' )
    
    # # Round the location values
    # xh = round( xh.m, 0 )
    # yh = round( yh.m, 0 )
    
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

        # # Compute Height (km AGL)
        # zh_agl = zh - ( zs.values * units( 'km' ) )
    
    # If terrain is not included    
    else:
        
        # Height (AGL) is the same as the model grid height
        zs = zh
    #-----------------------------------------------
    
    # # Logic to determine if x-axis should be peak-relative or not
    # #------------------------------------------------------------
    # if( peak_rel == True ):
        
    #     # Convert x-coordinates from grid-relative to peak-relative coordinates
    #     xh = xh - peak_pos 
        
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
    # Begin: Compute Sounding Parameters
    #-----------------------------------------------------------------------------
    
    # Calculate LCL height
    lcl_pressure, lcl_temperature = mpcalc.lcl( pressure[0], temp[0], td[0] )
    
    # Calculate LFC height
    lfc_pressure, lfc_temperature = mpcalc.lfc( pressure, temp, td )
    
    # Calculate Equilibrium Level
    el_pressure, el_temperature = mpcalc.el( pressure, temp, td )
    
    # Compute SBCAPE/CIN
    sbcape, sbcin = mpcalc.surface_based_cape_cin( pressure, temp, td )
   
    # Compute MLCAPE/CIN
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin( pressure, temp, td )
   
    # Compute MUCAPE/CIN 
    mucape, mucin = mpcalc.most_unstable_cape_cin( pressure, temp, td )
    
    # 3CAPE
    cape3km, cin3km = mpcalc.surface_based_cape_cin( pressure[0:km3], temp[0:km3], td[0:km3] )
    
    
    
    # Compute 0-1 km shear components
    shear_u_1km, shear_v_1km = mpcalc.bulk_shear(
                                                 pressure = pressure,
                                                 u = u,
                                                 v = v,
                                                 height = zh,
                                                 depth = units.Quantity( 1000, 'm' )
                                                )
            
    # Get 0-1 km Shear magnitude & direction
    shear1km = np.sqrt( shear_u_1km **2 + shear_v_1km ** 2 )
    shear_1km_dir = mpcalc.wind_direction( u = shear_u_1km, v = shear_v_1km )
    
    
    # Compute 0-3 km shear components
    shear_u_3km, shear_v_3km = mpcalc.bulk_shear(
                                                 pressure = pressure,
                                                 u = u,
                                                 v = v,
                                                 height = zh,
                                                 depth = units.Quantity( 3000, 'm' )
                                                )
    
    # Get 0-3 km Shear magnitude & direction
    shear3km = np.sqrt( shear_u_3km **2 + shear_v_3km ** 2 )
    shear_3km_dir = mpcalc.wind_direction( u = shear_u_3km, v = shear_v_3km )
    
    
    # Compute 0-6 km shear components
    shear_u_6km, shear_v_6km = mpcalc.bulk_shear(
                                                 pressure = pressure,
                                                 u = u,
                                                 v = v,
                                                 height = zh,
                                                 depth = units.Quantity( 6000, 'm' )
                                                )
    
    # Get 0-6km shear magnitude & direction
    shear6km = np.sqrt( shear_u_6km **2 + shear_v_6km ** 2 )
    shear_6km_dir = mpcalc.wind_direction( u = shear_u_6km, v = shear_v_6km )
    
    
    
    # Compute bunkers magnitude
    rm_mag = mpcalc.wind_speed( u = rm[0], v = rm[1] )
    lm_mag = mpcalc.wind_speed( u = lm[0], v = lm[1] )
    
    # Compute bunkers directions
    rm_dir = mpcalc.wind_direction( u = rm[0], v = rm[1] )
    lm_dir = mpcalc.wind_direction( u = lm[0], v = lm[1] )
    
    # Get mean wind magnitude & direction
    mean_wind = mpcalc.wind_speed( mean_uv[0], mean_uv[1] )
    mean_uv_dir = mpcalc.wind_direction( u = mean_uv[0], v = mean_uv[1] )
    
    # Compute 0-0.5 km SRH
    srh_500m = mpcalc.storm_relative_helicity(
                                             u = u,
                                             v = v,
                                             height = zh,
                                             depth = units.Quantity( 500, 'm' ),
                                             storm_u = rm[0],
                                             storm_v = rm[1]
                                            )
    
    # Compute 0-1 km SRH
    srh_1km = mpcalc.storm_relative_helicity(
                                             u = u,
                                             v = v,
                                             height = zh,
                                             depth = units.Quantity( 1000, 'm' ),
                                             storm_u = rm[0],
                                             storm_v = rm[1]
                                            )
    
    # Compute 0-3 km SRH
    srh_3km = mpcalc.storm_relative_helicity(
                                             u = u,
                                             v = v,
                                             height = zh,
                                             depth = units.Quantity( 3000, 'm' ),
                                             storm_u = rm[0],
                                             storm_v = rm[1]
                                            )

    # Compute Critical Angle
    crit_angle = mpcalc.critical_angle(
                                       pressure = pressure,
                                       u = u,
                                       v = v,
                                       height = zh,
                                       u_storm = rm[0],
                                       v_storm = rm[1]
                                      )
        
    #-----------------------------------------------------------------------------
    # End: Compute Sounding Parameters
    #----------------------------------------------------------------------------
    
    # Return the axes object 
    return(
           zs, lcl_pressure, lfc_pressure, el_pressure, 
           sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3km,
           shear1km, shear_1km_dir, shear3km, shear_3km_dir, shear6km, shear_6km_dir,
           rm_mag, rm_dir, lm_mag, lm_dir, mean_wind, mean_uv_dir,
           srh_500m, srh_1km, srh_3km, crit_angle
          )

#-----------------------------------------------------------------------------
# End: Define cm1_sounding_plotter
#-----------------------------------------------------------------------------