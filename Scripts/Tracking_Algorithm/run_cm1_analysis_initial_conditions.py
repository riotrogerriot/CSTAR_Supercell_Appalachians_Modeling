#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 14:21:26 2022

@author: roger
"""

#-----------------------------------------------------------------------------
# Function: cm1_initial_conditions
#
#   Update Records:
#       (3/20/22) Script Created
#
#-----------------------------------------------------------------------------
#
# Summary:  
#
#   This function pulls out variables of interest from the initial conditions
#   netcdf dataset  from CM1. This function may require modifications over time
#   to include variables required for future applications.
# 
# Arguments:
#
#   ds: A given CM1 model output xarray dataset
#   sx: The X-index for the given sounding location
#   sy: the Y-index for the given sounding location
#
# Returns:
#   mean_uv: A pint array containing the mean wind components
#   rm: A pint array containing the Bunker's Right components
#   
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Load external libraries
#-----------------------------------------------------------------------------
from metpy.units import units
import metpy.calc as mpcalc


#-----------------------------------------------------------------------------
# Begin: Define cm1_initial_conditions Function
#-----------------------------------------------------------------------------
def cm1_initial_conditions( ds, sx, sy ):
    
    # Pull out the wind-components for the initial sounding location 
    u = ds.metpy.parse_cf( 'uinterp' ).isel( time = 0, yh = sy, xh = sx )
    v = ds.metpy.parse_cf( 'vinterp' ).isel( time = 0, yh = sy, xh = sx )
    
    # Pull out pressure along the vertical axis of the given sounding location
    pressure = ds.metpy.parse_cf( 'prs' ).isel( time = 0, yh = sy, xh = sx )
    
    # Pull out the height values
    height = ds.coords[ 'zh' ].values * units( 'km' )
    
    # Compute Bunkers Storm Motions and Mean Wind
    rm, lm, mean_uv = mpcalc.bunkers_storm_motion( pressure, u, v, height )
    
    # Return the computed storm-motion vectors    
    return( mean_uv, rm, lm )

#-----------------------------------------------------------------------------
# End: Define cm1_initial_conditions Function
#-----------------------------------------------------------------------------

