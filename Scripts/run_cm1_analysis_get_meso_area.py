#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 17:02:30 2022

@author: rriggin
"""

#-----------------------------------------------------------------------------
# Function: cm1_get_meso_area.py
#
#   Update Records:
#       (6/28/22) Script Created
#
# Notes:
#
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

import pint
import pint_xarray
import cf_xarray


#-----------------------------------------------------------------------------
# Begin: Define cm1_uh_tracker
#-----------------------------------------------------------------------------
def cm1_get_meso_area( 
                      idx, ds, nx, ny, nz, dx, dy, uh_x, uh_y, z, 
                      meso_window = 20.0, track_on = 24, 
                      w_thresh = 10.0, zeta_thresh = 0.01, w_zeta_thresh = 0.1, dbz_thresh = 30.0,
                      m500 = 4, km1 = 6, km3 = 13, km5 = 18, km8 = 24, km10 = 28,
                      terrain = False,
                     ):
    
    #-----------------------------------------------------------------------------
    # Begin: Get variables from xarray dataset
    #-----------------------------------------------------------------------------
    
    # Data Variables
    #-----------------------------------------------------------------------------
    
    # Get the w-component of the wind (m/s)   
    w = ds.metpy.parse_cf( 'winterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the vertical vorticity (1/s)   
    zvort = ds.metpy.parse_cf( 'zvort' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get the z-axis coordinates (km)      
    zh = ds.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' )
    
    # Get the radar reflectivity (dBZ)
    dbz = ds.metpy.parse_cf( 'dbz' ).isel( time = 0, zh = z, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Get potential temperature pertubation at surface (K)
    thpert = ds.metpy.parse_cf( 'thpert' ).isel( time = 0, zh = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
    
    # Drop units to speed up looping
    w = w.values
    zvort = zvort.values
    thpert = thpert.values
    
    # Compute rotating updraft metric (m/s^2)
    w_zeta = w * zvort

    #-----------------------------------------------------------------------------
    # Begin: Compute Mesocyclone Area Estimates
    #-----------------------------------------------------------------------------
    
    # Define the search radius about the uh max as a callable integer value for looping purposes (Units: m)
    meso_radius = int( ( (meso_window / 2.0 ) * 1000.0 ) / dx )
    
    # Initialize area values for each level 
    
    thpert_surface = 0.0
    
    w_area_500m = 0.0
    zeta_area_500m = 0.0
    w_zeta_area_500m = 0.0
    
    w_area_1km = 0.0
    zeta_area_1km = 0.0
    w_zeta_area_1km = 0.0
    
    w_area_3km = 0.0
    zeta_area_3km = 0.0
    w_zeta_area_3km = 0.0
    
    w_area_5km = 0.0
    zeta_area_5km = 0.0
    w_zeta_area_5km = 0.0
    
    w_area_8km = 0.0
    zeta_area_8km = 0.0
    w_zeta_area_8km = 0.0
    
    dbz_area = 0.0
    
    # Get the looping indices 
    i_start = uh_x - meso_radius
    i_end = uh_x + meso_radius
    j_start = uh_y - meso_radius
    j_end = uh_y + meso_radius
    
    if( i_end >= nx ):
        i_end = nx -1
        
    if( j_end  >= ny ):
        j_end = ny - 1
    
    # Put mandatory levels in an array for looping purposes
    surface = 0
    levels = [ m500, km1, km3, km5, km8, surface ]
    
    # Begin: Loop through mandatory levels
    #-----------------------------------------------------------------------------
    for k in range( 0, len(levels) ):
        
        # Begin: Loop through pre-defined updraft diameter in the x-direction
        #-----------------------------------------------------------------------------
        for i in range( i_start, i_end ):
        
            # Begin: Loop through pre-defined updraft diameter in the y-direction
            #-----------------------------------------------------------------------------
            for j in range( j_start, j_end ):
            
                # Get meso values at current grid point
                w_current = w[ k, j, i ]
                zeta_current = zvort[ k, j, i ]
                w_zeta_current = w_zeta[ k, j, i ]
                thpert_current = thpert[j,i]
                dbz_current = dbz[j,i]
                
                # Determine if current gridpoint exceeds threshold or not
                if( k == 0 and w_current >= w_thresh ):
                    w_area_500m += 1
                if( k == 0 and zeta_current >= zeta_thresh ):
                    zeta_area_500m += 1
                if( k == 0 and w_zeta_current >= w_zeta_thresh ):
                    w_zeta_area_500m += 1
                    
                if( k == 1 and dbz_current >= dbz_thresh ):
                    dbz_area += 1
                if( k == 1 and w_current >= w_thresh ):
                    w_area_1km += 1
                if( k == 1 and zeta_current >= zeta_thresh ):
                    zeta_area_1km += 1
                if( k == 1 and w_zeta_current >= w_zeta_thresh ):
                    w_zeta_area_1km += 1    
                
                if( k == 2 and w_current >= w_thresh ):
                    w_area_3km += 1
                if( k == 2 and zeta_current >= zeta_thresh ):
                    zeta_area_3km += 1
                if( k == 2 and w_zeta_current >= w_zeta_thresh ):
                    w_zeta_area_3km += 1
            
                if( k == 3 and w_current >= w_thresh ):
                    w_area_5km += 1
                if( k == 3 and zeta_current >= zeta_thresh ):
                    zeta_area_5km += 1
                if( k == 3 and w_zeta_current >= w_zeta_thresh ):
                    w_zeta_area_5km += 1
                    
                if( k == 4 and w_current >= w_thresh ):
                    w_area_8km += 1
                if( k == 4 and zeta_current >= zeta_thresh ):
                    zeta_area_8km += 1
                if( k == 4 and w_zeta_current >= w_zeta_thresh ):
                    w_zeta_area_8km += 1
                    
                if( k == 5 and thpert_current <= -1.0 ):
                    thpert_surface += 1
                    
                
                    
            
            # End: Loop through pre-defined updraft diameter in the y-direction
            #-----------------------------------------------------------------------------
        
        # End: Loop through pre-defined updraft diameter in the x-direction
        #-----------------------------------------------------------------------------
        
    # End: Loop through mandatory levels
    #-----------------------------------------------------------------------------
    
    # Get the area of a single gridpont (m^2)
    area = dx * dy
    
    # Compute total area satisfying thresholds at each level & convert to km^2
    w_area_500m_total =  ( w_area_500m * area ) * (1e-6) 
    zeta_area_500m_total = ( zeta_area_500m * area ) * (1e-6) 
    w_zeta_area_500m_total = ( w_zeta_area_500m * area ) * (1e-6) 
    
    w_area_1km_total = ( w_area_1km * area ) * (1e-6) 
    zeta_area_1km_total = ( zeta_area_1km * area ) * (1e-6) 
    w_zeta_area_1km_total = ( w_zeta_area_1km * area ) * (1e-6) 
    
    w_area_3km_total = ( w_area_3km * area ) * (1e-6) 
    zeta_area_3km_total = ( zeta_area_3km * area ) * (1e-6) 
    w_zeta_area_3km_total = ( w_zeta_area_3km * area ) * (1e-6) 
    
    w_area_5km_total = ( w_area_5km * area ) * (1e-6) 
    zeta_area_5km_total = ( zeta_area_5km * area ) * (1e-6) 
    w_zeta_area_5km_total = ( w_zeta_area_5km * area ) * (1e-6)
    
    w_area_8km_total = ( w_area_8km * area ) * (1e-6) 
    zeta_area_8km_total = ( zeta_area_8km * area ) * (1e-6) 
    w_zeta_area_8km_total = ( w_zeta_area_8km * area ) * (1e-6)
    
    thpert_surface_total = ( thpert_surface * area ) * (1e-6)
    dbz_total = ( dbz_area * area ) * (1e-6)
    
    #-----------------------------------------------------------------------------
    # End: Compute Mesocyclone Area Estimates
    #-----------------------------------------------------------------------------


    #-----------------------------------------------------------------------------
    # Begin: Compute Mesocyclone Depth Estimates
    #-----------------------------------------------------------------------------
    
    # Assume meso is centered on uh_maxima and find all vertical gridpoints that satisfy meso threshold between 0-10 km
    meso = np.where( zvort[:km10, uh_y, uh_x] >= zeta_thresh )
    
    # Requires try/except clause to properly execute without error when meso_depth is zero
    try:
        # Grab the lower bound of the meso
        meso_lower = meso[0][0]
        
        # Grab the upper bound of the meso
        meso_upper = meso[0][-1]
        
        # Compute the mesocyclone depth at the uh_maxima location
        meso_depth = zh[ meso_upper ].m - zh[ meso_lower ].m
    
    except IndexError:
        meso_depth = 0.0
    
    #-----------------------------------------------------------------------------
    # End: Compute Mesocyclone Depth Estimates
    #-----------------------------------------------------------------------------


    #-----------------------------------------------------------------------------
    # Begin: Get minimum thpert value within mesocyclone box
    #-----------------------------------------------------------------------------
    
    if( idx > 1 ):
        thpert_surface_min = np.amin( thpert[ j_start:j_end, i_start:i_end ] )
    else:
        thpert_surface_min = 0.0
        
    print( 'CP Area: {} km^2 \tCP Intensity: {} K'.format( thpert_surface_total, thpert_surface_min ) )
    
    #-----------------------------------------------------------------------------
    # End: Get minimum thpert value within mesocyclone box
    #-----------------------------------------------------------------------------

    # Return area box indices and computed values
    return(
           i_start, i_end, j_start, j_end,
           w_area_500m_total, w_area_1km_total, w_area_3km_total, w_area_5km_total, w_area_8km_total,
           zeta_area_500m_total, zeta_area_1km_total, zeta_area_3km_total, zeta_area_5km_total, zeta_area_8km_total,
           w_zeta_area_500m_total, w_zeta_area_1km_total, w_zeta_area_3km_total, w_zeta_area_5km_total, w_zeta_area_8km_total,
           meso_depth, thpert_surface_total, thpert_surface_min, dbz_total
          )

#-----------------------------------------------------------------------------
# End: Define cm1_uh_tracker
#-----------------------------------------------------------------------------




