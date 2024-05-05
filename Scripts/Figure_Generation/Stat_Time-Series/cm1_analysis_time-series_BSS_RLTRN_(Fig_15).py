#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 21:05:04 2023

@author: roger
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:47:28 2023

@author: roger
"""

# Import external modules
#-------------------------------------------------------------------
import sys
import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import pandas as pd   
from matplotlib import cm
from scipy import stats
import scipy
from matplotlib import rcParams
import matplotlib.lines as mlines
import matplotlib.patheffects as PathEffects

# Begin Read-In Data Function
#-------------------------------------------------------------------

def read_stats_data( filename ):
    
    # Report program status to the terminal
    print( '\n\t Now opening input file: {} '.format( filename ) )
    
    # Store data columns
    ( 
     mode, time, x1, y1, x2, y2, zs, 
     w500m, w1km, w3km, w5km,
     zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
     uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
     mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
     shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
     srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
     lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
     zA500m, zA1km, zA3km, zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth,
     cp_intensity, cp_area                                                    
    ) = np.genfromtxt(
                      fname = filename,
                      delimiter = ',',
                      skip_header = 3,
                      usecols = np.arange( 0, 70, 1 ),
                      unpack = True
                    )
        
    # Report program status to terminal
    print( '\n\t Raw data successfully stored into dataframe structure...')

    # Return read in data arrays
    return ( 
            mode, time, x1, y1, x2, y2, zs, 
            w500m, w1km, w3km, w5km,
            zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
            uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
            mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
            shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
            srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
            lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
            zA500m, zA1km, zA3km, zA5km, zA8km, wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth,
            cp_intensity, cp_area                                                        
           )

# End Read-In Data Function
#-------------------------------------------------------------------


# Begin Moving Average Function
#-------------------------------------------------------------------

def gauss_smoother( arr, smooth = 6 ) :
    ret_arr = scipy.ndimage.gaussian_filter1d( arr, smooth )
    return( ret_arr )

# End Moving Average Function
#-------------------------------------------------------------------



#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Begin Main Script
#-------------------------------------------------------------------
#-------------------------------------------------------------------

# Store program name
program = 'CM1_Analysis_Cross_RLTRN_Variable-State_Time-Series.py'

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( '\nBegin Program: {}'.format ( program ) )


# Simulation-Specific Parameters
#------------------------------------------------------------------------------
sim = [ 'bss_nc_rltrn', 'bss_cs_rltrn' ]        # Sim names
color = [ 'tab:blue', 'tab:red'  ]              # Sim colors
sup_start = [ 43, 40 ]                          # Supercell Start Index
sup_end = [ 62, 55 ]                            # Supercell End Index
dis_end = [ 79, 58 ]                            # Dissipation Time Index
linear = [ 0, 0 ]                               # Upscale Growth Index
#------------------------------------------------------------------------------
        
        
# Smoothing Parameter
smooth = 2.0

# Line widths
lwidth = 3.0

# Set tick size (Must be before calling plot object)
plt.rcdefaults()
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)

landscape = False

# Figure orientation logic
if( landscape == False ):
    figure_mosaic = """
                    AABB
                    CCDD
                    EEFF
                    GGHH
                    """
    figsize = ( 18, 22.5 )
else:
    figure_mosaic = """
                    AABBCCDD
                    EEFFGGHH
                    """
    figsize = ( 35, 12.5 )
    
# Create Figure and axes object                
fig, axes =  plt.subplot_mosaic( mosaic = figure_mosaic, figsize = figsize, tight_layout = True )

# Create callable axes objects
ax_A = fig.add_subplot( axes['A'] )
ax_B = fig.add_subplot( axes['B'] )
ax_C = fig.add_subplot( axes['C'] )
ax_D = fig.add_subplot( axes['D'] )
ax_D2 = ax_D.twinx()
ax_E = fig.add_subplot( axes['E'] )
ax_F = fig.add_subplot( axes['F'] )
ax_G = fig.add_subplot( axes['G'] )
ax_H = fig.add_subplot( axes['H'] )
ax_H2 = ax_H.twinx()

# Set titles
ax_A.set_title( label = 'a: Mesocyclone Intensity', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_B.set_title( label = 'b: Mesocyclone Depth', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_C.set_title( label = 'c: Mesocyclone Area', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_D.set_title( label = 'd: Inflow Thermodynamics', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_E.set_title( label = 'e: Inflow Kinematics', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_F.set_title( label = 'f: Inflow Storm-Relative Helicity', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_G.set_title( label = 'g: Inflow Sounding Surface Elevation', loc = 'left', fontsize = 18, fontweight = 'bold' )
ax_H.set_title( label = 'h: Cold Pool Metrics', loc = 'left', fontsize = 18, fontweight = 'bold' )


# Set axes labels
ax_G.set_xlabel( 'Integration Time after CI (min.)', fontsize = 18, fontweight = 'bold'  )
ax_H.set_xlabel( 'Integration Time after CI (min.)', fontsize = 18, fontweight = 'bold'  )
ax_A.set_ylabel( '2-5 km UH (m$^{2}$ s$^{-2}$)', fontsize = 18, fontweight = 'bold'  )
ax_B.set_ylabel( 'Depth (km)', fontsize = 18, fontweight = 'bold'  )
ax_C.set_ylabel( 'w\u03b6 > 0.1 m s$^{-2}$ Area (km$^{2}$)', fontsize = 18, fontweight = 'bold'  )
ax_D.set_ylabel( 'MLCAPE (J kg$^{-1}$)', fontsize = 18, fontweight = 'bold'  )
ax_D2.set_ylabel( 'MLCIN (J kg$^{-1}$)', fontsize = 18, fontweight = 'bold'  )
ax_E.set_ylabel( 'Shear (m s$^{-1}$)', fontsize = 18, fontweight = 'bold'  )
ax_F.set_ylabel( 'SRH (m$^{2}$ s$^{-2}$)', fontsize = 18, fontweight = 'bold'  )
ax_G.set_ylabel( 'Z (m)', fontsize = 18, fontweight = 'bold'  )
ax_H2.set_ylabel( 'Min $\u03b8_{pert}$ (K)', fontsize = 18, fontweight = 'bold'  )
ax_H.set_ylabel( '$\u03b8_{pert}$ < -1 K Area (km$^{2}$)', fontsize = 18, fontweight = 'bold' )


# Special line
# ax_A.axhline( y = 150.0, color = 'k' ) 
ax_H.axhline( y = 0, color = 'k' ) 

ax_A.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_A.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_A.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_B.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_B.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_B.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_C.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_C.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_C.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_D.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_D.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_D.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_E.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_E.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_E.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_F.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_F.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_F.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_G.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_G.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_G.axvline( x = 360.0, color = 'k', linestyle = '--' )

ax_H.axvline( x = 180.0, color = 'k', linestyle = '--' )
ax_H.axvline( x = 300.0, color = 'k', linestyle = '--' )
ax_H.axvline( x = 360.0, color = 'k', linestyle = '--' )


txta_1 = ax_A.annotate( 'BSS0', xy = (110, 475), fontsize = 12, fontweight = 'bold' )
txta_2 = ax_A.annotate( 'BSS0 -> BSS1', xy = (215, 475), fontsize = 12, fontweight = 'bold' )
txta_3 = ax_A.annotate( 'BSS1 -> BSS2', xy = (301.5, 475), fontsize = 12, fontweight = 'bold' )
txta_4 = ax_A.annotate( 'BSS_Start', xy = (177.5, 375), rotation = 90, fontsize = 14, fontweight = 'bold' )
txta_5 = ax_A.annotate( 'BSS1', xy = (297.5, 400), rotation = 90, fontsize = 14, fontweight = 'bold' )
txta_6 = ax_A.annotate( 'BSS_End', xy = (365, 400), rotation = 90, fontsize = 14, fontweight = 'bold' )

txtb_1 = ax_B.annotate( 'BSS0', xy = (110, 9.75), fontsize = 12, fontweight = 'bold' )
txtb_2 = ax_B.annotate( 'BSS0 -> BSS1', xy = (215, 9.75), fontsize = 12, fontweight = 'bold' )
txtb_3 = ax_B.annotate( 'BSS1 -> BSS2', xy = (301.5, 9.75), fontsize = 12, fontweight = 'bold' )
txtb_4 = ax_B.annotate( 'BSS_Start', xy = (177.5, 7.75), rotation = 90, fontsize = 14, fontweight = 'bold' )
txtb_5 = ax_B.annotate( 'BSS1', xy = (297.5, 8), rotation = 90, fontsize = 14, fontweight = 'bold' )
txtb_6 = ax_B.annotate( 'BSS_End', xy = (365, 8), rotation = 90, fontsize = 14, fontweight = 'bold' )

# Create a halo around North Arrow
plt.setp( [txta_1, txta_2, txta_3, txta_4, txta_5, txta_6, txtb_1, txtb_2, txtb_3, txtb_4, txtb_5, txtb_6 ], path_effects =[ PathEffects.withStroke( linewidth = 6, foreground = 'w', alpha = 0.75 ) ] )


# Fixed y-limits
# Fixed y-limits
ax_A.set_ylim( -10, 500 )
ax_A.set_xlim( 60, 380 )
ax_B.set_ylim( -0.25, 10.25 )
ax_B.set_xlim( 60, 380  )
ax_C.set_ylim( -0.25, 11.25  )
ax_C.set_xlim( 60, 380 )
ax_D.set_ylim( -25, 1500 )
ax_D2.set_ylim( 0, 200 )
ax_D.set_xlim( 60, 380 )
ax_E.set_ylim( 10, 40 )
ax_E.set_xlim( 60, 380 )
ax_F.set_ylim( -10, 425  )
ax_F.set_xlim( 60, 380 )
ax_G.set_ylim( -50, 1200 )
ax_G.set_xlim( 60, 380 )
ax_H2.set_ylim( -8, 0 )
ax_H.set_ylim( -25, 800 )
ax_H.set_xlim( 60, 380 )

# Begin: Loop through each simulation
#------------------------------------------------------------------------------
for i in range( 0, len( sim ) ):
    
    # Construct filename string (!!! Requires appropriate directory containing BSS IDTRN CSVs !!!)
    os_type = 0
    if( os_type == 0 ):
        wdir = '/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR_Modeling_Project/Simulations/Variable_State_Sims/Stats_Spreadsheets/'
    else:
        wdir = 'C:/Users/rriggin/OneDrive - University of North Carolina at Charlotte/CSTAR_Modeling_Project/Simulations/Variable_State_Sims/Stats_Spreadsheets/'
    filename = wdir + str(sim[i]) + '_model_output_stats.csv'
    
    
    # Read in current dataset
    ( mode, time, x1, y1, x2, y2, zs, w500m, w1km, w3km, w5km,
      zvort_surface, zvort_500m, zvort_1km, zvort_3km, zvort_5km,
      uh_max, uh_area, uhx, uhy, uhi, uhj, sx, sy, sxi, syj,
      mean_uv, mean_uv_dir, rm, rm_dir, lm, lm_dir,
      shear1km, shear1km_dir, shear3km, shear3km_dir, shear6km, shear6km_dir,
      srh500m, srh1km, srh3km, sbcape, sbcin, mlcape, mlcin, mucape, mucin, cape3,
      lcl, lfc, delta_lfc_lcl, el, wA500m, wA1km, wA3km, wA5km, wA8km,
      zA500m, zA1km, zA3km, zA5km, zA8km,
      wzA500m, wzA1km, wzA3km, wzA5km, wzA8km, meso_depth,
      cp_intensity, cp_area ) = read_stats_data( filename )
    
    # Run appropriate parameters through guassian smoother
    time = gauss_smoother( time, smooth )
    
    # Correct indices for 2 hr CI delay with RLTRN
    start_time = 0
    time = time - 120.0
    sup_start[i] = sup_start[i] - 24
    sup_end[i] = sup_end[i] - 24
    dis_end[i] = dis_end[i] - 24
    if( linear[i] > 0.0 ):
        linear[i] = linear[i] - 24

    
    uh_max = gauss_smoother( uh_max, smooth )
    uh_area = gauss_smoother( uh_area, smooth )
    
    meso_depth = gauss_smoother( meso_depth, smooth )
    
    zA3km = gauss_smoother( wzA3km, smooth )
    zA5km = gauss_smoother( wzA5km, smooth )
    
    sbcape = gauss_smoother( sbcape, smooth )
    sbcin = gauss_smoother( sbcin, smooth )
    mlcape = gauss_smoother( mlcape, smooth )
    mlcin = gauss_smoother( -1 * mlcin, smooth )
    
    # shear1km = gauss_smoother( shear1km, smooth )
    shear3km = gauss_smoother( shear3km, smooth )
    shear6km = gauss_smoother( shear6km, smooth)
    
    # srh500m = gauss_smoother( srh500m, smooth )
    srh1km = gauss_smoother( srh1km, smooth )
    srh3km = gauss_smoother( srh3km, smooth )
    
    cp_intensity = gauss_smoother( cp_intensity, smooth)
    
    # cp_area = cp_area / 100.0
    cp_area = gauss_smoother( cp_area, smooth )

    # Run terrain profile through guassian smoother    
    zs_alt = gauss_smoother( zs * 1000, smooth )
    
    # Find peak altitude location
    peak = np.where( zs_alt == np.amax( zs_alt[ sup_start[i]:dis_end[i] ] ) )
    peak = peak[0][0] 

    for ax in axes:
        axes[ax].axvline( x = time[peak], linestyle = '-', linewidth = 1.5*lwidth, color = color[i], alpha = 0.5 )

    zline = mlines.Line2D( [], [], color = 'k', linewidth = 1.5*lwidth, linestyle = '-', label = 't = Peak Elevation', alpha = 0.5 )
    
    # Begin Plotting Data
    #--------------------------------------------------------------------------
    
    # Panel A Line Plot
    ax_A.plot( time[ start_time:dis_end[i] ], uh_area[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-',label = sim[i].upper() )
    
    # Demarcate start/end supercell mode
    ax_A.scatter( time[ sup_start[i] ], uh_area[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
    ax_A.scatter( time[ sup_end[i] ], uh_area[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

    # Demarcate start/end linear mode
    if( linear[i] > 0 ):
        ax_A.scatter( time[ sup_end[i] + 1 ], uh_area[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        ax_A.scatter( time[ linear[i] ], uh_area[ linear[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )

    # Add legend & Grid
    ax_A.legend( title = 'Simulations', prop = {'size': 14}, loc = 'lower right', facecolor = 'lightgrey' )
    ax_A.grid()



    # Same for Panel B
    ax_B.plot( time[ start_time:dis_end[i] ], meso_depth[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper() )
    
    ax_B.scatter( time[ sup_start[i] ], meso_depth[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
    ax_B.scatter( time[ sup_end[i] ], meso_depth[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
    
    if( linear[i] > 0 ):
        ax_B.scatter( time[ sup_end[i] + 1 ], meso_depth[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        ax_B.scatter( time[ linear[i] ], meso_depth[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10 )
        
    ax_B.legend( title = 'Simulations', prop = {'size': 14}, loc = 'lower right', facecolor = 'lightgrey' )
    ax_B.grid()
    
    
    
    # Same for Panel C
    ax_C.plot( time[ start_time:dis_end[i] ], zA3km[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '--', label = sim[i].upper() + ' \u03b6-3km' )
    ax_C.plot( time[ start_time:dis_end[i] ], zA5km[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper() + ' \u03b6-5km' )
    
    ax_C.scatter( time[ sup_start[i] ], zA3km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_C.scatter( time[ sup_end[i] ], zA3km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_C.scatter( time[ sup_start[i] ], zA5km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_C.scatter( time[ sup_end[i] ], zA5km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    if( linear[i] > 0 ):
        ax_C.scatter( time[ sup_end[i] + 1 ], zA3km[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_C.scatter( time[ linear[i] ], zA3km[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_C.scatter( time[ sup_end[i] + 1 ], zA5km[ sup_end[i] + 1 ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_C.scatter( time[ linear[i] ], zA5km[ linear[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    # Create artificial lines for custom legend
    lineC1 = mlines.Line2D( [], [], color = 'k', linestyle = '--', label = ' w\u03b6-3km' )
    lineC2 = mlines.Line2D( [], [], color = 'k', linestyle = '-', label = ' w\u03b6-5km' )
    
    ax_C.legend( handles = [lineC1, lineC2], prop = {'size': 14}, loc = 'upper right', facecolor = 'lightgrey' )
    ax_C.grid()
    
    
    
    # Same for Panel D
    ax_D2.plot( time[ start_time:dis_end[i] ], mlcin[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '--', label = sim[i].upper() + ' MLCIN' )
    ax_D.plot( time[ start_time:dis_end[i] ], mlcape[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper() + ' MLCAPE' )
    
    ax_D2.scatter( time[ sup_start[i] ], mlcin[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_D2.scatter( time[ sup_end[i] ], mlcin[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_D.scatter( time[ sup_start[i] ], mlcape[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_D.scatter( time[ sup_end[i] ], mlcape[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    if( linear[i] > 0 ):
        ax_D2.scatter( time[ sup_end[i] + 1 ], mlcin[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_D2.scatter( time[ linear[i] ], mlcin[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_D.scatter( time[ sup_end[i] + 1 ], mlcape[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_D.scatter( time[ linear[i] ], mlcape[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    lineD1 = mlines.Line2D( [], [], color = 'k', linestyle = '--', label = ' MLCIN' )
    lineD2 = mlines.Line2D( [], [], color = 'k', linestyle = '-', label = ' MLCAPE' )
    
    ax_D.legend( handles = [lineD1, lineD2], prop = {'size': 14}, loc = 'upper right', facecolor = 'lightgrey' )
    ax_D.grid()
    ax_D2.grid( alpha = 0.25 )
    
    
    
    # Same for Panel E
    # ax_E.plot( time[ 20:dis_end[i] ], shear1km[ 20:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = ':', label = sim[i].upper() + ' 0-1km Shear' ) 
    ax_E.plot( time[ start_time:dis_end[i] ], shear3km[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '--', label = sim[i].upper() + ' 0-3km Shear' )
    ax_E.plot( time[ start_time:dis_end[i] ], shear6km[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper() + ' 0-6km Shear' )
    
    # ax_E.scatter( time[ sup_start[i] ], shear1km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    # ax_E.scatter( time[ sup_end[i] ], shear1km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_E.scatter( time[ sup_start[i] ], shear3km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_E.scatter( time[ sup_end[i] ], shear3km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_E.scatter( time[ sup_start[i] ], shear6km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_E.scatter( time[ sup_end[i] ], shear6km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    if( linear[i] > 0 ):
        # ax_E.scatter( time[ sup_end[i] + 1 ], shear1km[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        # ax_E.scatter( time[ linear[i] ], shear1km[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_E.scatter( time[ sup_end[i] + 1 ], shear3km[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_E.scatter( time[ linear[i] ], shear3km[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_E.scatter( time[ sup_end[i] + 1 ], shear6km[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_E.scatter( time[ linear[i] ], shear6km[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    # lineE1 = mlines.Line2D( [], [], color = 'k', linestyle = ':', label = ' 0-1 km Shear' )
    lineE2 = mlines.Line2D( [], [], color = 'k', linestyle = '--', label = ' 0-3 km Shear' )
    lineE3 = mlines.Line2D( [], [], color = 'k', linestyle = '-', label = ' 0-6 km Shear' )
    
    ax_E.legend( handles = [ lineE2, lineE3 ], prop = {'size': 14}, loc = 'upper right', facecolor = 'lightgrey' )
    ax_E.grid()
    
    
    
    # Same for Panel F
    # ax_F.plot( time[ 20:dis_end[i] ], srh500m[ 20:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = ':', label = sim[i].upper() + ' 0-500m SRH' )
    ax_F.plot( time[ start_time:dis_end[i] ], srh1km[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '--', label = sim[i].upper() + ' 0-1km SRH' )
    ax_F.plot( time[ start_time:dis_end[i] ], srh3km[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper() + ' 0-3km SRH' )
    
    # ax_F.scatter( time[ sup_start[i] ], srh500m[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    # ax_F.scatter( time[ sup_end[i] ], srh500m[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_F.scatter( time[ sup_start[i] ], srh1km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_F.scatter( time[ sup_end[i] ], srh1km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_F.scatter( time[ sup_start[i] ], srh3km[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_F.scatter( time[ sup_end[i] ], srh3km[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    if( linear[i] > 0 ):
        # ax_F.scatter( time[ sup_end[i] + 1 ], srh500m[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        # ax_F.scatter( time[ linear[i] ], srh500m[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_F.scatter( time[ sup_end[i] + 1 ], srh1km[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_F.scatter( time[ linear[i] ], srh1km[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_F.scatter( time[ sup_end[i] + 1 ], srh3km[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_F.scatter( time[ linear[i] ], srh3km[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    # lineF1 = mlines.Line2D( [], [], color = 'k', linestyle = ':', label = ' 0-500 m SRH' )
    lineF2 = mlines.Line2D( [], [], color = 'k', linestyle = '--', label = ' 0-1 km SRH' )
    lineF3 = mlines.Line2D( [], [], color = 'k', linestyle = '-', label = ' 0-3 km SRH' )
    
    ax_F.legend( handles = [ lineF2, lineF3 ], prop = {'size': 14}, loc = 'upper right', facecolor = 'lightgrey' )
    ax_F.grid()
    
    
    # # Run terrain profile through guassian smoother   
    # zs_alt = gauss_smoother( zs * 1000, smooth )
    
    # Same plotting methods for panel G
    ax_G.plot( time[ start_time:dis_end[i] ], zs_alt[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper()  )
    
    ax_G.scatter( time[ sup_start[i] ], zs_alt[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_G.scatter( time[ sup_end[i] ], zs_alt[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    if( linear[i] > 0 ):
        ax_G.scatter( time[ sup_end[i] + 1 ], zs_alt[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_G.scatter( time[ linear[i] ], zs_alt[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
    ax_G.grid()
    
    
    
    # Same for panel H
    ax_H2.plot( time[ start_time:dis_end[i] ], cp_intensity[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '--', label = sim[i].upper() + '$\u03b8_{pert}$ Min ' )
    ax_H.plot( time[ start_time:dis_end[i] ], cp_area[ start_time:dis_end[i] ], color = color[i], linewidth = lwidth, linestyle = '-', label = sim[i].upper() + ' $\u03b8_{pert}$ Area' )
 
    
    ax_H2.scatter( time[ sup_start[i] ], cp_intensity[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_H2.scatter( time[ sup_end[i] ], cp_intensity[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    
    ax_H.scatter( time[ sup_start[i] ], cp_area[ sup_start[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
    ax_H.scatter( time[ sup_end[i] ], cp_area[ sup_end[i] ], marker = 'v', color = color[i], edgecolor = 'k', s = 100,  zorder = 10  )
    
    if( linear[i] > 0 ):
        ax_H2.scatter( time[ sup_end[i] + 1 ], cp_intensity[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_H2.scatter( time[ linear[i] ], cp_intensity[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        
        ax_H.scatter( time[ sup_end[i] + 1 ], cp_area[ sup_end[i] + 1 ], marker = 's', color = color[i], edgecolor = 'k', s = 100, zorder = 10  )
        ax_H.scatter( time[ linear[i] ], cp_area[ linear[i] ], marker = 's', color = color[i], edgecolor = 'k', s = 100,  zorder = 10  )
    
    lineH1 = mlines.Line2D( [], [], color = 'k', linestyle = '-', label = ' $\u03b8_{pert}$ Area' )
    lineH2 = mlines.Line2D( [], [], color = 'k', linestyle = '--', label = ' $\u03b8_{pert}$ Min' )
    
    ax_H.legend( handles = [lineH1, lineH2], prop = {'size': 14}, loc = 'upper right', facecolor = 'lightgrey' )
    ax_H.grid()
    ax_H2.grid( alpha = 0.25 )

    # End Plotting Data
    #--------------------------------------------------------------------------

# End: Loop through each simulation
#------------------------------------------------------------------------------

ax_A.grid()
ax_B.grid()
ax_C.grid()
ax_D.grid()
ax_E.grid()
ax_F.grid()
ax_G.grid()
ax_H.grid()

# Create storm mode custom legend for panel G    
sup_lab = mlines.Line2D( [], [], color = 'k', marker = 'v', ms = 12, label = 'Supercell' )
linear_lab = mlines.Line2D( [], [], color = 'k', marker = 's', ms = 12, label = 'Linear' )
g_leg1 = ax_G.legend( title = 'Storm Mode (Start/End)', handles = [sup_lab, linear_lab, zline], loc = 'upper right', prop = {'size': 14}, facecolor = 'lightgrey' )
# g_leg2 = ax_G.legend( title = 'Simulations', prop = {'size': 14}, loc = 'upper left', facecolor = 'lightgrey' )
ax_G.add_artist( g_leg1 )
# ax_G.add_artist( g_leg2 )

# Save figure
fig.savefig(
            fname = 'BSS_RLTRN_Time-Series.jpeg',
            dpi = 300,
            bbox_inches = "tight"
            )


# Report script runtime
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )

