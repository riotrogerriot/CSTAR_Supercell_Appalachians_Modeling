#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 22:17:07 2021

@author: roger
"""

#-----------------------------------------------------------------------------
# composites_2_cm1.py
#
# Roger Riggin          # UNCC CSTAR Research Assistantship
#
#    Latest Update: 2/25/2021
# 
# Program Summary:
#
#   This script opens an excel spreadsheet containing the composite sounding 
#   data and reformats it to match the required format to use as a CM1 initial
#   conditions sounding.
#-----------------------------------------------------------------------------


# Import external modules
#-----------------------------------------------------------------------------
import os
import pandas as pd
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units


# Set Working Directory
#-----------------------------------------------------------------------------

# Path to my CM1 directory
composite_dir = "/Users/roger/Desktop/School/CSTAR/CM1/Init_Sounding"

# Change the working directory to the provided path
os.chdir( composite_dir )
    
# Retrieve the current directory
cwd = os.getcwd()

# Report working directory to user
print( '\nCurrent Working Directory: \n {}'.format(cwd) )


#-----------------------------------------------------------------------------
# Begin main script
#-----------------------------------------------------------------------------

# Excel Filename
filename = "all_composites.xlsx"

# Store sheet names in a list for looping
sheet = [ "Up-Non", "Up-Cross", "Peak-Non", 
          "Peak-Cross", "Down-Non", "Down-Cross" ]

# Index list for relevant non-crossing parameters
nc_index = [ 0, 37, 74, 111, 148, 185, 222, 259 ]

# Index list for relevant Up-Cross & Peak-Cross parameters
c_index12 = [ 0, 24, 48, 72, 96, 120, 144, 168 ]

# Index list for relevant Down-Cross parameters
c_index3 = [ 0, 23, 46, 69, 92, 115, 138, 161 ]


# Loop through each excel sheet
#------------------------------
for j in range( len( sheet ) ):
    
    # Note: Cross sheets are odd and Non-Cross are Even indices in sheet var    
    
    # Begin indice list logic
    #------------------------------
    
    # Select non-cross indice list when j is even
    if ( j == 0 or j % 2 == 0 ) :
        index = nc_index
    
    # Down-Cross indices do not match up & peak!    
    elif ( j == 5 ):
        
        # Select proper indice list for Down-Cross sheet
        index = c_index3
            
    # Select crosser index when j is odd  
    else:
        index = c_index12

    # End indice list logic
    #------------------------------
    
    
    # Report to user which sounding is being converted
    print( "\nBegin converting '" + sheet[j] + "' to CM1 compatible input sounding...")
    
    # Open the excel file with pandas
    data = pd.read_excel( filename, sheet_name = sheet[j] )
    
    # Pull relevant sounding data out and convert to np array
    composite = pd.DataFrame( data.iloc[ :, index ]  ).to_numpy()
    
    
    # Begin Surface Data
    #-----------------------------------------------------------------------------
    
    # Pull surface info from composite array
    surface = composite[1, 1:4]
    
    # Subset and add units for future metpy computations
    surfaceP  = units.Quantity( surface[0], "hPa"   )                 
    surfaceT  = units.Quantity( surface[1], "degC"  )                    
    surfaceTd = units.Quantity( surface[2], "degC"  )
    
    # Compute surface theta
    surface_th = mpcalc.potential_temperature( 
                                              pressure = surfaceP,
                                              temperature = surfaceT 
                                             )
    
    # Report surface theta
    #print( surface_th.magnitude )
    
    # Compute surface RH
    surface_rh = mpcalc.relative_humidity_from_dewpoint( 
                                                        temperature = surfaceT,
                                                        dewpoint = surfaceTd 
                                                       )
    # Report surface RH
    #print( surface_rh.magnitude )
    
    # Use surface RH to compute surface Qv
    surface_qv = mpcalc.mixing_ratio_from_relative_humidity( 
                                                            pressure = surfaceP,
                                                            temperature = surfaceT,
                                                            relative_humidity = surface_rh
                                                              ).to('g/kg')
    
    # Report surface Qv
    #print( surface_qv.magnitude )
    
    # End Surface Data
    #-----------------------------------------------------------------------------
    
    
    # Begin Altitude Data
    #-----------------------------------------------------------------------------
    
    # Pull out altitude data from composites array 
    # Note: Start w/ 2 so we don't repeat surface obs
    alt   = composite[ 2:, 0 ]
    temp  = composite[ 2:, 2 ] 
    td    = composite[ 2:, 3 ] 
    press = composite[ 2:, 1 ] 
    u     = composite[ 2:, 5 ] 
    v     = composite[ 2:, 6 ] 
    
    # Create an empty array for computed data
    th = units.Quantity( np.zeros( len( temp ) ), "kelvin"        )      
    rh = units.Quantity( np.zeros( len( td   ) ), "dimensionless" )
    qv = units.Quantity( np.zeros( len( rh   ) ), "dimensionless" )
    
    # Add units to pressure, temp, and dewpoint
    temp  = units.Quantity( temp,  "degC" )
    td    = units.Quantity( td,    "degC" )
    press = units.Quantity( press, "hPa"  )
    
    # End Altitude data
    #-----------------------------------------------------------------------------
    

    # Create the textfile for CM1 input sounding to be written to
    init_sounding = open( r"CM1_Init_Sounding_" + sheet[j] + ".txt", "w+" )
    
    # Format each surface variable to CM1 format
    surfaceP =   round( surfaceP.magnitude,   4 )
    surface_th = round( surface_th.magnitude, 4 )
    surface_qv = round( surface_qv.magnitude, 4 )
    
    # Place variables for first string into a list
    surface_str = [ surfaceP, surface_th, surface_qv ] 
    
    # Construct and print the formatted first line
    # print( 
    #       "\n{:>12.4f}\t{:>12.4f}\t{:>12.4f}".format(
    #                             surface_str[0],
    #                             surface_str[1],
    #                             surface_str[2]
    #                            )
    #       )
    

    # Begin loop for alt data
    #-----------------------------------
    
    # Loop through the entire length of the dataset
    for i in range( 0, len( temp ) ):
        
        # Compute theta from current alt temp
        th[i] = mpcalc.potential_temperature( 
                                             pressure = press[i],
                                             temperature = temp[i] 
                                            )
        
        
        # Compute relative humidity from current alt temp & td
        rh[i] = mpcalc.relative_humidity_from_dewpoint( 
                                                       temperature = temp[i],
                                                       dewpoint = td[i]
                                                      )
        
        
        # Compute mixing-ratio from RH at current alt
        qv[i] = mpcalc.mixing_ratio_from_relative_humidity( 
                                                           pressure = press[i],
                                                           temperature = temp[i],
                                                           relative_humidity = rh[i]
                                                          )
        
        # Place alt data for current line into a list
        alt_str = [ 
                   alt[i],                              # Altitude (m)
                    th[i].magnitude,                    # Theta (K)
                    qv[i].to("g/kg").magnitude,         # Mixing-Ratio (g/kg)
                     u[i],                              # U-comp (m/s)
                     v[i]                               # V-comp (m/s)
                  ]
        
        # Construct and print current line using format descriptors
        # print( 
        #       "\n{:>12.0f}\t{:>12.4f}\t{:>12.4f}\t{:>12.4f}\t{:>12.4f}"
        #           .format(
        #                   alt_str[0],
        #                   alt_str[1],
        #                   alt_str[2],
        #                   alt_str[3],
        #                   alt_str[4] 
        #                   )
        #       )
        
        if ( i == 0 ):
            # Write the first line (Surface) to the file
            init_sounding.write( 
                                " {:>9.4f}\t{:>9.4f}\t{:>9.4f}".format(
                                                      surface_str[0],
                                                      surface_str[1],
                                                      surface_str[2]
                                                     )
                                  )
                
        # Write each alt line to the sounding file
        init_sounding.write( 
                            "\n {:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}"
                                .format(
                                        alt_str[0],
                                        alt_str[1],
                                        alt_str[2],
                                        alt_str[3],
                                        alt_str[4] 
                                       )
                           )
        
    # Report conversion of current file is complete
    print( "\nConversion of '" + sheet[j] + "' to CM1 compatible input sounding is complete! ")
    
    # End alt data loop
    #-----------------------------------
    
    
    # Create estimated data between 15-20 km 
    #---------------------------------------
    
    # Boolean logic to add artifical data for 20 km depth
    increase_depth = True
    
    # Run this code if increase depth is turned on
    if( increase_depth == True ):
        
        # Create an array containing every 100 m interval between 15-20 km
        add_alt = np.arange( 15100.0, 20100.0, 100.0 )
        
        # Pull out the 15 km theta value
        add_th = th[i].magnitude
        
        # Loop through the additional UL depth
        for j in range( 0, len( add_alt )  ):
            
            # Increase theta by 10 K per km ( Writing data for every 100 m so th + 1 )
            add_th = add_th + 1.0
            
            # Place alt data for current line into a list
            alt_str = [ 
                       add_alt[j],                       # Altitude (m)
                       add_th,                           # Theta (K) 
                       qv[i].to("g/kg").magnitude,       # Mixing-Ratio (g/kg) (Assumed Constant)
                       u[i],                             # U-comp (m/s) (Assumed Constant)
                       v[i]                              # V-comp (m/s) (Assumed Constant)
                      ]
            
            # Write each alt line to the sounding file
            init_sounding.write( 
                                "\n {:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}"
                                    .format(
                                            alt_str[0],
                                            alt_str[1],
                                            alt_str[2],
                                            alt_str[3],
                                            alt_str[4] 
                                           )
                               )
            
    # End estimated data between 15-20 km 
    #---------------------------------------
    
    
    # Close the newly written sounding text file
    init_sounding.close

# End excel sheet loop
#------------------------------


# End of program flag
print( "\nAll conversion have been completed... Goodbye!" )