#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------------
# 
#   Script Name: corrected_sounding_to_20km.py
# 
#   Roger Riggin    CSTAR Research Assistantship
#
#      Last Update: 6/16/2021
#
#   Program Summary:
#      
#
#
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Import external libraries
#-----------------------------------------------------------------------------
import os
import numpy as np
import pandas as pd
from datetime import datetime



#-----------------------------------------------------------------------------
# Custom Functions
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# Begin Main Program
#-----------------------------------------------------------------------------

# Record Script start time
startTime = datetime.now()

# Report program status to the terminal
print( '\nBegin program...')



#-----------------------------------------------------------------------------
# Set Working Directory
#-----------------------------------------------------------------------------

# Path to soundings
directory = "/Users/roger/Desktop/School/CSTAR/CM1/Init_Sounding/Eastin-Corrected"

# Path to Eastin's orginial files
direct = directory + "/Eastin-Originals"

# Report to terminal
print( "\nChanging working directory to original sounding location..." )

# Change the working directory to the provided path
os.chdir( direct )
    
# Retrieve the current directory
cwd = os.getcwd()

# Report working directory to user
print('\nCurrent Working Directory: \n {}'.format( cwd ) )

# Begin file name loops
#----------------------

# Store each file name in a list
files = [ name for name in os.listdir( '.' ) if os.path.isfile( name ) ]

# Report program status to the terminal
print( "\nReporting list of file names in cwd: " )

# Print out the contents of the file name list 
for i in range( len( files ) ):
    print( "\tFile #{}: {}".format( i + 1 , files[i] ) )

# End file name loop
#----------------------



#-----------------------------------------------------------------------------
# Open the input files and read in the data
#-----------------------------------------------------------------------------


# Begin Sounding Loop
#-----------------------

for j in range( 0, len(files) ):
        
    # String containing input file name and type
    infile = files[j]
    
    # Report program status to the terminal
    print( '\n\nNow opening input file: {} '.format( infile ) )
    
    # Open the sounding text file in read-only mode
    file = open( infile, 'r' )
    
    # Report program status to the terminal
    print( '\nPulling Surface data from input file: {} '.format( infile ) )
    
    # Pull the surface data from the first line in the input sounding
    ( sfc_press, sfc_th, sfc_qv ) = np.genfromtxt( 
                                                  fname = file,         # Input sounding textfile
                                                  delimiter = "\t",     # Column delimiting type
                                                  skip_header = 0,      # Num of lines in header info to skip
                                                  skip_footer = 150,    # Num of lines in footer info to skip
                                                  unpack = True         # Force data to be returned as structure arrays
                                                  )
    
    # Close the input sounding file
    file.close()
    
    
    # Re-open the sounding text file in read-only mode
    file = open( infile, 'r' )
    
    # Report program status to the terminal
    print( '\nPulling Height data from input file: {}... '.format( infile ) )
    
    # Read the height data in as a pandas data frame
    df = pd.read_csv( file, delimiter = "\t", header = None, error_bad_lines=False, skiprows=(1) )
    
    # Report program status to the terminal
    print( '\nConverting height data into numpy array...' )
    
    # Convert the dataframe into a numpy array
    df = df.to_numpy()
    
    # Report program status to the terminal
    print( '\nStoring height data into arrays...' )
    
    # Store height variables in separate arrays
    hgt = df[:,0]               # Height (m)
    th  = df[:,1]               # Theta (K)
    qv  = df[:,2]               # Mixing-Ratio (g/kg)
    u   = df[:,3]               # U-comp (m/s)
    v   = df[:,4]               # V-comp (m/s)
    
    # Report program status to the terminal
    print( '\nClosing input file: {}... '.format( infile ) )
    
    # Close the input sounding file
    file.close()
    
    # Report program status to the terminal
    print( '\nData successfully read from input file: {}'.format( infile ) )


#-----------------------------------------------------------------------------
# Construct new sounding with additional height data
#-----------------------------------------------------------------------------


    # Logic to create new directory if needed
    #----------------------------------------
    
    # Create a new directory for plots if it does not exist
    if( os.path.isdir( directory + "/Modified-Soundings" ) == False ):
        
        # Report to terminal
        print( "\nCreating directory to store modified soundings..." )
        
        # Create directory
        plot_dir = os.mkdir( directory + "/Modified-Soundings" )
        
    # End logic for new directory creation
    #----------------------------------------
        
    
    # Report to terminal
    print( "\nChanging working directory to modified sounding location..." )
    
    # Change the working directory to the provided path
    os.chdir( directory + "/Modified-Soundings" )
        
    # Retrieve the current directory
    cwd = os.getcwd()
    
    # Report working directory to user
    print('\nCurrent Working Directory: \n {}'.format( cwd ) )
    
    # Create the textfile for CM1 input sounding to be written to
    mod_sounding = open( files[j], "w+" )
    
    # Place variables for first string into a list
    surface_str = [ sfc_press, sfc_th, sfc_qv ] 
    
    # Report program status to terminal
    print( "\nBegin writing modified sounding data to {}...".format( files[j] ) )
    
    # Report contents of the surface data to terminal
    print( 
          "\n {:>12.4f}\t{:>12.4f}\t{:>12.4f}"
                  .format(
                          surface_str[0],
                          surface_str[1],
                          surface_str[2],
                          )
         )
    
    # Sounding construction loop
    #-------------------------------
    
    for i in range( 0, len( hgt ) ):
        
        # Place alt data for current line into a list
        alt_str = [ 
                   hgt[i],     # Altitude (m)
                    th[i],     # Theta (K)
                    qv[i],     # Mixing-Ratio (g/kg)
                     u[i],     # U-comp (m/s)
                     v[i]      # V-comp (m/s)
                  ]
        
        # Construct and print current line using format descriptors
        print( 
              " {:>12.4f}\t{:>12.4f}\t{:>12.4f}\t{:>12.4f}\t{:>12.4f}"
                  .format(
                          alt_str[0],
                          alt_str[1],
                          alt_str[2],
                          alt_str[3],
                          alt_str[4] 
                          )
              )
        
        if ( i == 0 ):
            # Write the first line (Surface) to the file
            mod_sounding.write( 
                                "{:>9.4f}\t{:>9.4f}\t{:>9.4f}"
                                    .format(
                                            surface_str[0],
                                            surface_str[1],
                                            surface_str[2]
                                           )
                              )
                
        # Write each alt line to the sounding file
        mod_sounding.write( 
                           "\n{:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}"
                                .format(
                                        alt_str[0],
                                        alt_str[1],
                                        alt_str[2],
                                        alt_str[3],
                                        alt_str[4] 
                                       )
                           )
        
    
#-----------------------------------------------------------------------------    
# Create estimated data between 15-20 km     
#-----------------------------------------------------------------------------
    
    # Boolean logic to add artifical data for 20 km depth
    increase_depth = True
    
    # Run this code if increase depth is turned on
    if( increase_depth == True ):
        
        # Create an array containing every 100 m interval between 15-20 km
        add_alt = np.arange( 15100.0000, 20100.0000, 100.0000 )
        
        # Pull out the 15 km theta value
        add_th = th[i]
        
        # Loop through the additional UL depth
        for k in range( 0, len( add_alt )  ):
            
            # Increase theta by 10 K per km ( Writing data for every 100 m so th + 1 )
            add_th = add_th + 1.0
            
            # Place alt data for current line into a list
            alt_str = [ 
                       add_alt[k],                       # Altitude (m)
                       add_th,                           # Theta (K) 
                       qv[ len(qv) - 1 ],                            # Mixing-Ratio (g/kg) (Assumed Constant)
                        u[ len(u) - 1 ],                             # U-comp (m/s) (Assumed Constant)
                        v[ len(v) - 1 ]                              # V-comp (m/s) (Assumed Constant)
                      ]
            
            # Construct and print current line using format descriptors
            print( 
                  " {:>12.4f}\t{:>12.4f}\t{:>12.4f}\t{:>12.4f}\t{:>12.4f}"
                      .format(
                              alt_str[0],
                              alt_str[1],
                              alt_str[2],
                              alt_str[3],
                              alt_str[4] 
                              )
                  )
                
            mod_sounding.write( 
                                "\n{:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}\t{:>9.4f}"
                                    .format(
                                            alt_str[0],
                                            alt_str[1],
                                            alt_str[2],
                                            alt_str[3],
                                            alt_str[4] 
                                           )
                               )

    # Report program status to terminal
    print( "\n{} has been succesfully modified to 20 km depth! ".format( files[j]) )
    
    # Report program status to terminal
    print( "\nNow closing {}\n\n".format( files[j] ) )
    
    # Close newly written sounding file
    mod_sounding.close()
    
    
    # Logic to determine if need to change working directory
    #-------------------------------------------------------
    
    # Change the working directory back to original location if additional soundings need to be modded
    if( j < len(files) - 1 ):
    
        # Return to Eastin-Original directory
        os.chdir( direct )
            
        # Retrieve the current directory
        cwd = os.getcwd()
        
        # Report working directory to user
        print('\nChanging Working Directory to: \n {}\n'.format( cwd ) )
    
    # Notify user that all modifications are complete    
    else:
        
        # Report program status to terminal
        print( "\nAll soundings have been successfully modified!" )
        
    # End change directory logic
    #--------------------------------------------------------
        
      
# End Sounding Loop
#-----------------------       
   

#-----------------------------------------------------------------------------
# End of main program
#-----------------------------------------------------------------------------

# Report the time required to run the full script
print( "\nScript Runtime: {}".format( datetime.now() - startTime ) )

# End of program flag
print( '\nScript Complete... Goodbye!' )   


