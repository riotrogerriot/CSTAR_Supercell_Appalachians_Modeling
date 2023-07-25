#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#-----------------------------------------------------------------------------
# Script: generate_realistic_terrain_from_dem_CM1.py
#
#   Update Records:
#       (4/12/2022) Script Created & Adpoted from Branden Katona's Jupyter Notebook
#                   (Realistic CM1 Terrain FINAL.ipynb)
#
#-----------------------------------------------------------------------------
#
# Summary:  
#   Takes USGS DEM data and turns it into a useable CM1 format with appropriate
#   filtering to reduce influence of small-scale gravity waves
#   Methodology modeled from Solderholm et al., 2014; Katona and Markowski 2021
# 
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: Load external libraries
#-----------------------------------------------------------------------------
from osgeo import gdal 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pyproj import Transformer
import scipy.interpolate as interp
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as patches
from datetime import datetime
import cartopy.feature as cfeat

#-----------------------------------------------------------------------------
# End: Load external libraries
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: User Defined Variables
#-----------------------------------------------------------------------------

# Program name
program = 'generate_realistic_terrain_from_dem_CM1.py'

# Define DEM latitudes (Must be integer in descending order)
lat_start = 40          # units: degrees N
lat_end = 34            # units: degrees N
file_lats = np.arange( lat_start, lat_end -1, -1 )

# Define DEM longitudes (Used negatives as standard convention for Western Hemisphere)
lon_start = -87         # units: degrees W
lon_end = -78           # units: degrees W
file_lons = np.arange( lon_start, lon_end + 1, 1 )

# # Desired UTM Zone to Convert from Lat/Lon to XY (meters) before interpolation
zone = 17

# Average sea-level height around terrain to reduce zs field by
reduc = 300.0

# Directory containing USGS DEM data
dir_to_data = '/Users/roger/Desktop/cm1_dem_data/'
# dir_to_data = '/Users/roger/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatCharlotte/CSTAR/Simulations/Realistic_Terrain/cm1_dem_data/'

# Requested horizontal gridspace resolution for CM1 interpolation (in meters)
cm1_dx = 1000.0 

#-----------------------------------------------------------------------------
# End: User Defined Variables
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Begin: Main Program
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( "\nBegin {}...".format( program ) )


#-----------------------------------------------------------------------------
# Begin: Read in USGS DEM Files (Assumed Northern/Western Hemisphere)
#-----------------------------------------------------------------------------

# Report status to terminal
print( 'Begin reading in DEM data...' )

# Initialize lat-counter variable
count_i = 0

#-----------------------------------------------------------------------------
# Begin: Loop through each lat in the provided USGS dataset
#-----------------------------------------------------------------------------
for i in file_lats:

    # Initialize lon-counter variable (Needs to restart for each new latitude!)
    count_j = 0

    #-----------------------------------------------------------------------------
    # Begin: Loop through each lon in the provided USGS dataset
    #-----------------------------------------------------------------------------    
    for j in file_lons:
        
        # Construct current file string 
        # (Note: look into what makes this string constructor different from .format() method! )
        file_string = dir_to_data + ( 'USGS_1_n%02dw%03d.tif' %( i, j * -1 ) )  
        
        # Report to terminal for verification
        # print( '\tCurrent DEM: {}'.format( file_string ) )
        
        # Open the current USGS DEM file using GDAL
        geo = gdal.Open( file_string )
        
        # Read in height daa from the current DEM file
        dem_heights_temp = geo.ReadAsArray()
        
        #Grab some info about the XY dimensions of the current DEM
        width = geo.RasterXSize
        height = geo.RasterYSize
        
        # Geo-reference the DEM image via Affine method
        # (See https://gdal.org/tutorials/geotransforms_tut.html for additional info)
        gt = geo.GetGeoTransform()
        
        # Use geo-referenced metadata to construct geographic BB_Box coordinates
        minx = gt[0]
        miny = gt[3] + width * gt[4] + height * gt[5] 
        maxx = gt[0] + width * gt[1] + height * gt[2]
        maxy = gt[3] 
        
        # Create arrays containing lat/lon of current DEM via BB_Box and Grid Spacing
        # which was collected during Geo-Referencing
        lon_tmp = np.arange( minx, maxx, gt[1] )
        lat_tmp = np.arange( miny, maxy, gt[5] * -1 )
        
        
        # Check to assure Geo-Referenced Longitude Data is same size as Original data
        # (Rounding errors could potential change the array size so we may need to 
        # trim to match the original shape)
        #-----------------------------------------------------------------------------
        if np.size( lon_tmp ) > width:
         
            # Geo-Referenced Lon Data > OG data so trim it down
            size_diff = np.size( lon_tmp ) - width
            trim = -1 * ( size_diff )
            lon_tmp = lon_tmp[ 0:trim ]
        
        # During the first iteration construct an array to save all longitude height data too
        if count_j == 0:
            dem_heights_lons = dem_heights_temp
            lons = lon_tmp
            
        # During subsequent iterations append the current longitude height data to the final array
        if count_j > 0:
            
            #Append the current terrain field to the old one
            dem_heights_lons = np.append( dem_heights_lons, dem_heights_temp, axis = 1 )
            lons = np.append( lons, lon_tmp )
        
        # Add to the loop counters now that we are done with the current iteration
        count_j += 1
              
    #-----------------------------------------------------------------------------
    # End: Loop through each lon in the provided USGS dataset
    #-----------------------------------------------------------------------------  
        
        
    # Check to assure Geo-Referenced Latitude Data is same size as Original data
    # (Rounding errors could potential change the array size so we may need to 
    # trim to match the original shape)
    #-----------------------------------------------------------------------------   
    if np.size(lat_tmp) > height:
        
        # Geo-Referenced Lat Data > OG data so trim it down
        size_diff = np.size( lat_tmp ) - height
        trim = -1 * ( size_diff + 1 )
        lat_tmp = lat_tmp[ 0:trim ]
    
    #-----------------------------------------------------------------------------
    # Begin: Latitude Tiled DEM Array construction
    #-----------------------------------------------------------------------------
        
    # During the first iteration construct an array to save all latitude height data too
    if count_i == 0:
        dem_heights_lats = dem_heights_lons
        lats = lat_tmp
        
    # During subsequent iterations append the current latitude height data to the final array
    if count_i > 0:
        
        #Again, append the current terrain field to the old one
        old_lats = lats
        dem_heights_lats = np.append( dem_heights_lats, dem_heights_lons, axis = 0 )
        old_lats = lats
        lats = np.append( lat_tmp, old_lats )
        
    # Add to the loop counters now that we are done with the current iteration
    count_i += 1

    #-----------------------------------------------------------------------------
    # End: Latitude Tiled DEM Array construction
    #-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# End: Loop through each lat in the provided USGS dataset
#-----------------------------------------------------------------------------

# Report status to terminal
print( 'Read in DEM data complete!' )


#-----------------------------------------------------------------------------
# Begin: Reproject data from Lat/Lon to UTM
#-----------------------------------------------------------------------------

# Report status to terminal
print( 'Begin conversion from Lat/Lon to UTM...' )

# Now we need to flip the data because we appened from N -> S & need structure from S -> N
# (See filename printout to verify this!)
terrain_height_dem = dem_heights_lats
lats_dem = lats
lons_dem = lons

# Get Transformation info for Lat/Lon to requested UTM Zone 
trans = Transformer.from_crs( 
                              "epsg:4326",
                              "proj=utm +zone={} +ellps=WGS84".format( str(zone) ),
                              always_xy = True
                            )



# Use meshgrid to force Lat/Lon data to be placed on equal sized grid
lons_transf, lats_transf = np.meshgrid( lons_dem, lats_dem )
print( 'Start trans...' )
# Perform the reprojection from Lat/Lon to UTM
xx, yy = trans.transform( lons_transf, lats_transf )
print( 'End trans...' )
#Set UTM coordinates with 0,0 at the lower left corner of the domain
xx_0 = xx[ 0, : ] - xx[ 0, 0 ]
yy_0 = yy[ :, 0 ] - yy[ 0, 0 ]

# Report status to terminal
print( 'End conversion from Lat/Lon to UTM!' )

#-----------------------------------------------------------------------------
# End: Reproject data from Lat/Lon to UTM
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: Interpolate Tiled Heights to CM1 Grid
#-----------------------------------------------------------------------------

# Report status to terminal
print( 'Begin Bilinear Interpolation to CM1 Grid...' )

print( 'Creating CM1 Grid...' )

#Set up a grid with CM1 dx/dy. This will set up the interpolation to these values
interp_x_vals = np.arange( 0, xx_0[-1], cm1_dx )
interp_y_vals = np.arange( 0, yy_0[-1], cm1_dx )

XX,YY = np.meshgrid( xx_0, yy_0 )

print ( 'Performing Interpolation...' )
 
f = interp.interp2d( xx_0, yy_0, terrain_height_dem[ 0:-5 ] )
# f = interp.RectBivariateSpline( xx_0, yy_0, terrain_height_dem[0:-5] )

print( 'Interpolating terrain to CM1 grid...' )

# Perform the actual bilinear interpolation at the CM1 coordinates
interp_zs = f( interp_x_vals, interp_y_vals )

# Report status to terminal
print( 'End Bilinear Interpolation to CM1 Grid!' )

#-----------------------------------------------------------------------------
# End: Interpolate Tiled Heights to CM1 Grid
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: 6dx Gaussian Filtering (See Solderholm et al., 2014 for more details)
#-----------------------------------------------------------------------------

# Create meshgrids for plotting for both the interp vals and the original grid
X,Y = np.meshgrid( interp_x_vals, interp_y_vals )
XX,YY = np.meshgrid( xx_0, yy_0 )

# Rearrange the array so it's right for CM1
zs_tmp = np.flipud( interp_zs )

# This gets the "standard deviation" for the gaussian filter to filter out waves 6*dx and lower
# This is taken from a MATLAB program I downloaded called FILT2D which is a geospatial filter
sigma_guess = 6.0 / ( 2 * np.pi )

# Perform the low pass filtering
zs_filt = gaussian_filter( zs_tmp, sigma = ( sigma_guess, sigma_guess ) )

#-----------------------------------------------------------------------------
# End: 6dx Gaussian Filtering (See Solderholm et al., 2014 for more details)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: Relief Preservation w/ rpt Sea-Level
#-----------------------------------------------------------------------------

# Lower the terrain height by 300 m  
zs_lower = zs_filt - reduc

# If the terrain is lower than 0, set it to 0
zs_lower = np.where( zs_lower < 0, 0, zs_lower )

#-----------------------------------------------------------------------------
# End: Relief Preservation w/ rpt Sea-Level
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: Linear Decay outward from Ellipse
#-----------------------------------------------------------------------------

# This is where we define the parameters of the ellipse. The values here are in meters.
# I just eyeballed the values until they looked correct and did some trial and error.
# If you want to change these, just make some plots with the ellipse plotted on top and you can
# sort of adjust it as needed
g_ell_center = ( 425000, 320000 )        # Center of elipse
g_ell_width = 250000                     # Width of ellipse
g_ell_height = 125000                    # height of ellipse
angle = 35.0                             # Angle of ellipse
cos_angle = np.cos( np.radians( 180.0 - angle ) )
sin_angle = np.sin( np.radians( 180.0 - angle ) )

xc = X - g_ell_center[0] # Center coordinates on the ellipse
yc = Y - g_ell_center[1] 

xct = xc * cos_angle - yc * sin_angle
yct = xc * sin_angle + yc * cos_angle 

#Calculates the radius of the ellipse from the center points
rad_cc = ( xct**2 / ( g_ell_width/2.0)**2 ) + ( yct**2 / ( g_ell_height/2.0 )**2 )

# Keep the terrain within the calculated ellipse (where radius <= 1)
terr_ellipse = np.where( rad_cc <= 1.0, zs_lower, 0 )

#Between 1 and 3 ellipse radii, linearly smooth the terrain down to 0
terr_ellipse = np.where( ( rad_cc >= 1 ) & ( rad_cc <= 3 ), zs_lower * ( (3 - rad_cc ) / ( (3 - 1) ) ), terr_ellipse )

#-----------------------------------------------------------------------------
# End: Linear Decay outward from Ellipse
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: Figure Generation
#-----------------------------------------------------------------------------

z_min = 0
z_max = 2050
z_int = 50 

# Adjust these numbers to modify domain box on plots
cm1_x1 = 50000
cm1_y1 = 75000
cm1_xlen = 600000               # Units: meter
cm1_ylen = 450000               # Units: meter

# Create a 4-panel mosaic
figure_mosaic = """
                 AABB
                 CCDD
                """ 

# Create figure and axes objects
fig, axes = plt.subplot_mosaic( 
                               mosaic = figure_mosaic, 
                               figsize = ( 12, 9.5 ), 
                               tight_layout = True
                              )


# Panel A: Actual Topography of Southern Appalachians
#----------------------------------------------------

# Create the panel A axis object
ax4 = fig.add_subplot( axes['A'] )

# Add panel annotation to upper-left corner
ax4.set_title( 'a: DEM Topography Scaled to CM1 Grid', loc = 'left', fontsize = 16, fontweight = 'bold' )

# Use a filled contour to map the actual terrain to panel A
cax4 = ax4.contourf(
                    X,Y, 
                    zs_tmp,
                    levels = np.arange( z_min, z_max, z_int ),
                    cmap = mpl.cm.copper.reversed(),
                    alpha = 0.30
                   )

# Plot a rectangle representative of the CM1 Model Domain
rect = patches.Rectangle(
                         (cm1_x1, cm1_y1),
                         cm1_xlen, cm1_ylen,
                         fill = False,
                         edgecolor = 'black',
                         linestyle = '--',
                         linewidth = 3,
                         label = 'Model Domain'
                        )

# Add model domain rectangle to panel A
ax4.add_patch(rect)

# # Plot an ellipse around terrain feature of interest
# g_ellipse = patches.Ellipse(
#                             g_ell_center,
#                             g_ell_width,
#                             g_ell_height,
#                             angle = angle,
#                             fill = False,
#                             edgecolor = 'black',
#                             linewidth = 3,
#                             label = 'AOI' 
#                            )

# # Add the ellipse to the plot
# ax4.add_patch(g_ellipse)

# ax4.set_yticklabels( labels = np.arange( 34, 40, 0.75 ), fontsize = 14, fontweight = 'bold' )
# ax4.set_xticklabels( np.arange( 87, 78, -1.5 ), fontsize = 14, fontweight = 'bold' )

# ax4.set_ylabel( 'Latitude (N)', fontsize = 14, fontweight = 'bold' )
# ax4.set_xlabel( 'Longitude (W)', fontsize = 14, fontweight = 'bold' )


ax4.grid()
ax4.axes.xaxis.set_visible( False )
ax4.axes.yaxis.set_visible( False )


# Panel B: 6dx Filtering
#----------------------------------------------------

# Create the panel A axis object
ax3 = fig.add_subplot( axes['B'] )

# Add panel annotation to upper-left corner
ax3.set_title( 'b: 6dx Filtering & Reduction to Sea-Level', loc = 'left', fontsize = 16, fontweight = 'bold' )

# Use a filled contour to map the actual terrain to panel A
cax3 = ax3.contourf(
                    X,Y, 
                    zs_lower,
                    levels = np.arange( z_min + 0.01, z_max, z_int ),
                    cmap = mpl.cm.copper.reversed(),
                    alpha = 0.30
                   )

# Plot a rectangle representative of the CM1 Model Domain
rect = patches.Rectangle(
                         (cm1_x1, cm1_y1),
                         cm1_xlen, cm1_ylen,
                         fill = False,
                         edgecolor = 'black',
                         linestyle = '--',
                         linewidth = 3,
                         label = 'Model Domain'
                        )

# Add model domain rectangle to panel A
ax3.add_patch(rect)

# Plot an ellipse around terrain feature of interest
g_ellipse = patches.Ellipse(
                            g_ell_center,
                            g_ell_width,
                            g_ell_height,
                            angle = angle,
                            fill = False,
                            edgecolor = 'black',
                            linewidth = 3,
                            label = 'Area of Interest' 
                           )

# Add the ellipse to the plot
ax3.add_patch(g_ellipse)

ax3.grid()
ax3.axes.xaxis.set_visible( False )
ax3.axes.yaxis.set_visible( False )


# Panel C: Linear Decay outward from Ellipse
#----------------------------------------------------

# Plot to show the terrain we are keeping for the CM1 terrain file
ax2 = fig.add_subplot( axes['C'] )

# Add panel annotation to upper-left corner
ax2.set_title( 'c: Outward Linear Decay', loc = 'left', fontsize = 16, fontweight = 'bold' )

# Use a filled contour to map the actual terrain to panel A
cax2 = ax2.contourf(
                    X,Y, 
                    terr_ellipse,
                    levels = np.arange( z_min + 0.01, z_max, z_int ),
                    cmap = mpl.cm.copper.reversed(),
                    alpha = 0.30
                   )

# Plot a rectangle representative of the CM1 Model Domain
rect = patches.Rectangle(
                         (cm1_x1, cm1_y1),
                         cm1_xlen, cm1_ylen,
                         fill = False,
                         edgecolor = 'black',
                         linestyle = '--',
                         linewidth = 3,
                         label = 'Model Domain'
                        )

# Add model domain rectangle to panel A
ax2.add_patch(rect)

# Plot an ellipse around terrain feature of interest
g_ellipse = patches.Ellipse(
                            g_ell_center,
                            g_ell_width,
                            g_ell_height,
                            angle = angle,
                            fill = False,
                            edgecolor = 'black',
                            linewidth = 3,
                            label = 'Area of Interest' 
                           )

# Add the ellipse to the plot
ax2.add_patch(g_ellipse)

ax2.grid()
ax2.axes.xaxis.set_visible( False )
ax2.axes.yaxis.set_visible( False )

# Add a legend representative for panels A-C
ax2.legend()

# Panel D: CM1-Ready Terrain
#----------------------------------------------------

#This defines the bounding box (put in meters from the original large domain) 
# Adjust these numbers to add to model domain
xbox_0 = 25000
xbox_1 = 625000
ybox_0 = 125000
ybox_1 = 625000

#Grab the coordinates of the bounding box
x_start = np.where(interp_x_vals==xbox_0)[0][0]
x_end = np.where(interp_x_vals==xbox_1)[0][0]
y_start = np.where(interp_y_vals==ybox_0)[0][0]
y_end = np.where(interp_y_vals==ybox_1)[0][0]
zs_save = terr_ellipse[y_start:y_end,x_start:x_end]

#Get the shape
zs_shape = np.shape(zs_save)

#Put the zs field into a writeout variable (Can't remember why I did this, but it's here. Might be unnecessary)
zs_writeout = zs_save

ax1 = fig.add_subplot( axes['D'] )

# Add panel annotation to upper-left corner
ax1.set_title( 'd: CM1 Realistic Terrain Field', loc = 'left', fontsize = 16, fontweight = 'bold' )

# Use a filled contour to map the actual terrain to panel A
cax1 = ax1.contourf(
                    zs_writeout,
                    levels = np.arange( z_min + 0.01, z_max, z_int ),
                    cmap = mpl.cm.copper.reversed(),
                    alpha = 0.30
                   )

# Add grid
ax1.grid()

# ax1.set_xticklabels( labels = np.arange( 0, 600, 120 ), fontsize = 14, fontweight = 'bold' )
# ax1.set_yticklabels( np.arange( 0, 500, 100 ), fontsize = 14, fontweight = 'bold' )

ax1.set_xlabel( 'Zonal Distance (km)', fontsize = 14, fontweight = 'bold' )
ax1.set_ylabel( 'Meridional Distance (km)', fontsize = 14, fontweight = 'bold' )

# Create universal colorbar axes
cb_ax = fig.add_axes( [1.01, 0.1, 0.02, 0.8] )

# Use panel A data to make universal colorbar
cbar = fig.colorbar( cax4, cax = cb_ax )

# Set colorbar label 
cbar.set_label( label = 'Elevation (m)', fontsize = 16, fontweight = 'bold' )

# Set colorbar axes labels
cb_ax.tick_params( axis = 'y', labelsize = 14 )

# Save the current plot
fig.savefig(
            fname = 'cm1_terrain_fig.jpeg',
            dpi = 300,
            bbox_inches = "tight"
           )

#-----------------------------------------------------------------------------
# End: Figure Generation
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: perts.dat Writeout (File to provide to CM1)
#-----------------------------------------------------------------------------
#Open the binary file we'll be writing to
file = open( 'perts_1000m.dat', 'wb' )

#Convert to float 32 because that's what CM1 wants (discovered by trial and error and Stack Overflow)
terrain_final = np.float32( zs_writeout )

#Get shape of this field
terr_shape = np.shape( terrain_final )

#Write the terrain field line by line to the file (this is fast)
for i in range( terr_shape[0] ):
    lineArr = terrain_final[i,:]
    lineArr.tofile(file)

#Close the file 
file.close()

#-----------------------------------------------------------------------------
# End: perts.dat Writeout (File to provide to CM1)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# End: Main Program
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# Confirm that script successfully ran
print( "\n{} successfully completed!".format( program ) )

# Report the time required to run the function
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )






