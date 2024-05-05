#-----------------------------------------------------------------------------
# Function: terrain_induced_pert_from_BS.py
#
#   Update Records:
#       (10/9/2023) Script Created
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
# Begin: Load external libraries
#-----------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.patches as patches
import numpy as np
import scipy
import os
import xarray as xr
import metpy
from metpy.units import units
import metpy.calc as mpcalc
from datetime import datetime
from srh_cy import srh
import pint
import pint_xarray
import cf_xarray
#-----------------------------------------------------------------------------
# End: Load external libraries
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: User-Defined Variables
#-----------------------------------------------------------------------------

# Script Name
program = 'terrain-induced_Perts_from_BS_w_FrH_test.py'

# Record Script start time
startTime = datetime.now()

# Report program status to terminal
print( "\nBegin {}...".format( program ) )

# Path to model data
path = '/scratch/rriggin/cm1r20.3_init3d_mod/dry_simulations/'

# Grid spacing (m) and initial sounding location
dx = 250.0 
dy = dx
sx = int( 50000.0 / dx )
sy = int( 50000.0 / dy )

# Common Height Indices (Must change if z-stretch setting are modified)
m500 = 4
km1 = 6
km3 = 13
km5 = 18
km8 = 24
km10 = 28

# Manually set wind vector intervals
skip_val = 75

# Wind vector scale
v_scale  = 0.60

# Peak position of IDTRN (km)
peak_pos = 350.0 

# Array of simulation names for looping purposes (0 == IDTRN 1 = MOD_IDTRN)
sims_bool = 0

if( sims_bool == 0 ):
    sims = [ 'cs_idtrn', 'nc_idtrn' ]
    terrain = 'IDTRN'
    
    # Output file parameters (for upcoming output file loop)
    i_start = 1                                     # Start Time
    dt = 5.0                                        # units: (min)
    t_total = 4.0                                   # units: (hrs)
    nt = int( ( t_total * 60.0 ) / dt )             # Total simulation timesteps
    t_int = int(nt/2)                               # File Looping interval
    i_end = nt + 1                                  # File No. to end loop

elif( sims_bool == 1 ):
    sims = [ 'cs_mod', 'nc_mod' ]
    terrain = 'MOD_IDTRN'
    
    # Output file parameters (for upcoming output file loop)
    i_start = 1                                     # Start Time
    dt = 5.0                                        # units: (min)
    t_total = 4.0                                   # units: (hrs)
    nt = int( ( t_total * 60.0 ) / dt )             # Total simulation timesteps
    t_int = int(nt/2)                               # File Looping interval
    i_end = nt + 1                                  # File No. to end loop

else:
    sims = [ 'cs_rltrn', 'nc_rltrn' ]
    terrain = 'RLTRN'
    
    # Output file parameters (for upcoming output file loop)
    i_start = 25                                    # Start at 2 hr mark to allow for terrain adjusments
    dt = 5.0                                        # units: (min)
    t_total = 6.0                                   # units: (hrs)
    nt = int( ( t_total * 60.0 ) / dt )             # Total simulation timesteps
    t_int = 48                                      # File Looping interval
    i_end = nt + 1                                  # File No. to end loop

# Min/Max value ranges and increments for base-state var plots
### Note: CIN values need to be reversed in order to plot in line with negative CIN -> Bad for Storms
if( terrain != 'RLTRN' ):
    cape_min = 500
    cin_min = -25
    srh_min = 50
    FrH_min = 0.0
    
    cape_max = 1500
    cin_max = 0
    srh_max = 150
    FrH_max = 1.5
    
    cape_int = 50
    cin_int = 2.5
    srh_int = 10
    Frh_int = 0.1
    
    # Min/Max value ranges and increments for pert var plots
    capep_min = -500
    cinp_min = -50
    srhp_min = -100
    FrHp_min = -0.5
    
    capep_max = 500
    cinp_max = 50
    srhp_max = 100
    FrHp_max = 0.5
    
    capep_int = 25
    cinp_int = 5
    srhp_int = 5
    FrHp_int = 0.05
    
else:
    cape_min = 500
    cin_min = -25
    srh_min = 50
    FrH_min = 0.25
    
    cape_max = 2000
    cin_max = 0
    srh_max = 200
    FrH_max = 1.5
    
    cape_int = 100
    cin_int = 2.5
    srh_int = 10
    Frh_int = 0.1
    
    # Min/Max value ranges and increments for pert var plots
    capep_min = -500
    cinp_min = -50
    srhp_min = -100
    FrHp_min = -0.5
    
    capep_max = 500
    cinp_max = 50
    srhp_max = 100
    FrHp_max = 0.5
    
    capep_int = 25
    cinp_int = 5
    srhp_int = 5
    FrHp_int = 0.05    

#-----------------------------------------------------------------------------
# End: User-Defined Variables
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Set-Up Working and Plotting Directories
#-----------------------------------------------------------------------------

# Change the working directory to the provided path
os.chdir( path )

# Report to terminal
print( "\nChanging working directory to model output location... \n\t{}".format( path ) )

plot_dir = path + "/plots/"

# Create a new directory for plots if it does not exist
#------------------------------------------------------
if( os.path.isdir( plot_dir ) == False ):
    
    # Report to terminal1
    print( "\nCreating new directory to store plots...\n\t{}".format( plot_dir ) )
    
    # Create directory
    os.mkdir( plot_dir )
#------------------------------------------------------

#-----------------------------------------------------------------------------
# End Set-Up Working and Plotting Directories
#-----------------------------------------------------------------------------


# Create figure
#-----------------------------------------------------------------------------

# Representative of the individual plots that compose the figure 
figure_mosaic = """
                ADGJM
                BEHKN
                CFILO
                PQRST
                """ 

# Set tick size (Must be before calling plot object)
plt.rcdefaults()
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

# Create figure and axes objects
fig, axes = plt.subplot_mosaic( mosaic = figure_mosaic, figsize = ( 40, 28 ), tight_layout = True, )
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Begin: Simulation loop
#-----------------------------------------------------------------------------
for j in range( 0, len(sims) ):
    
    # Change the working directory to appropriate simulation directory
    #-----------------------------------------------------------------------------
    cwd = path + "{}/".format( sims[j] )
    os.chdir( cwd )
    
    # Report terrain type and CWD to terminal
    print( '\n\tTerrain Type: {} \tCWD: {}'.format( terrain, cwd ) )
    #-----------------------------------------------------------------------------
    
    
    # Define Model Domain and plotting extents based on terrain type 
    #-----------------------------------------------------------------------------
    if( terrain == 'IDTRN' or terrain == 'MOD_IDTRN' ):
        
        # Domain (Num of gridpoints)
        nx = 2400
        ny = 1600
        nz = 48
        
        # Plotting Extent (km)
        x1 = -200
        x2 = 200
        y1 = 50
        y2 = 250
        
        # Windward averaging locations
        wwX1 = 700
        wwX2 = 1100
        wwY1 = 300
        wwY2 = 900
        
        # Leeward averaging locations
        lwX1 = 1700
        lwX2 = 2100
        lwY1 = 300
        lwY2 = 900
        
        # Peak Elevation & Terrain Contour Interval (m)
        inc = 150.0
        peak_elv = 750.0 
    
    # RLTRN
    else:
        
        # Domain (Num of gridpoints)
        nx = 2400
        ny = 2000
        nz = 48
        
        # Plotting Extent (km)
        x1 = 175
        x2 = 525
        y1 = 150
        y2 = 350
        
        # # Windward averaging locations
        wwX1 = 800
        wwX2 = 1200
        wwY1 = 700
        wwY2 = 1100
        
        # # Leeward averaging locations
        lwX1 = 1600
        lwX2 = 2000
        lwY1 = wwY1 + 100
        lwY2 = wwY2 + 100
        
        # Peak Elevation & Terrain Contour Interval (m)
        inc = 500.0
        peak_elv = 1500.0 

    #-----------------------------------------------------------------------------  
        
    
    # Initialize base-state arrays
    #-----------------------------------------------------------------------------
    print( '\n\tInitializing Base-state arrays...' )
    uinterp_bs = np.zeros( shape = ( km1, ny, nx ) )
    vinterp_bs = np.zeros( shape = ( km1, ny, nx ) )
    cape_bs = np.zeros( shape = ( ny, nx ) )
    cin_bs = np.zeros( shape = ( ny, nx ) )
    srh1km_bs = np.zeros( shape = ( ny, nx ) )
    FrH_bs = np.zeros( shape = ( ny, nx ) )


    # Initialize perturbation from base-state arrays
    #-----------------------------------------------------------------------------
    print( '\n\tInitializing perturbation arrays...' )
    uinterp_pert = np.zeros( shape = ( nt, km1, ny, nx ) )
    vinterp_pert = np.zeros( shape = ( nt, km1, ny, nx ) )
    cape_pert = np.zeros( shape = ( nt, ny, nx ) )
    cin_pert = np.zeros( shape = ( nt, ny, nx ) )
    srh1km_pert = np.zeros( shape = ( nt, ny, nx ) )
    FrH_pert = np.zeros( shape = ( nt, ny, nx ) )
    
    # Initialize windward and leeward averages
    #-----------------------------------------------------------------------------
    ww_cape_bs_avg = 0
    ww_cin_bs_avg = 0
    ww_srh_bs_avg = 0
    ww_FrH_bs_avg = 0
    lw_cape_bs_avg = 0
    lw_cin_bs_avg = 0
    lw_srh_bs_avg = 0
    lw_FrH_bs_avg = 0
    ww_cape_pert_avg = 0
    ww_cin_pert_avg = 0
    ww_srh_pert_avg = 0
    ww_FrH_pert_avg = 0
    lw_cape_pert_avg = 0
    lw_cin_pert_avg = 0
    lw_srh_pert_avg = 0
    lw_FrH_pert_avg = 0
    
    #-----------------------------------------------------------------------------  
    
    
    #-----------------------------------------------------------------------------
    # Begin: Output loop
    #-----------------------------------------------------------------------------
    for i in range( i_start, i_end + 1, t_int ):

        #-----------------------------------------------------------------------------
        # Begin: Data Collection and Manipulation
        #-----------------------------------------------------------------------------
        
        print( '\n\t\tOutput Loop: i: {} Time: {} min...'.format( i, str((i-1) * dt) ) )
        
        # Simulation time logic
        if( i == i_start):
            t0 = (i -1) * dt
        else:
            mtime = (i-1) * dt
            
        # Logic to open terrain-interpolated netCDF
        #-----------------------------------------------------------------------------
       
        # Define the non-interpolated filename as a string
        if ( i < 10 ):
            filename2 = 'cm1out_00000' + str(i) + '.nc'
            
            # Report program status to terminal
            print( "\n\t\tOpening " + filename2 + "..." )
            
        elif( i < 100 and i >= 10 ):
            filename2 = 'cm1out_0000' + str(i) + '.nc'
            
            # Report program status to terminal
            print( "\n\t\tOpening " + filename2 + "..."  )
            
        else:
            filename2 = 'cm1out_000' + str(i) + '.nc'
            print( "\n\t\tOpening " + filename2 + "...")
    
        # Open the current terrain-interpolated netCDF file with xarray 
        DS = xr.open_dataset( filename2, engine = "netcdf4", decode_cf = True )
        #-----------------------------------------------------------------------------
        
                        
        # Define the terrain-interpolated filename as a string
        #-----------------------------------------------------------------------------
        if ( i < 10 ):
            filename = 'cm1out_00000' + str(i) + '_i.nc'
            
            # Report program status to terminal
            print( "\n\t\tOpening " + filename + "..." )
            
        elif( i < 100 and i >= 10 ):
            filename = 'cm1out_0000' + str(i) + '_i.nc'
            
            # Report program status to terminal
            print( "\n\t\tOpening " + filename + "..."  )
            
        else:
            filename = 'cm1out_000' + str(i) + '_i.nc'
            print( "\n\t\tOpening " + filename + "...")
    
        # Open the current terrain-interpolated netCDF file with xarray 
        ds = xr.open_dataset( filename, engine = "netcdf4", decode_cf = True )
        #-----------------------------------------------------------------------------
        
        
        # Read-in current output variables
        #-----------------------------------------------------------------------------
        print( '\n\t\t\tPull required vars from output file...' )

        # Get the u-component of the wind (m/s)    
        u = DS.metpy.parse_cf( 'uinterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
        u = u.metpy.quantify()
        
         # Get the v-component of the wind (m/s)   
        v = DS.metpy.parse_cf( 'vinterp' ).isel( time = 0, zh =  slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
        v = v.metpy.quantify()
        
        # Get pressure (Pa) and convert to hPa
        P = DS.metpy.parse_cf( 'prs' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
        P = (P/100.0) * units('hPa')
        
        # Get potential temperature perturbation (K) 
        th = DS.metpy.parse_cf( 'th' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
        th = th.metpy.quantify()
        
        # Get water-vapor mixing ratio (kg/kg) and convert to g/kg
        qv = DS.metpy.parse_cf( 'qv' ).isel( time = 0, zh = slice( 0, nz ), yh = slice( 0, ny ), xh = slice( 0, nx ) )
        qv = qv.metpy.quantify()
        qv = qv.pint.to( 'g/kg' )
        
        # Get CAPE (J/kg)
        cape = ds.metpy.parse_cf( 'cape' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
        
        # Get CIN (J/kg)
        cin = ds.metpy.parse_cf( 'cin' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
        cin = -1 * cin
        
        # Get terrain (m)
        zs = ds.metpy.parse_cf( 'zs' ).isel( time = 0, yh = slice( 0, ny ), xh = slice( 0, nx ) )
        zs = zs.values

        # Coordinate Variables (x,y,z in km)
        xh = ds.coords[ 'xh' ].isel( xh = slice( 0, nx ) ).values * units( 'km' )      
        yh = ds.coords[ 'yh' ].isel( yh = slice( 0, ny ) ).values * units( 'km' )     
        zh = ds.coords[ 'zh' ].isel( zh = slice( 0, nz ) ).values * units( 'km' )
        
        # Point variables from init sounding location
        u_sm = DS.metpy.parse_cf( 'uinterp' ).isel( time = 0, yh = sy, xh = sx )
        v_sm = DS.metpy.parse_cf( 'vinterp' ).isel( time = 0, yh = sy, xh = sx )
        pressure = DS.metpy.parse_cf( 'prs' ).isel( time = 0, yh = sy, xh = sx )
        #-----------------------------------------------------------------------------
        
        
        # Close the current file
        #-----------------------------------------------------------------------------
        print( "\n\t\tClosing " + filename + "...")
        print( "\n\t\tClosing " + filename2 + "...")
        DS.close()
        ds.close()
        #-----------------------------------------------------------------------------
        
        
        #-----------------------------------------------------------------------------
        # Begin: Calculations
        #-----------------------------------------------------------------------------
        
        # Compute Bunkers Storm Motions and Mean Wind for SRH Calculations
        # Assuming constant storm-motion due to 1D limitations of Metpy Bunkers Calculation
        # Would probably be better to let storm-motion vary with grid point due to flow changes
        # within mountainous areas...
        print( f'\n\t\t\tComputing storm-motion... {sims[j]}, {(i-1)*dt} min' )
        rm, lm, mean_uv = mpcalc.bunkers_storm_motion( pressure, u_sm, v_sm, zh )
        
        # Convert zh into 3D array required for SRH function
        print( f'\n\t\t\tConverting zh into 3D grid... {sims[j]}, {(i-1)*dt} min' )
        zh_grid = np.zeros( shape = (nz, ny, nx) )
        
        # Interpolate heights to a 3D domain required for SRH function
        #-------------------------------------------------------------
        for k in range( 0, nz ):
            zh_grid[k,:,:] = zh[k].m
            
            # Convert to meters
            zh_grid[k,:,:] = zh_grid[k,:,:] * 1000.0
        #-------------------------------------------------------------
        
        # Compute 0-1 km SRH using Cython function (Must convert from pint to np.array using .values method)
        print( f'\n\t\t\tComputing SRH... {sims[j]}, {(i-1)*dt} min' )
        srh1km = srh( u = u.values, v = v.values,                 # units: m/s
                      zh = zh_grid, zs = zs,                      # units: m
                      u_storm = rm[0].m, v_storm = rm[1].m,       # units: m/s
                      layer = 1000.0 )                            # units: m
        
        
        # Compute theta-v via metpy
        temp = mpcalc.temperature_from_potential_temperature( pressure, th ) 
        temp = temp.metpy.quantify()
        temp = temp.pint.to( 'degC' )
        th_v = mpcalc.virtual_potential_temperature( pressure, temp, qv )
        
        
        # Get 0-1 km moist Brunt-Vaisala Frequency using Th-V across grid and average the results over 0-1 km z-dim
        print( f'\n\t\t\tComputing Brunt-Vaisala Frequency... {sims[j]}, {(i-1)*dt} min' )
        N = mpcalc.brunt_vaisala_frequency( zh_grid[0:km1,:,:]*units('m'), th_v[0:km1,:,:], vertical_dim = 0 )
        N2 = mpcalc.brunt_vaisala_frequency_squared( zh_grid[0:km1,:,:]*units('m'), th_v[0:km1,:,:], vertical_dim = 0 )
        ### Note: Have to go up one grid from the surface to prevent nans with NC environment
        N_avg = np.mean( N[1:], axis = 0 )
        
        # Calculate the mountain froude number (FrH = U_avg/(N_avg*H) ) [Markowski and Richarson 2010]
        print( f'\n\t\t\tComputing Mountain Froude Number... {sims[j]}, {(i-1)*dt} min' )
        u_avg = np.mean( u[1:km1,:,:], axis = 0 ) 
        FrH = ( u_avg ) / ( N_avg * ( (np.amax(zs) + 1) * units( 'm' ) ) )

    
        #-----------------------------------------------------------------------------
        # End: Calculations
        #-----------------------------------------------------------------------------
    
        
        # Collect Base-State Variables
        #-----------------------------------------------------------------------------
        if( i == i_start ):
            print( f'\n\t\tStoring base-state values... {sims[j]}, {(i-1)*dt} min' )
            uinterp_bs = np.copy(u[0:km1,:,:])                   # units: (m/s)
            vinterp_bs = np.copy(v[0:km1,:,:])                   # units: (m/s)
            cape_bs = np.copy(cape.values)                       # units: (J/kg)
            cin_bs = np.copy(cin.values)                         # units: (J/kg)
            srh1km_bs = np.copy(srh1km.base)                     # units: (m2/s2)
            FrH_bs = np.copy(FrH.values)                         # units: (n/a)
            
            # Calculate windward averages 
            #-----------------------------------------------------------------------------
            print( f'\n\tCalculate Windward Base-State Averages... {sims[j]}, {(i-1)*dt} min' )
            ww_cape_bs_avg = np.mean( cape_bs[wwY1:wwY2, wwX1:wwX2] )
            ww_cin_bs_avg = np.mean( cin_bs[wwY1:wwY2, wwX1:wwX2] )
            ww_srh_bs_avg = np.mean( srh1km_bs[wwY1:wwY2, wwX1:wwX2] )
            ww_FrH_bs_avg = np.nanmean( FrH_bs[wwY1:wwY2, wwX1:wwX2] )
            print( f'\t\tWindward Average: Frh BS = {ww_FrH_bs_avg}' )
            
            # Calculate leeward averages 
            #-----------------------------------------------------------------------------
            print( f'\n\tCalculate Leeward Base-State Averages... {sims[j]}, {i*dt} min' )
            lw_cape_bs_avg = np.mean( cape_bs[lwY1:lwY2, lwX1:lwX2] )
            lw_cin_bs_avg = np.mean( cin_bs[lwY1:lwY2, lwX1:lwX2] )
            lw_srh_bs_avg = np.mean( srh1km_bs[lwY1:lwY2, lwX1:lwX2] )
            lw_FrH_bs_avg = np.nanmean( FrH_bs[lwY1:lwY2, lwX1:lwX2] )
            print( f'\t\tLeeward Average: Frh BS = {lw_FrH_bs_avg}' )
            
        #-----------------------------------------------------------------------------
        
        
        # Compute current time-step perturbation from base-state for subsequent output files
        #-----------------------------------------------------------------------------
        if( i > i_start ):
            
            print( f'\n\t\tComputing perterbations from base-state values... {sims[j]}, {(i-1)*dt} min' )
            # U wind pert (m/s)
            u_pert = u[0:km1,:,:].values - uinterp_bs
            uinterp_pert[i-2,:,:] = u_pert
            
            # V wind pert (m/s)
            v_pert = v[0:km1,:,:].values - vinterp_bs
            vinterp_pert[i-2,:,:] = v_pert
            
            # CAPE pert (J/kg)
            cape_p = cape.values - cape_bs
            cape_pert[i-2,:,:] = cape_p
            
            # CIN pert (J/kg)
            cin_p = cin.values - cin_bs
            cin_pert[i-2,:,:] = cin_p
    
            # SRH pert (m2/s2)
            srh1km_p = srh1km.base - srh1km_bs
            srh1km_pert[i-2,:,:] = srh1km_p
            
            # FrH pert (n/a)
            FrH_p = FrH.values - FrH_bs
            FrH_pert[i-2,:,:] = FrH_p
       
            # Calculate windward perturbation averages 
            #-----------------------------------------------------------------------------
            print( '\n\t\tCalculate Windward Perturbation Averages...' )
            ww_cape_pert_avg = np.mean( cape_p[wwY1:wwY2, wwX1:wwX2] )
            ww_cin_pert_avg = np.mean( cin_p[wwY1:wwY2, wwX1:wwX2] )
            ww_srh_pert_avg = np.mean( srh1km_p[wwY1:wwY2, wwX1:wwX2] )
            ww_FrH_pert_avg = np.nanmean( FrH_p[wwY1:wwY2, wwX1:wwX2] )
            print( f'\t\tWindward Average: Frh Pert = {ww_FrH_pert_avg}' )
            
            # Calculate leeward perturbation averages 
            #-----------------------------------------------------------------------------
            print( '\n\t\tCalculate Leeward Perturbation Averages...' )
            lw_cape_pert_avg = np.mean( cape_p[lwY1:lwY2, lwX1:lwX2] )
            lw_cin_pert_avg = np.mean( cin_p[lwY1:lwY2, lwX1:lwX2] )
            lw_srh_pert_avg = np.mean( srh1km_p[lwY1:lwY2, lwX1:lwX2] )
            lw_FrH_pert_avg = np.nanmean( FrH_p[lwY1:lwY2, lwX1:lwX2] )
            print( f'\t\tLeeward Average: Frh Pert = {lw_FrH_pert_avg}' )
        
        #-----------------------------------------------------------------------------
        # End: Data Collection and Manipulation
        #-----------------------------------------------------------------------------
    
    #-----------------------------------------------------------------------------
    # End: Output loop
    #-----------------------------------------------------------------------------
        
    
    #-----------------------------------------------------------------------------
    # Begin: Plotting Logic (After perterbations have been computed)
    #-----------------------------------------------------------------------------
        
        
    print( f'\n\t\tGenerating summary plot...\ti={i}' )
    #-----------------------------------------------------------------------------
    
    # Mask out lowest level of terrain (Needed for Realistic Terrain, Does not affect Idealized Terrain)
    zs_smooth = scipy.ndimage.gaussian_filter1d( zs, 10 )
    trn_mask = np.ma.array( zs_smooth, mask = zs < 0.01 )
    
    # Skip value array to subset wind vectors when plotting
    skip = ( slice( None, None, skip_val ), slice( None, None, skip_val ) ) 
    
    # Logic to determine if x-axis should be peak-relative or not
    #-----------------------------------------------------------------------------
    if( terrain == 'IDTRN' or terrain == 'MOD_IDTRN' ):
        
        # Convert x-coordinates from grid-relative to peak-relative coordinates
        print( '\n\t\t\tConverting to peak-relative coordinates for IDTRN...' )
        xh = xh - peak_pos * units( 'km' )
        
        # Reference vector location
        qx = -125
        qy =  235
        
        # Panel annotation location
        anX = 125
        anY = 212.5
        
        # Average annotation location
        avgX1 = -175
        avgX2 = 75
        avgY1 = 55
        avgY2 = 55
        
        # Change x-axis label to reflect coordinate shift
        xlabel_str = "Distance from Peak (km)" 
        
    # Label the axis as grid-relative    
    else:
        xlabel_str = "Zonal Distance (km)"
        
        # Reference vector location
        qx = 235
        qy = 325
        
        # Panel annotation location
        anX = 475
        anY = 320
        
        # Average annotation location
        avgX1 = 205
        avgX2 = 405
        avgY1 = 155
        avgY2 = avgY1 + 25

    #-----------------------------------------------------------------------------
    
    # Create a plotting grid based on the x,y coordinates (Has to come after conversion to plot correctly!)
    X,Y = np.meshgrid( xh.m, yh.m )
    
    #-----------------------------------------------------------------------------
    # Begin: Axes Plotting Loop
    #-----------------------------------------------------------------------------
    for ax in axes:
        
        print( f'\n\t\tAxis Plotting: i = {i}, j = {j}, ax = {ax}, {sims[j]}, {(i-1)*dt} min' ) 
        
        if( ax == 'M' or ax == 'N' or ax == 'O' or ax == 'T' ):
            axes[ax].axis( 'off' )
        
        # Universal axes configurations
        #-----------------------------------------------------------------------------
            
        # Add a grid
        axes[ax].grid( alpha = 0.5 )
        
        # Set a title for the two top panels
        if( ax == 'A' ):
            if( j == 0 ):
                    axes[ax].set_title( f'{sims[j].upper()}: Base-State (t= {t0} min)', fontsize = 22, fontweight = 'bold' )
        if(  ax == 'G' ):
            if( j == 1 ):
                axes[ax].set_title( f'{sims[j].upper()}: Base-State (t= {t0} min)', fontsize = 22, fontweight = 'bold' )
        
        if( ax == 'D' ):
            if( j == 0 ):
                    axes[ax].set_title( f'{sims[j].upper()}: Perturbations (t= {mtime} min)', fontsize = 20, fontweight = 'bold' )
        if( ax == 'J' ):
            if( j == 1 ):
                    axes[ax].set_title( f'{sims[j].upper()}: Perturbations (t= {mtime} min)', fontsize = 20, fontweight = 'bold' )
            
        # Set axes labels for left-hand panels
        if( ax == 'A' or ax == 'B' or ax == 'C' or ax == 'P' ):
            axes[ax].set_ylabel( "Meridional Distance (km)", fontsize = 20, fontweight = 'bold' )
        if( ax == 'P' or ax == 'Q' or ax == 'R' or ax == 'S' ):
            axes[ax].set_xlabel( xlabel_str, fontsize = 20, fontweight = 'bold' )
        
        # Remove inner ticklabels
        if( ax == 'A' or ax == 'B' or ax == 'C' or ax == 'D' or ax == 'E' or 
            ax == 'F' or ax == 'G' or ax == 'H' or ax == 'I' or ax == 'J' or 
            ax == 'K' or ax == 'L' ):
            axes[ax].xaxis.set_ticklabels([])
        if( ax == 'D' or ax == 'E' or ax == 'F' or ax == 'G' or ax == 'H' or 
            ax == 'I' or ax == 'J' or ax == 'K' or ax == 'L' or ax == 'M' or
            ax == 'N' or ax == 'O' or ax == 'Q' or ax == 'R' or ax == 'S'):
            axes[ax].yaxis.set_ticklabels([])
        #-----------------------------------------------------------------------------
            
               
        # Plot variables required on all panels
        #-----------------------------------------------------------------------------
        if( ax == 'A' or ax == 'B' or ax == 'C' or ax == 'D' or ax == 'E' or 
            ax == 'F' or ax == 'G' or ax == 'H' or ax == 'I' or ax == 'J' or 
            ax == 'K' or ax == 'L' or ax == 'P' or ax == 'Q' or ax == 'R' or ax == 'S' ):
            
            # Set plotting extent
            axes[ax].set_xlim( x1, x2 )
            axes[ax].set_ylim( y1, y2 )
            
            # Plot terrain contours 
            c_trn = axes[ax].contour( X,Y, trn_mask, levels = np.arange( 0, peak_elv + inc, inc ),
                                      colors = 'tab:brown', linewidths = 3, alpha = 0.75 )
            
            # Create averaging window patches 
            ww_avg_area = patches.Rectangle( xy = (xh[wwX1].m, yh[wwY1].m), width = (xh[wwX2].m - xh[wwX1].m),
                                             height = (yh[wwY2].m - yh[wwY1].m), fill = False,
                                             linestyle = '--', linewidth = 1.5,
                                             edgecolor = 'k', alpha = 0.75, zorder = 100  )
            
            lw_avg_area = patches.Rectangle( xy = (xh[lwX1].m, yh[lwY1].m), width = (xh[lwX2].m - xh[lwX1].m),
                                             height = (yh[lwY2].m - yh[lwY1].m),  fill = False,
                                             linestyle = '--', linewidth = 1.5,
                                             edgecolor = 'k', alpha = 0.75, zorder = 100  )
            
            # Add averaging window patches to first panel only
            axes[ax].add_patch( ww_avg_area )
            axes[ax].add_patch( lw_avg_area ) 
            
            if( terrain != 'RLTRN' ):
            
                # Plot terrain labels for IDTRN only (NEEDS TWEAKING)
                c_trn_lab = axes[ax].clabel( CS = c_trn, levels = np.arange( 0, peak_elv + inc, inc ),
                                              fontsize = 16, use_clabeltext = True, )
                
                # Add halos 
                plt.setp( [ c_trn_lab ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        #-----------------------------------------------------------------------------
                         
        
        
        
        # Plot variables required on all base-state panels
        #-----------------------------------------------------------------------------
        if( ax == 'A' or ax == 'B' or ax == 'C' or ax == 'P' ):
            if( j == 0 ):
            
                print( f'\n\t\tVector Base-State Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot wind vectors
                q = axes[ax].quiver( X[skip], Y[skip], uinterp_bs[0,:,:][skip], vinterp_bs[0,:,:][skip],
                                      color = 'k', units = 'xy', angles = 'uv', scale = v_scale,
                                      scale_units = 'xy', pivot = 'tail', alpha = 0.75, zorder = 100 )
                
                # Add wind reference vector just to panel A (NEEDS ADJUSTING to make larger and add halo)
                if( ax == 'A' ):
                    q_lab = axes[ax].quiverkey( Q = q, coordinates = 'data', X = qx, Y = qy, U = 10.0, label = "V: 10 m $s^{-1}$",
                                                labelpos = 'N', labelsep = 0.075, color = 'black', fontproperties = {'size':20} )
                    
                    # Add halos 
                    plt.setp( [ q_lab ], path_effects = 
                              [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
   
        #-----------------------------------------------------------------------------
        

        # Plot variables required on all perterbation panels
        #-----------------------------------------------------------------------------
        if( ax == 'D' or ax == 'E' or ax == 'F' or ax == 'Q'):
            if( j == 0 ):
           
                print( f'\n\t\tVector Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' )    
            
                # Plot perturbation wind vectors
                q2 = axes[ax].quiver( X[skip], Y[skip], uinterp_pert[i-2, 0,:,:][skip], vinterp_pert[i-2, 0,:,:][skip],
                                      color = 'k', units = 'xy', angles = 'uv', scale = v_scale,
                                      scale_units = 'xy', pivot = 'tail', alpha = 0.75, zorder = 100 )
                
                # Add perturbation wind reference vector just to panel D (NEEDS ADJUSTING to make larger and add halo)
                if( ax == 'D' ):
                    q_lab2 = axes[ax].quiverkey( Q = q2, coordinates = 'data', X = qx, Y = qy, U = 10.0, label = "$V_{Pert}$: 10 m $s^{-1}$",
                                                labelpos = 'N', labelsep = 0.075, color = 'black', fontproperties = {'size':20} )
                    
                    # Add halos 
                    plt.setp( [ q_lab2 ], path_effects = 
                              [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )

        #-----------------------------------------------------------------------------





        # Round 2 Plot variables required on all base-state panels
        #-----------------------------------------------------------------------------
        if(  ax == 'G' or ax == 'H' or ax == 'I' or ax == 'R' ):
            if( j == 1 ):
        
                print( f'\n\t\tVector Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot wind vectors
                q = axes[ax].quiver( X[skip], Y[skip], uinterp_bs[0,:,:][skip], vinterp_bs[0,:,:][skip],
                                      color = 'k', units = 'xy', angles = 'uv', scale = v_scale,
                                      scale_units = 'xy', pivot = 'tail', alpha = 0.75, zorder = 100 )
                
                # Add wind reference vector just to panel A (NEEDS ADJUSTING to make larger and add halo)
                if( ax == 'G' ):
                    q_lab = axes[ax].quiverkey( Q = q, coordinates = 'data', X = qx, Y = qy, U = 10.0, label = "V: 10 m $s^{-1}$",
                                                labelpos = 'N', labelsep = 0.075, color = 'black', fontproperties = {'size':20} )
                    
                    # Add halos 
                    plt.setp( [ q_lab ], path_effects = 
                              [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        
        #-----------------------------------------------------------------------------
        
        
        # Plot variables required on all perterbation panels
        #-----------------------------------------------------------------------------
        if( ax == 'J' or ax == 'K' or ax == 'L' or ax == 'S' ):
            if( j == 1 ):
            
                print( f'\n\t\tVector Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot perturbation wind vectors
                q2 = axes[ax].quiver( X[skip], Y[skip], uinterp_pert[i-2, 0,:,:][skip], vinterp_pert[i-2, 0,:,:][skip],
                                      color = 'k', units = 'xy', angles = 'uv', scale = v_scale,
                                      scale_units = 'xy', pivot = 'tail', alpha = 0.75, zorder = 100 )
                
                # Add perturbation wind reference vector just to panel D (NEEDS ADJUSTING to make larger and add halo)
                if( ax == 'J' ):
                    q_lab2 = axes[ax].quiverkey( Q = q2, coordinates = 'data', X = qx, Y = qy, U = 10.0, label = "$V_{Pert}$: 10 m $s^{-1}$",
                                                labelpos = 'N', labelsep = 0.075, color = 'black', fontproperties = {'size':20} )
                    
                    # Add halos 
                    plt.setp( [ q_lab2 ], path_effects = 
                              [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        
        #-----------------------------------------------------------------------------




        #-----------------------------------------------------------------------------
        # Begin: Round 1 Plot Base-State Variables (Panels A-F)               
        #-----------------------------------------------------------------------------

        # Panel A Plotting: Base-State CAPE   
        #-----------------------------------------------------------------------------
        if( ax == 'A' ):
            if( j == 0 ):
            
                print( f'\n\t\tCAPE Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot base-state CAPE
                cf_cape_bs = axes[ax].contourf( X,Y, cape_bs[:,:], levels = np.arange( cape_min, cape_max + cape_int, cape_int ),
                                                cmap = plt.cm.get_cmap( 'Reds' ), extend = 'both', alpha = 0.5 )
                
                # Add annotations
                an_A = axes[ax].annotate( 'a', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cape_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cape_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cape_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cape_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_A, ww_cape_ann, lw_cape_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        
        #-----------------------------------------------------------------------------
        
        
        # Panel B Plotting: Base-State CIN 
        #-----------------------------------------------------------------------------
        if( ax == 'B' ):
            if( j == 0 ):
            
                print( f'\n\t\tCIN Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot base-state CIN
                cf_cin_bs = axes[ax].contourf( X,Y, cin_bs[:,:], levels = np.arange( cin_min, cin_max + cin_int, cin_int ),
                                                cmap = plt.cm.get_cmap( 'Blues_r' ), extend = 'max', alpha = 0.5 )
                
                # Add panel annotation
                an_B = axes[ax].annotate( 'b', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cin_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cin_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cin_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cin_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_B, ww_cin_ann, lw_cin_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
                
        #-----------------------------------------------------------------------------
        
        
        # Panel C Plotting: Base-State SRH 
        #-----------------------------------------------------------------------------
        if( ax == 'C' ):
            if( j == 0 ):
            
                print( f'\n\t\tSRH Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot base-state SRH
                cf_srh1km_bs = axes[ax].contourf( X,Y, srh1km_bs[:,:], levels = np.arange( srh_min, srh_max + srh_int, srh_int ),
                                                cmap = plt.cm.get_cmap( 'Purples' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_C = axes[ax].annotate( 'c', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_srh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_srh_bs_avg,2))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_srh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_srh_bs_avg,2))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_C, ww_srh_ann, lw_srh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
                
        #-----------------------------------------------------------------------------
        
        # Panel P Plotting: Base-State FrH
        #-----------------------------------------------------------------------------
        if( ax == 'P' ):
            if( j == 0 ):
                
                print( f'\n\t\tFrH Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot base-state SRH
                cf_FrH_bs = axes[ax].contourf( X,Y, FrH_bs[:,:], levels = np.arange( FrH_min, FrH_max + Frh_int, Frh_int ),
                                                cmap = plt.cm.get_cmap( 'Greys' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_D = axes[ax].annotate( 'd', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_frh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_FrH_bs_avg,2))} ',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_frh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_FrH_bs_avg,2))} ',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_D, ww_frh_ann, lw_frh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
                    
        #-----------------------------------------------------------------------------
        
        #-----------------------------------------------------------------------------
        # End: Round 1 Plot Base-State Variables (Panels A-C)               
        #-----------------------------------------------------------------------------
        
        
        
        #-----------------------------------------------------------------------------
        # Begin: Round 1 Plot Perturbation Variables (Panels D-F)               
        #-----------------------------------------------------------------------------
        
        # Panel D Plotting: Perturbation CAPE  
        #-----------------------------------------------------------------------------
        if( ax == 'D' ):
            if( j == 0 ):
            
                print( f'\n\t\tCAPE Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' )  
                
                # Plot pert. CAPE
                cf_cape_pert = axes[ax].contourf( X,Y, cape_pert[i-2,:,:], levels = np.arange( capep_min, capep_max + capep_int, capep_int ),
                                                  cmap = plt.cm.get_cmap( 'coolwarm' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_E = axes[ax].annotate( 'e', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cape_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cape_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cape_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cape_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_E, ww_cape_ann, lw_cape_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        #-----------------------------------------------------------------------------
        
        
        # Panel E Plotting: Perturbation CIN   
        #-----------------------------------------------------------------------------
        if( ax == 'E' ):
            if( j == 0 ):
            
                print( f'\n\t\tCIN Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot pert. CIN
                cf_cin_pert = axes[ax].contourf( X,Y, cin_pert[i-2,:,:], levels = np.arange( cinp_min, cinp_max + cinp_int, cinp_int ),
                                                  cmap = plt.cm.get_cmap( 'coolwarm_r' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_F = axes[ax].annotate( 'f', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cin_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cin_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cin_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cin_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_F, ww_cin_ann, lw_cin_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        #-----------------------------------------------------------------------------
             
        
        # Panel F Plotting: Perturbation SRH 
        #-----------------------------------------------------------------------------
        if( ax == 'F'):
            if( j == 0 ):
            
                print( f'\n\t\tSRH Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
            
                # Plot pert. srh
                cf_srh1km_pert = axes[ax].contourf( X,Y, srh1km_pert[i-2,:,:], levels = np.arange( srhp_min, srhp_max + srhp_int, srhp_int ),
                                                cmap = plt.cm.get_cmap( 'coolwarm' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_G = axes[ax].annotate( 'g', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_srh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_srh_pert_avg,1))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_srh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_srh_pert_avg,1))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_G, ww_srh_ann, lw_srh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        #-----------------------------------------------------------------------------
        
        
        # Panel Q Plotting: Perturbation FrH
        #-----------------------------------------------------------------------------
        if( ax == 'Q' ):
            if( j == 0 ):
            
                print( f'\n\t\tFrH Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' )  
            
                # Plot pert FrH
                cf_FrH_pert = axes[ax].contourf( X,Y, FrH_pert[i-2,:,:], levels = np.arange( FrHp_min, FrHp_max + FrHp_int, FrHp_int ),
                                                cmap = plt.cm.get_cmap( 'coolwarm' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_H = axes[ax].annotate( 'h', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_frh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_FrH_pert_avg,2))} ',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_frh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_FrH_pert_avg,2))} ',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_H, ww_frh_ann, lw_frh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        #-----------------------------------------------------------------------------
        
        #-----------------------------------------------------------------------------
        # End: Round 1 Plot Perturbation Variables (Panels G-L)               
        #-----------------------------------------------------------------------------
        
        
        
        
        # ROUND 2: Required to ensure double plotting doesn't occur (I was to lazy to do this in a more pythonic way)
        #-----------------------------------------------------------------------------


        #-----------------------------------------------------------------------------
        # Begin: Plot Base-State Variables (Panels A-F)               
        #-----------------------------------------------------------------------------

        # Panel G Plotting: Base-State CAPE   
        #-----------------------------------------------------------------------------
        if( ax == 'G' ):
            if( j == 1 ):
            
                print( f'\n\t\tCAPE Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot base-state CAPE
                cf_cape_bs = axes[ax].contourf( X,Y, cape_bs[:,:], levels = np.arange( cape_min, cape_max + cape_int, cape_int ),
                                                cmap = plt.cm.get_cmap( 'Reds' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_I = axes[ax].annotate( 'i', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cape_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cape_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cape_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cape_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_I, ww_cape_ann, lw_cape_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
            
        if( ax == 'M' and j == 1 ):
            
            # Construct axis objects for base-state CAPE colorbar
            cbarA = fig.colorbar( cf_cape_bs, ax = axes[ax], fraction = 0.05, pad = 0.15 )

            # Set base-state CAPE colorbar label 
            cbarA.set_label( label = 'SBCAPE (J kg$^{-1}$)', fontproperties = { 'size':20, 'weight':'bold' } )
        
        #-----------------------------------------------------------------------------
        
        
        # Panel H Plotting: Base-State CIN 
        #-----------------------------------------------------------------------------
        if( ax == 'H' ):
            if( j == 1 ):
                
                print( f'\n\t\tCIN Base-State Plotting: {sims[j]}, {(i-1)*dt} min' )     
            
                # Plot base-state CIN
                cf_cin_bs = axes[ax].contourf( X,Y, cin_bs[:,:], levels = np.arange( cin_min, cin_max + cin_int, cin_int ),
                                                cmap = plt.cm.get_cmap( 'Blues_r' ), extend = 'min', alpha = 0.5 )
                
                # Add panel annotation
                an_J = axes[ax].annotate( 'j', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cin_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cin_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cin_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cin_bs_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_J, ww_cin_ann, lw_cin_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
            
        if( ax == 'N' and j == 1 ):
            
            # Construct axis objects for base-state CIN colorbar
            cbarB = fig.colorbar( cf_cin_bs, ax = axes[ax], fraction = 0.05, pad = 0.15 )

            # Set base-state CIN colorbar label 
            cbarB.set_label( label = 'SBCIN (J kg$^{-1}$)', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
        
        
        # Panel I Plotting: Base-State SRH 
        #-----------------------------------------------------------------------------
        if( ax == 'I' ):
            if( j == 1 ):
                
                print( f'\n\t\tSRH Base-State Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot base-state SRH
                cf_srh1km_bs = axes[ax].contourf( X,Y, srh1km_bs[:,:], levels = np.arange( srh_min, srh_max + srh_int, srh_int ),
                                                cmap = plt.cm.get_cmap( 'Purples' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_K = axes[ax].annotate( 'k', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_srh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_srh_bs_avg,1))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_srh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_srh_bs_avg,1))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_K, ww_srh_ann, lw_srh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
                
        if( ax == 'O' and j == 1  ):
            
            # Construct axis objects for base-state SRH colorbar
            cbarC = fig.colorbar( cf_srh1km_bs, ax = axes[ax], fraction = 0.05, pad = 0.15 )

            # Set base-state SRH colorbar label 
            cbarC.set_label( label = 'SRH1KM (m$^{2}$ s$^{-2}$)', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
        
        
        # Panel R Plotting: Base-State FrH
        #-----------------------------------------------------------------------------
        if( ax == 'R' ):
            if( j == 1 ):
               
                print( f'\n\t\tFrH Base-State Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot base-state SRH
                cf_FrH_bs = axes[ax].contourf( X,Y, FrH_bs[:,:], levels = np.arange( FrH_min, FrH_max + Frh_int, Frh_int ),
                                                cmap = plt.cm.get_cmap( 'Greys' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_L = axes[ax].annotate( 'l', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_frh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_FrH_bs_avg,2))} ',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_frh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_FrH_bs_avg,2))} ',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_L, ww_frh_ann, lw_frh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
            
        if( ax == 'T' and j == 1  ):
            
            # Construct axis objects for pert. srh colorbar
            cbarQ = fig.colorbar( cf_FrH_bs, ax = axes[ax], fraction = 0.05, pad = 0.15 )

            # Set pert. srh colorbar label 
            cbarQ.set_label( label = 'FrH1KM', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
        
        #-----------------------------------------------------------------------------
        # End: Plot Base-State Variables (Panels A-C)               
        #-----------------------------------------------------------------------------
        
        
        
        #-----------------------------------------------------------------------------
        # Begin: Plot Perturbation Variables (Panels D-F)               
        #-----------------------------------------------------------------------------
        
        # Panel J Plotting: Perturbation CAPE  
        #-----------------------------------------------------------------------------
        if( ax == 'J' ):
            if( j == 1 ):
            
                print( f'\n\t\tCAPE Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot pert. CAPE
                cf_cape_pert = axes[ax].contourf( X,Y, cape_pert[i-2,:,:], levels = np.arange( capep_min, capep_max + capep_int, capep_int ),
                                                  cmap = plt.cm.get_cmap( 'coolwarm' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_M = axes[ax].annotate( 'm', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cape_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cape_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cape_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cape_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_M, ww_cape_ann, lw_cape_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
            
        if( ax == 'M' and j == 1 ):
            
            # Construct axis objects pert. CAPE colorbar
            cbarD = fig.colorbar( cf_cape_pert, ax = axes[ax], fraction = 0.5, pad = 0.05 )

            # Set pert. CAPE colorbar label 
            cbarD.set_label( label = '$SBCAPE_{Pert}$ (J kg$^{-1}$)', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
        
        
        # Panel K Plotting: Perturbation CIN   
        #-----------------------------------------------------------------------------
        if( ax == 'K' ):
            if( j == 1 ):
            
                print( f'\n\t\tCIN Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot pert. CIN
                cf_cin_pert = axes[ax].contourf( X,Y, cin_pert[i-2,:,:], levels = np.arange( cinp_min, cinp_max + cinp_int, cinp_int ),
                                                  cmap = plt.cm.get_cmap( 'coolwarm_r' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_N = axes[ax].annotate( 'n', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_cin_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_cin_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_cin_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_cin_pert_avg,1))} ' + '(J kg$^{-1}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_N, ww_cin_ann, lw_cin_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        
        if( ax == 'N' and j == 1  ):
            
            # Construct axis objects pert. CIN colorbar
            cbarE = fig.colorbar( cf_cin_pert, ax = axes[ax],  fraction = 0.5, pad = 0.05 )

            # Set pert. CIN colorbar label 
            cbarE.set_label( label = '$SBCIN_{Pert}$ (J kg$^{-1}$)', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
             
        
        # Panel L Plotting: Perturbation SRH 
        #-----------------------------------------------------------------------------
        if( ax == 'L' ):
            if( j == 1 ):
            
                print( f'\n\t\tSRH Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
            
                # Plot pert. srh
                cf_srh1km_pert = axes[ax].contourf( X,Y, srh1km_pert[i-2,:,:], levels = np.arange( srhp_min, srhp_max + srhp_int, srhp_int ),
                                                cmap = plt.cm.get_cmap( 'coolwarm' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_O = axes[ax].annotate( 'o', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_srh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_srh_pert_avg,1))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_srh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_srh_pert_avg,1))} ' + '(m$^{2}$ s$^{-2}$)',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_O, ww_srh_ann, lw_srh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
            
        if( ax == 'O' and j == 1  ):
            
            # Construct axis objects for pert. srh colorbar
            cbarF = fig.colorbar( cf_srh1km_pert, ax = axes[ax], fraction = 0.5, pad = 0.05 )

            # Set pert. srh colorbar label 
            cbarF.set_label( label = '$SRH1KM_{Pert}$ (m$^{2}$ s$^{-2}$)', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
        
        
        # Panel S Plotting: Perturbation FrH
        #-----------------------------------------------------------------------------
        if( ax == 'S' ):
            if( j == 1 ):
                
                print( f'\n\t\tFrH Perturbation Plotting: {sims[j]}, {(i-1)*dt} min' ) 
                
                # Plot pert. srh
                cf_FrH_pert = axes[ax].contourf( X,Y, FrH_pert[i-2,:,:], levels = np.arange( FrHp_min, FrHp_max + FrHp_int, FrHp_int ),
                                                cmap = plt.cm.get_cmap( 'coolwarm' ), extend = 'both', alpha = 0.5 )
                
                # Add panel annotation
                an_P = axes[ax].annotate( 'p', xy = (anX, anY),  fontsize = 72, fontweight = 'bold' )
                ww_frh_ann = axes[ax].annotate( f'Windward Avg \n{str(round(ww_FrH_pert_avg,2))} ',
                                                 xy = (avgX1, avgY1), fontsize = 18, fontweight = 'bold' )
                lw_frh_ann = axes[ax].annotate( f'Leeward Avg \n{str(round(lw_FrH_pert_avg,2))} ',
                                                 xy = (avgX2, avgY2), fontsize = 18, fontweight = 'bold' )
                
                # Add halos 
                plt.setp( [ an_P, ww_frh_ann, lw_frh_ann ], path_effects = 
                          [ PathEffects.withStroke( linewidth = 5, foreground = 'white', alpha = 0.75 ) ] )
        
        if( ax == 'T' and j == 1  ):
            
            # Construct axis objects for pert. srh colorbar
            cbarQ = fig.colorbar( cf_FrH_pert, ax = axes[ax],  fraction = 0.5, pad = 0.05 )

            # Set pert. srh colorbar label 
            cbarQ.set_label( label = '$FrH1KM_{Pert}$', fontproperties = { 'size':20, 'weight':'bold' } )
        #-----------------------------------------------------------------------------
        
        
        #-----------------------------------------------------------------------------
        # End: ROUND 2 Plot Perturbation Variables (Panels G-L)               
        #-----------------------------------------------------------------------------

        
    #-----------------------------------------------------------------------------
    # End: Axes Plotting Loop
    #-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# End: Simulation loop
#-----------------------------------------------------------------------------

# Save the current plot
fig.savefig( fname = plot_dir + "/{}_env_pert_panel_plot_w_FrH".format( terrain ), dpi = 300 )
            # bbox_inches = "tight" )

# Confirm that script successfully ran
print( "\n{} successfully completed!".format( program ) )

# Report the time required to run the function
print( "\nScript Total Runtime: {}".format( datetime.now() - startTime ) )
