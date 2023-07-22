#Branden Katona
#
#srh_cy.pyx
#Cython implementation of code to calcualte SRH since it involves a for-loop
#heavy approach due to terrain following coordinates and underlying terrain.
#
#u = 3d u winds
#v = 3d v winds
#zh = 3d field of height coordinates for each model grid
#zs = height of lower boundary
#u_storm = u_storm motion input by user
#v_storm = v component of storm motion input by user
#layer = depth over which we calculate SRH
#
#output 
# srh_out = srh over some layer specified by user
#
#


def srh (float[:,:,:] u, float[:,:,:] v, double[:,:,:] zh, float[:,:]zs, double u_storm, double v_storm, float layer):
    #Define loop control parameters based on input shapes
    cdef int ydim = u.shape[1]
    cdef int xdim = u.shape[2]
    cdef int zdim = u.shape[0]
    cdef int i,j,k
    import numpy as np
    #Pre-define return array for speed
    cdef double[:,:] srh_out = np.zeros((ydim,xdim))
    #Subtract out terrain field from height coords to get SRH AGL instead of 
    #above mean sea level
    for i in range(0,ydim):
        for j in range(0,xdim):
            for k in range(0,zdim):
                #If next model level is greater than the depth we want to 
                #integrate over, then interpolate to the final layer depth
                #and calculate the last part of the SRH calculation
                flag = 0
                if (zh[k+1,i,j]-zs[i,j])  > layer :
                    u[k+1,i,j] = u[k,i,j]+ (layer - (zh[k,i,j]-zs[i,j])) * (u[k+1,i,j] - u[k,i,j]) / (zh[k+1,i,j] - zh[k,i,j])
                    v[k+1,i,j] = v[k,i,j]+ (layer - (zh[k,i,j]-zs[i,j])) * (v[k+1,i,j] - v[k,i,j]) / (zh[k+1,i,j] - zh[k,i,j])
                    flag = 1

    
                #perform srh calculation here based on Markowski + Richardson 
                #2010
                srh_out[i,j] += ((u[k+1,i,j]-u_storm)*(v[k,i,j] - v_storm)) - ((u[k,i,j] - u_storm)*(v[k+1,i,j]-v_storm))
                if flag == 1:
                    break
    return srh_out
