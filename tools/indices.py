#!/usr/bin/env python


# Load libraries ---------------------------------------------------------
import numpy as np
from skimage import measure
import math
import scipy
import sys
import os
import scipy.spatial
from utility_functions import stats, nne_size, boot, boot_obj_size


##########################################################################
##########################################################################

def calculate_Ishape(cluster_field, full_output = False):
 

    '''
    SHAPE INDEX - I.shape is defined as sqrt(A)/(0.282*P) and 
    0<= I_shape <= 1.0. I_shape=0 is the shape of a 
    line which has zero area and I_shape=1 is the shape of a perfect circle 
    => 0.5 is something like an ellipsoid

    The area A and perimeter of the objects are computed as follows:
    A=object's number pixels * pixel's area 
    P=contour line through the centers of border pixels * pixel's length. 
    Pixels are assumed to be squared.
 

    INPUT
    =====
    cluster_field: labeled / segmented 2d field
    
    OUTPUT
    ======
    I_shape_mean: average shape index
    I_shape: shape index for each object in case full_output=True
    
    '''        

    # Use a build-in function from scikit
    properties = measure.regionprops(cluster_field)

    # object's number of pixels (A) and the contour line through the centers of border pixels using a 4-connectivity (peri) 
    A = [prop.area for prop in properties]  
    peri = [prop.perimeter for prop in properties] 
     

    # Shape index 
    factor= 1./(2.*(math.pi)**(0.5)) # 0.282
    I_shape = np.sqrt(np.asarray(A))/(factor*np.asarray(peri))

    # Clip values larger than 1
    I_shape[np.where(I_shape > 1.0)] = 1.0

    # Mean over all objects
    I_shape_mean = np.mean(I_shape)
    
    if full_output:
        return I_shape_mean, I_shape
    else:
        return I_shape_mean


##########################################################################
##########################################################################

def calculate_Iorg(x_clm, y_clm, r_equ, grid, 
                   max_dist = 1000. , 
                   dist_interval = 1.,
                   number_of_repetitions = 100, 
                   full_output = False):
         
    '''
    Calculates the organization index (Tompkins and Semi (2017), Nair et al. (1998))
    
    INPUT
    =====
    x_clm: x-coordinates of objects in meters 
    y_clm: y-coordiantes of objects in meters
    r_equ: equivalent radii of objects
    grid: meshgrid  where objects are defined in
    max_dist: maximal distance in km for which the NNCDF is computed 
    dist_interval: bin width in km for the distance vector
    number_of_repetitions: number of repetitions for the bootstraping
    full_output: full_output=True gives NNCDFs
    

    OUTPUT
    ======
    Iorg: organization index
    UP: Iorg + 1.96STD
    LW: Iorg - 1.96STD
    dist: distances in km for which the NNCDF is computed
    Fo: observed nearest neighbor cumulative distribution function
    Fs: simulated nearest neighbor cumulative distribution function
    Fm: mean of nreps NNCDFs (the inhibition NNCDF)
    Fmh: Fm + 1.96STD
    Fml: Fm - 1.96STD
    '''    

    # number of repetitions for the bootstrapping
    nreps = number_of_repetitions 

    # Number of clusters
    N_c = len( x_clm )

    # get external grid
    x_lon, y_lat = grid
    
    # Array of shape (N_cl,2), containing the coordinates in meters
    coords_clm = np.transpose(np.vstack((x_clm, y_clm)))

    # Euclidean Distance between each pair of the two input arrays
    # Output: distance of pairs in a symmetric and trace-less matrix as 
    # two input arrays are identical
    # The distance is given in meters
    dist_clm = scipy.spatial.distance.cdist(coords_clm,coords_clm, 'euclidean')


    # Change dist_clm from m to km
    dist_clm = dist_clm/1000.
   


    # Adapt the distance to acount for the size of the 
    # objects (Iedas function nne_size reduces the distance based on centers
    # of mass by the equivalent radii of the two objects considered)
    # nne_size also replaces trace of the matrix with nan, trace elements contain 
    # zero distance as they reflect the distance between the same object
    dist_size = nne_size(dist_clm, r_equ)
    
     
    # nne_size can lead to negative distances in case the equivalent radius is larger than the distance 
    # occurs if two objects are closer to each other than their equivalent radii
    # These distances are then replaced by NAN. 
    dist_size[np.where(dist_size <= 0.)] = np.nan 


    if (len(np.where(dist_size <= 0.)[0]) > 0.):
       sys.exit('dist_size is <= zero !')
    

    # Number of elements in dist_size not containing nan
    n_dist_size = np.count_nonzero(~np.isnan(dist_size))


    # Distances (in km) for which the NNCDF is computed
    # Bin width = dist_interval km

    dist = np.arange( 0., max_dist, dist_interval) 


    # Calculate NNCDF only in case that dist_size array contains at least two objects
    # and not only NaN entries

    if n_dist_size > 0:

        # NN distance
        nnd=np.zeros(N_c)
	nnd[:]=np.nan
	# only NN distance is retained
	for ob in range(N_c):
	    nnd[ob]=np.sort(dist_size[ob,:])[0]
	    nnaux=nnd[ob]
            


        # NNCDF of the objects based on observations (fractional number of objects whose NN distance is <=r)
        frac=[]
        for fr in dist:
            fra = np.count_nonzero(nnd <=fr)/(float(nnd.shape[0]))
            frac.append(fra)

        # Observed NNCDF     
        Fo = np.asarray(frac)    

        # ####################################################################################  
        # BOOTSTRAPING: statistical method to resample
        # ####################################################################################

        print '   BOOTSTRAPING'
       
        print 'number of repetitions = ', nreps 

        frac_sim = []
        reps_count = 0
        new_N = 0


        for i in range(nreps):
          
            # Avoid to have only one single object in dist_sim, so draw new sample in that case.
            # Might occur in case total number of objects is very low (2-5)
            while new_N < 2:
                # draws samples from a Poisson distribution with lambda=N_c, 
                # lambda=variance and expected value
                new_N = np.random.poisson(N_c) # take out one sample

             
            # All simulations with the same number of objects
            #new_N=N_c

            # generates a random sample of equivalent radius (siz) out of the r_equ array 
            # (takes out new_N values from the r_equ-array)
            # (one object might be taken out several times)
            siz = np.random.choice(r_equ,new_N)  


            # The centers of mass of the objects are randomly distributed in space without overlap. 
            # Each object in siz gets a randomly chosen new position (meaning new lat/lon coordinates) 
            # inside the RADOLAN domain
            # Compute the distance between all possible pairs of objects based on the new positions
            # Replace distances = 0 by NAN (distance between the same objects)
            # Compute distance considering the equivalent radii of the objects

            dist_size_sim = boot_obj_size(x_lon.flatten(),y_lat.flatten(),new_N, siz)
            

            # Care for potential negative distances due to overlap of objects' areas
            # Replace negative distances by NAN
            dist_size_sim[np.where(dist_size_sim<=0.)] = np.nan   

            if (len(np.where(dist_size_sim <= 0.)[0]) > 0.):
               sys.exit('dist_size_sim is <= zero !')


            # Number of elements in dist_sim not containing nan
            n_dist_size_sim = np.count_nonzero(~np.isnan(dist_size_sim))
         

            # ################################################################################   
            # Inhibition NNCDF with 95% confidence interval
            # "inhibition" comes from Poisson point statistics (point process 
            # for which the presence of one point decreases the probability of finding 
            # neighbouring points in its vicinity (Weger et al.)
            # ################################################################################

            # Calculate NNCDF only in case that dist_size array contains at least two objects
            # and not only NaN entries
            

            if n_dist_size_sim > 0:

                #NN distance  
                Snnd=np.zeros(new_N)
		Snnd[:]=np.nan
		#only NN distance is retained
		for ob in range(new_N):
		    Snnd[ob]=np.sort(dist_size_sim[ob,:])[0]
		    Snnaux=Snnd[ob]
                   


                # NNCDF of the objects based on bootstrapping (fractional number of objects 
                # whose NN distance is <=r(=dist))
                for fr in dist:
                    fra = np.count_nonzero(Snnd <=fr)/(float(Snnd.shape[0]))
                    frac_sim.append(fra)

                reps_count += 1
           
            else: 

                for fr in dist:
                    fra = np.nan
                    frac_sim.append(fra)
                   
                reps_count += 1
 
        nreps_final = reps_count
        if nreps_final != number_of_repetitions:
            print("reps count ",reps_count)

        # Store all nreps NNCDFs in one array
        frac_sim = np.reshape(np.asarray(frac_sim),(nreps_final,np.shape(dist)[0]))
        
        # Save the simulated NNCDFs 
        Fs=frac_sim

        # Mean of nreps NNCDFs is the inhibition NNCDF
        Fm = np.mean(frac_sim, axis=0)

        # 95% confidence interval: 
        # with a certainty of 95% the sample lies in between +- 1.96*standard deviation
        Fmh=Fm + 1.96*np.std(frac_sim,axis=0)/(np.sqrt(float(nreps_final)))
        Fml=Fm - 1.96*np.std(frac_sim,axis=0)/(np.sqrt(float(nreps_final))) 

        # ####################################################################################  
        # ORGANISATION INDEX - I.org
        # ####################################################################################  
        print '   ORGANISATION INDEX'
        # Integrate the area under the curve of inhibition NNCDF vs. "observed" NNCDF
        Ior = np.trapz(Fo, Fm)
        UP = np.trapz(Fo, Fmh)
        LW = np.trapz(Fo, Fml)


        I2 = '{:04.2f}'.format(Ior)
        print '   I_org = ', I2 
        print '   UP = ', UP
        print '   LW = ', LW

    else:
        Ior = np.nan
        UP = np.nan
        LW = np.nan
        Fm=np.empty([len(dist)])
        Fm[:]=np.nan
        Fml=np.empty([len(dist)])
        Fml[:]=np.nan
        Fmh=np.empty([len(dist)])
        Fmh[:]=np.nan 
        Fo= np.empty([len(dist)])
        Fo[:]=np.nan
        Fs= np.empty([nreps,len(dist)])         
        Fs[:,:]=np.nan

    print 'shape Fo =', Fo.shape 
    print 'shape Fs =', Fs.shape
   
    if full_output:
        return Ior, UP, LW, dist, Fo, Fs, Fm, Fmh, Fml
    else:
        return Ior


#####################################################################
#####################################################################

def calculate_SCAI(x_clm, y_clm,  N_max = 100, L = 100):
    
    '''
    SCAI according to Tobin et al., 2012
    Based on the all nearest-neighbour (NN)distances between 
    the centers of mass for each pair of objects

    INPUT
    =====
    x_clm: x-positions of the objects
    y_clm: y-positions of the objects
    N_max: Potential maximum number of objects in the domain
    L    : Length in km of the characteristic domain size
    
    OUTPUT
    ======
    indx: SCAI index
    d0_clm: geometric mean of distances
    d1_clm: arithmetic mean of distances
    '''

    # Number of clusters
    N_c = len( x_clm )

    # Array of shape (N_c,2), containing the coordinates in meters
    coords_clm = np.transpose(np.vstack((x_clm, y_clm)))

    # Euclidean Distance between each pair of the two input arrays
    # Output: distance of pairs in a symmetric and trace-less matrix as 
    # two input arrays are identical
    # The distance is given in meters
    dist_clm = scipy.spatial.distance.cdist(coords_clm,coords_clm, 'euclidean')

    # Change dist_clm from m to km
    dist_clm = dist_clm/1000.
    
    print('   SCAI')

    # Extract the distances from the dist_clm symmetric matrix to form a vector
    # Indices corresponding to the upper triangular + diagonal of dist_clm
    ind_up_di_dist_clm = np.triu_indices_from(dist_clm)

    # Writing the distances in a vector
    dist_clm_vec = np.asarray(dist_clm[ind_up_di_dist_clm])

    # Removing the diagonal elements (distances = 0)
    dist_clm_vec = np.delete(dist_clm_vec,np.where(dist_clm_vec == 0))

    # Compute order-zero diameter D0 and order-one diameter D1 in km
    # which is the geometrical/arithmetical mean od distances in pairs
    d0_clm = scipy.stats.mstats.gmean(dist_clm_vec)
            
    d1_clm = np.mean(dist_clm_vec)
            
    print '   DO = {} km'.format(d0_clm)
    print '   D1 = {} km'.format(d1_clm)

    # Calculate SCAI
    indx = float(N_c)/float(N_max)* (d0_clm/L)*1000.            

    print '   SCAI = {}'.format(indx)

    return indx, d0_clm, d1_clm

##############################################################################
##############################################################################

def calculate_COP_fast(ObjCX , ObjCY , ObjArea):
    
    '''
    COP = Convective organization potential according to White et al., 2018, JAS
    Based on the all nearest-neighbour (NN) distances between the centers of mass for each pair of objects
    Impromvment of SCAI because the area of each object in taken into account to address the "interaction potenial"

    
    INPUT
    =====
    INPUT
    ObjCX:   x-positions of the objects [km]
    ObjCY:   y-positions of the objects [km]
    ObjArea: Object area [km]

    
    OUTPUT
    ======
    COP: COP index  []
    
    
    HISTORY
    =======
    2018-03-08: first version by Matthias Brueck
    2018-03-09: modifications (vectorization) by Fabian Senf
    '''

    # Number of clusters
    N_c = len( ObjCX )

    # Array of shape (N_c,2), containing the coordinates in 
    coords_clm = np.transpose(np.vstack((ObjCX, ObjCY))) # [km]

    # Euclidean Distance between each pair of the two input arrays
    # Output: distance of pairs in a symmetric and trace-less matrix as 
    # two input arrays are identical
    distance = scipy.spatial.distance.cdist(coords_clm, coords_clm, 'euclidean') # [km]
    
    
    # define the Radius of Objects
    radius = np.sqrt( ObjArea  / np.pi )
    
    
    # the White et al. (2018) COP Defintion has a summed adius matrix in the nominator
    radius_matrix, radius_matrix_T = np.meshgrid( radius, radius)
    
    # per-pair potential
    per_pair_potential = np.ma.divide( radius_matrix + radius_matrix_T, distance ) # masked out distance zero
    
    
    # define the upper triangle 
    #  - the diagonal includes invalid "self" pairing
    #  - the lower triangle has pairs already counted in the upper triangle 
    m = np.triu_indices(N_c, k = 0)

    
    # now we just take the average
    COP = np.ma.mean( per_pair_potential[m] )
    
    return COP





