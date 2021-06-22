
#####################################################################
#####################################################################

import numpy as np
import scipy.spatial

#####################################################################
#####################################################################


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