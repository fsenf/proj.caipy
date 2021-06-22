#####################################################################
#####################################################################

def calculate_COP(ObjCX , ObjCY , ObjArea):
    
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
    '''

    # Number of clusters
    N_c = len( ObjCX )

    # Array of shape (N_c,2), containing the coordinates in 
    coords_clm = np.transpose(np.vstack((ObjCX, ObjCY))) # [km]

    # Euclidean Distance between each pair of the two input arrays
    # Output: distance of pairs in a symmetric and trace-less matrix as 
    # two input arrays are identical
    dist_clm = scipy.spatial.distance.cdist(coords_clm,coords_clm, 'euclidean') # [km]
    NumPermu = N_c * (N_c-1) / 2

    V=np.zeros(NumPermu)
    V[:]=np.nan
    Ind=0

    for i in range(0, N_c):
        for j in range(i+1, N_c):
            k      = (N_c*(N_c-1)/2) - (N_c-i)*((N_c-i)-1)/2 + j - i - 1 # index of upper diagonal elements https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
            Vij    = (np.sqrt(ObjArea[i])+np.sqrt(ObjArea[j])) / ( dist_clm[i,j] * np.sqrt(np.pi) )
            V[Ind] = Vij
            Ind    =Ind+1

    COP=np.sum(V)/NumPermu
    return COP
