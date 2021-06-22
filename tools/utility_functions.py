#!/usr/bin/env python

# MIUB
# Auf dem Huegel 20
# 53121 Bonn

# Author: Ieda Pscheidt
# Last update: 19.01.2018

###################################################################################
# FUNCTIONS USED IN THE COMPUTATION OF ORGANISATION INDICES
# #################################################################################
import numpy as np
import scipy 
 
# ##################################################################################
# BOOTSTRAPPING
# ##################################################################################

def boot(lon,lat,n):
  
    '''
    Objects are randomly distributed in space (area covered by radolan) without overlap.
    
    USAGE:
    position = boot(lon,lat,n)
    
    INPUT:
    lon= longitude of the objects
    lat= latitude of the objects
    n = number of objects
    
    OUTPUT:
    position= new positions (lon/lat) for the objects    
    
    '''
 
    x1= np.random.choice(a=range(len(lon)), size=n, replace=False)   
    new_lon=lon[x1]
    new_lat=lat[x1]
    coords=np.vstack((new_lon,new_lat)).T
    
    return(coords)
    
    
# ###################################################################################
# NEAREST NEIGHBOR DISTANCE
# ###################################################################################
def nne_size(dist,re):

    '''
    Distance between all possible pairs of objects considering the effective radius of each object.
    
    USAGE:
    distance = nne_size(dist,re)
    
    INPUT:
    dist= distance between centers of mass of the objects (considered as points)
    re= effrective radius (radius of a circle with the same area as the object)

    
    OUTPUT:
    distance= distance considering the size of the objects
    
    '''

    x=dist.shape[1]
    y=dist.shape[0]
    n=len(re)
    
    new_dist=[]
    for i in range(y):
       for j in range(x):
          if dist[i,j] != 0:
             new=dist[i,j] - re[i] -re[j]
          else:
             new=np.nan
             
          new_dist.append(new)  
          
    return(np.asarray(new_dist).reshape(y,x))

# ###################################################################################
# BOOTSTRAPPING AND COMPUTATION OF NN DISTANCE CONSIDERING R_EQU OF OBJECTS  
# ###################################################################################

def boot_obj_size(lon,lat,n,siz):

    '''
    Objects are randomly distributed in space (area covered by radolan) without overlap.
    Compute the distance between all possible pairs of objects based on the new positions
    and considering the r_equ of the objects
    
    
    USAGE:
    position = boot_obj_sitze(lon,lat,n,siz)
    
    INPUT:
    lon= longitude of the objects
    lat= latitude of the objects
    n = number of objects
    siz = size of the objects
    
    OUTPUT:   
    dist_size_sim = distance between all possible pairs 
                    of objects based on the new positions
                    and on the size of the objects
    '''
   
     
    # redistribution of the objects in space 
    x1= np.random.choice(a=range(len(lon)), size=n, replace=False)
    new_lon=lon[x1]
    new_lat=lat[x1]
    coords=np.vstack((new_lon,new_lat)).T

    # distance between all possible pairs of objects based on the new positions
    dist = scipy.spatial.distance.cdist(coords,coords,'euclidean')

    # change distance from m to km
    dist_km = dist/1000.

    # Replace distances = 0 by NAN (distance between the same objects)
    dist_sim = dist_km
    dist_sim[np.where(dist_sim == 0.)]=np.nan

    if (len(np.where(dist_sim <0.)[0]) > 0.):
       sys.exit('dist_sim is negative !')

    # distance considering the size of the objects (nan is already set in nne)
    dist_size_sim = nne_size(dist_sim, siz)


    return(dist_size_sim)


# ####################################################################################
# SUM OF THE VALUES IN AN ARRAY
# ####################################################################################

def sum_label(input, labels=None, index=None):
  
    '''
    Calculate the sum of the values of the array.
   
    Parameters
    ----------
    input : array_like
    Values of `input` inside the regions defined by `labels` are summed together.
    labels : array_like of ints, optional
    Assign labels to the values of the array. Has to have the same shape as `input`.
    index : array_like, optional
    A single label number or a sequence of label numbers of the objects to be measured.

    Returns
    -------
    sum : ndarray or scalar
    An array of the sums of values of `input` inside the regions defined by `labels` with the same shape as `index`. 
    If 'index' is None or scalar, a scalar is returned.

    Examples
    --------
    >>> input = [0,1,2,3]
    >>> labels = [1,1,2,2]
    >>> sum(input, labels, index=[1,2])
    [1.0, 5.0]
    >>> sum(input, labels, index=1)
    1
    >>> sum(input, labels)
    6
    '''
    
    count, sum = stats(input, labels, index)
    return sum
  
# ##################################################################################
# NUMBER OF VALUES/PIXELS IN EACH REGION LABELED BY INDEX
# ##################################################################################

def count_label(input, labels=None, index=None):
    
    '''
    Calculate the number of points/pixels of the different
       regions labeled by the labels in `index`.

    Parameters
    ----------
    input : array_like
       Array on which to compute the mean of elements over distinct
       regions.
    labels : array_like, optional
       Array of labels of same shape, or broadcastable to the same shape as
      `input`. All elements sharing the same label form one region over
       which the mean of the elements is computed.
    index : int or sequence of ints, optional
       Labels of the objects over which the mean is to be computed.
       Default is None, in which case the mean for all values where label is
       greater than 0 is calculated.

    Returns
    -------
    count : list
       Sequence of same length as `index`, with the number of points/pixels of the different
       regions labeled by the labels in `index`.

    Examples
    --------
    >>> a = np.arange(25).reshape((5,5))
    >>> labels = np.zeros_like(a)
    >>> labels[3:5,3:5] = 1
    >>> index = np.unique(labels)
    >>> labels
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [0, 0, 0, 1, 1],
           [0, 0, 0, 1, 1]])
    >>> index
    array([0, 1])
    >>> ndimage.mean(a, labels=labels, index=index)
    [10.285714285714286, 21.0]
    '''
    
    count, sum = stats(input, labels, index)
 #   return (sum / np.asanyarray(count).astype(np.float), count)
    return (count)

# ###############################################################################
#  COUNT AND SUM OF INPUT BY LABEL
# ###############################################################################
 
 
def stats(input, labels=None, index=None, centered=False):
    
    '''
    Count, sum, and optionally compute (sum - centre)^2 of input by label

    Parameters
    ----------
    input : array_like, n-dimensional
       The input data to be analyzed.
    labels : array_like (n-dimensional), optional
       The labels of the data in `input`. This array must be broadcast
       compatible with `input`; typically it is the same shape as `input`.
       If `labels` is None, all nonzero values in `input` are treated as
       the single labeled group.
    index : label or sequence of labels, optional
       These are the labels of the groups for which the stats are computed.
       If `index` is None, the stats are computed for the single group where
       `labels` is greater than 0.
    centered : bool, optional
       If True, the centered sum of squares for each labeled group is
       also returned. Default is False.

    Returns
    -------
    counts : int or ndarray of ints
       The number of elements in each labeled group.
    sums : scalar or ndarray of scalars
       The sums of the values in each labeled group.
    sums_c : scalar or ndarray of scalars, optional
       The sums of mean-centered squares of the values in each labeled group.
       This is only returned if `centered` is True.

    '''
    
    def single_group(vals):
       if centered:
          vals_c = vals - vals.mean()
          return vals.size, vals.sum(), (vals_c * vals_c.conjugate()).sum()
       else:
          return vals.size, vals.sum()

    if labels is None:
       return single_group(input)

    # ensure input and labels match sizes
    input, labels = np.broadcast_arrays(input, labels)

    if index is None:
       return single_group(input[labels > 0])

    if np.isscalar(index):
       return single_group(input[labels == index])
    
    def _sum_centered(labels):
       # `labels` is expected to be an ndarray with the same shape as `input`.
       # It must contain the label indices (which are not necessarily the labels
       # themselves).
       means = sums / counts
       centered_input = input - means[labels]
       # bincount expects 1d inputs, so we ravel the arguments.
       bc = np.bincount(labels.ravel(),weights=(centered_input *centered_input.conjugate()).ravel())
       return bc
    
    # Remap labels to unique integers if necessary, or if the largest
    # label is larger than the number of values.

    if (not _safely_castable_to_int(labels.dtype) or labels.min() < 0 or labels.max() > labels.size):
       # Use np.unique to generate the label indices. `new_labels` will
       # be 1-d, but it should be interpreted as the flattened n-d array of
       # label indices.
       unique_labels, new_labels = np.unique(labels, return_inverse=True)
       counts = np.bincount(new_labels)
       sums = np.bincount(new_labels, weights=input.ravel())
       if centered:
          # Compute the sum of the mean-centered squares.
          # We must reshape new_labels to the n-d shape of `input` before
          # passing it _sum_centered.
          sums_c = _sum_centered(new_labels.reshape(labels.shape))
       idxs = np.searchsorted(unique_labels, index)
       # make all of idxs valid
       idxs[idxs >= unique_labels.size] = 0
       found = (unique_labels[idxs] == index)
    else:
       # labels are an integer type allowed by bincount, and there aren't too
       # many, so call bincount directly.
       counts = np.bincount(labels.ravel())
       sums = np.bincount(labels.ravel(), weights=input.ravel())
       if centered:
          sums_c = _sum_centered(labels)
       # make sure all index values are valid
       idxs = np.asanyarray(index, np.int).copy()
       found = (idxs >= 0) & (idxs < counts.size)
       idxs[~found] = 0


    counts = counts[idxs]
    counts[~found] = 0
    sums = sums[idxs]
    sums[~found] = 0

    if not centered:
       return (counts, sums)
    else:
       sums_c = sums_c[idxs]
       sums_c[~found] = 0
       return (counts, sums, sums_c) 
     
     
# ###############################################################################     
#  FUNCTION USED BY STATS
# ###############################################################################
     
def _safely_castable_to_int(dt): 
  
   '''
   Test whether the numpy data type `dt` can be safely cast to an int.
   '''
   
   int_size = np.dtype(int).itemsize
   safe = ((np.issubdtype(dt, int) and dt.itemsize <= int_size) or 
           (np.issubdtype(dt, np.unsignedinteger) and dt.itemsize < int_size))
   
   return safe
        
    
