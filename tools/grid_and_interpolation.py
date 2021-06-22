#!/usr/bin/env python


# Copyright (c) 2012- 2016

# TROPOS,
# Permoserstr. 15
# 04318 Leipzig, Germany. 

# Author:
# ====== 
# Fabian Senf <senf@tropos.de>


# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as 
# published by the Free Software Foundation; either version 3 of 
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.


import sys, copy

import numpy as np
from scipy.interpolate import griddata

try:
    import SimpleITK
except:
    print 'SimpleITK not available.'


######################################################################
######################################################################
######################################################################
#####                          #######################################
##### GRID MANIPULATIONS, ETC. #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################



def COSMO_mean_vertical_levels(Nlev):

    ''' 
    Outputs the mean vertical levels of COSMO-DE taken from the 
    model description of COSMO (see reference below).


    USAGE:
    =====
    lev = COSMO_mean_vertical_levels(Nlev)


    INPUT:
    =====
    Nlev: number of vertical levels, i.e.
          Nlev = 51 for intermediate half-levels (e.g. for w)
          Nlev = 50 for main levels (e.g. for u)


    OUPUT:
    =====
    lev: mean vertical COSMO-DE levels 

           
    REFERENCE:
    =========
    M. Baldauf, J. Foerstner, S. Klink, T. Reinhardt, C. Schraff, A. Seifert und K. Stephan
    "Kurze Beschreibung des Lokal-Modells Kuerzestfrist COSMO-DE (LMK) und seiner Datenbanken
    auf dem Datenserver des DWD", Version 1.6, Stand: 31. 03. 2011
    
    Tabelle 1, Seite 29


    '''



    zm = np.array([22000.0000,21000.0000,20028.5703,19085.3594,18170.0000,17282.1406,16421.4297,\
                   15587.5000,14780.0000,13998.5703,13242.8594,12512.5000,11807.1367,11126.4297,\
                   10470.0000,9837.5000,9228.5703,8642.8594,8080.0000,7539.6367,7021.4297,\
                   6525.0000,6050.0000,5596.0664,5162.8594,4750.0000,4357.1367,3983.9299,\
                   3630.0000,3295.0000,2978.5701,2680.3601,2400.0000,2137.1399,1891.4299,\
                   1662.5000,1450.0000,1253.5698,1072.8599,907.5000,757.1399,621.4299,\
                   500.0000,392.5000,298.5698,217.8600,150.0000,94.6400,51.4300,20.0000,0.0000])

    if Nlev == 51:
        return zm / 1e3
    elif Nlev == 50:
        return 0.5e-3 * (zm[1:] + zm[:-1])


######################################################################
######################################################################

def ll2xy(lon, lat, lon0=10., lat0=50., R = 6380.):

    ''' 
    Transformation between longitude and latitude and local Cartesian coordinates
    in west-east direction (x) and south-north direction (y).


    USAGE:
    =====
    x,y =  ll2xy(lon, lat, lon0=10., lat0=50., R = 6380.)
    

    INPUT:
    =====
    lon: longitude
    lat: latitude
    lon0: arbitrary initial longitude for x = 0 (OPTIONAL)
    lat0: arbitrary initial latitude for y = 0 (OPTIONAL)
    R: approximation for the radius of Earth (OPTIONAL)
    

    OUPUT:
    =====
    x: Cartesian coordinate in west-east direction (increasing towards the East)
    y: Cartesian coordinate in south-north direction (increasing towards the North)

    
    COMMENT:
    ========
    The unit of Earth's radius determines the units of x and y. 
    Default (km).

    '''


    # arbritary initial value for longitude and latitude
    lam0 = np.deg2rad(lon0)
    phi0 = np.deg2rad(lat0)


    # conversion from degree to radiant
    lam = np.deg2rad(lon)
    phi = np.deg2rad(lat)
    
    # get x and y from lam and phi
    x = R  * np.cos(phi) * (lam - lam0)
    y = R * (phi - phi0)

    return x,y 

######################################################################
######################################################################

def xy2ll(x, y, lon0=10., lat0=50., R = 6380):
    ''' 
    Transformation between local Cartesian coordinates in west-east direction (x) 
    and south-north direction (y) and longitude and latitude.


    USAGE:
    =====
    lon, lat = xy2ll(x, y, lon0=10., lat0=50., R = 6380)
    

    INPUT:
    =====
    x: Cartesian coordinate in west-east direction (increasing towards the East)
    y: Cartesian coordinate in south-north direction (increasing towards the North)
    lon0: arbitrary initial longitude for x = 0 (OPTIONAL)
    lat0: arbitrary initial latitude for y = 0 (OPTIONAL)
    R: approximation for the radius of Earth (OPTIONAL)
    

    OUPUT:
    =====
    lon: longitude
    lat: latitude

    
    COMMENT:
    ========
    The unit of Earth's radius determines the units of x and y. 
    Default (km).

    '''


    # arbritary initial value for longitude and latitude
    lam0 = np.deg2rad(lon0)
    phi0 = np.deg2rad(lat0)

    # get phi and lam from x and y 
    phi = y / R + phi0
    lam = x / ( R * np.cos(phi) ) + lam0
    
    # conversion from radiant to degree
    lon = np.rad2deg(lam)
    lat = np.rad2deg(phi)
    

    return lon, lat

######################################################################
######################################################################
######################################################################
#####                          #######################################
##### HRV Upscaling            #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################


def make_hrv_upscaling(v):
    
    '''
    A 2d field is upscaled to 3-fold higher resolution using linear 
    interpolation. Edge values are determined by linear extrapolation.
   

    USAGE:
    =====
    vhigh = make_hrv_upscaling(v)

    
    INPUT
    =====
    v: 2d variable field (Nrow, Ncols)

    
    OUTPUT
    ======
    vhigh: 2d variable field with 3-times higher resolution (3*Nrows, 3*Ncols)
    
    '''

    var = make_add_edge(v)

    # x shift first ..................................................
    dvar_x = var[1:,:] - var[:-1,:]
    var0 = var[:-1,:]
    var1 = var0 + 1./3. * dvar_x
    var2 = var0 + 2./3. * dvar_x

    Nx,Ny = var0.shape

    varx = np.zeros((3*Nx,Ny))
    
    varx[0::3,:] = var0
    varx[1::3,:] = var1
    varx[2::3,:] = var2
    
    dvar_y = varx[:,1:] - varx[:,:-1]
    var0 = varx[:,:-1]
    var1 = var0 + 1./3. * dvar_y
    var2 = var0 + 2./3. * dvar_y

    Nx,Ny = var0.shape

    vary = np.zeros((Nx,3*Ny))
    
    vary[:,0::3] = var0
    vary[:,1::3] = var1
    vary[:,2::3] = var2

    vhigh =  vary[2:-1,2:-1]
    

    return vhigh



######################################################################
######################################################################


def make_add_edge(var):

    '''
    Adds egde region to 2d data field. Values in added left and right column, 
    and lower and upper row are linearily extrapolated from their neighbors.
   


    USAGE:
    =====
    var_new = make_add_edge(var)

    
    INPUT
    =====
    var: 2d variable field (Nrow, Ncols)

    
    OUTPUT
    ======
    var_new: 2d variable field with added edge values (Nrows + 2, Ncols + 2)
    
    '''
    
    N1, N2 = var.shape
    
    var_new = np.zeros((N1+2,N2+2))
    
    var_new[1:N1+1,1:N2+1] = var
    
    # upper edge
    dv = var[1,:] - var[0,:]
    var_new[0,1:N2+1] = var_new[1,1:N2+1] - dv

    # lower edge
    dv = var[-2,:] - var[-1,:]
    var_new[-1,1:N2+1] = var_new[-2,1:N2+1] - dv

    # left edge
    dv = var_new[:,2] - var_new[:,1]
    var_new[:,0] = var_new[:,1] - dv

    # right edge
    dv = var_new[:,-3] - var_new[:,-2]
    var_new[:,-1] = var_new[:,-2] - dv

    

    return var_new



######################################################################
######################################################################
######################################################################
#####                          #######################################
##### INTERPOLATION STUFF, ETC #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################




def get_index(p, lon, lat):
    
    '''
    Given a point p the nearest grid point is returned.

    USAGE:
    =====
    ir, ic = get_index(p, lon, lat)
    

    INPUT:
    =====
    p: point in the 2d grid
    lon: grid values of the 1st dimension (can be longitude, x, etc.)
    lat: grid values of the 2nd dimension (can be latitude, y, etc.)
    

    OUPUT:
    =====
    ir: row index (the 1st)
    ic: column index (the 2nd)
    '''

    # get shapes of fields 
    Nr, Nc = lon.shape
    
    # start couting the number of rows and columns
    jr = np.arange(Nr)
    jc = np.arange(Nc)
    
    # make a mesh of row and column counts
    ic, ir = np.meshgrid(jc, jr)

    # transfrom the point in a shape used in the interpolation routine
    p = np.asarray([p]).T
    
    # do the interpolation in the row and column mesh
    irp = interpolate_field_to_hor_pos(p, lon, lat, ir, method='nearest')
    icp = interpolate_field_to_hor_pos(p, lon, lat, ic, method='nearest')
    
    return irp, icp

######################################################################
######################################################################


def interpolate_field_to_hor_pos(pos, lon, lat, field, **kwargs):

    '''
    Given a sequence of points (or just one) a field (2d or 3d) given a
    grid is interpolated to the new positions.


    USAGE:
    =====
    f = interpolate_field_to_hor_pos(pos, lon, lat, field, **kwargs)
    

    INPUT:
    =====
    pos: sequence of points given as np.array((x,y))
         where x and y can be n-dimensional 
    lon: grid values of the 1st dimension (can be longitude, x, etc.)
    lat: grid values of the 2nd dimension (can be latitude, y, etc.)
    field: 2d or 3d field to be interpolated to pos

    **kwargs: keyword arguments for the griddata routine of
              the scipy.interpolate package can be passed through  

    OUPUT:
    =====
    f: field interpolated to the position sequence pos,
       f has the dimension n x Nlev, 
       where n is the dimension of pos, Nlev the vertical dimension

    '''

      
    # set dimensions .................................................
    Nhor = lon.flatten().shape[0]

    if len(field.shape) == 3:
        Nlev = field.shape[2]
    else:
        Nlev = 1


    # set input coordinates ..........................................
    Pin = np.zeros((Nhor,2))
    Pin[:,0] = lon.flatten()
    Pin[:,1] = lat.flatten()
    
    # set input field ................................................
    fin = field.flatten().reshape(Nhor, Nlev)
    
    
    # output coordinates
    NP = pos.flatten().shape[0]/2
    Pout = pos.T
    
    # horizontal interpolation ...........................................
#    print 'Do horizontal interpolation'
    if not kwargs.has_key('method'):
        kwargs['method'] = 'nearest'
    
    fout = griddata(Pin,fin,Pout, **kwargs)
    
    if NP == 1:
        profile = fout.squeeze()
    else:
        profile = fout
    
    return profile

######################################################################
######################################################################



def make_index_set(nrows, ncols):

    '''
    This simple routine makes an 2d index set for given row and column
    number of a 2d field.

    INPUT
    =====
    nrows: number of rows
    ncols: number of columns
    
    OUTPUT
    ======
    ii: 2d field of rows indices
    jj: 2d field of column indices
    '''
    
    i = np.arange(0, nrows)
    j = np.arange(0, ncols)
    ii, jj = np.meshgrid(i, j, indexing = 'ij')

    return ii, jj

######################################################################
######################################################################



def create_interpolation_index(lon1, lat1, lon2, lat2):

    '''
    Given georeference of two grids in lon and lat, an index is built to map 
    a field on grid 1 to grid 2.


    USAGE:
    =====
    ir, ic = create_interpolation_index(lon1, lat1, lon2, lat2)
    

    INPUT:
    =====
    lon1:  longitude of input grid 1
    lat1:  latitude of input grid 1
    lon2:  longitude of target grid 2
    lat2:  latitude of target grid 2


    OUPUT:
    =====
    ir:  index field for transforming rows
    ic:  index field for transforming columns
    '''


    Nr, Nc = lon1.shape
    s2 = lon2.shape

    ir = np.arange(Nr)
    ic = np.arange(Nc)

    ic1, ir1 = np.meshgrid(ic, ir)


    x1, y1 = ll2xy(lon1, lat1)
    x2, y2 = ll2xy(lon2, lat2)

    pos = np.array((x2.flatten(), y2.flatten()))
    
    ir2 = interpolate_field_to_hor_pos(pos, x1, y1, ir1).reshape(s2) 
    ic2 = interpolate_field_to_hor_pos(pos, x1, y1, ic1).reshape(s2)


    return ir2, ic2


######################################################################
######################################################################

def curve_flow_filter(f, numberOfIterations = 5):

    '''
    Smoothing filter depending on isoline curvature. Interface for 
    curvature flow filter from simpleITK toolkit.

    USAGE
    =====
    f_sm = curve_flow_filter(f, numberOfIterations = 5)


    INPUT
    =====
    f: 2d field
    numberOfIterations: (optional) number of iterations, increases smooting effect


    OUTPUT
    ======
    f_sm: smoothed 2d field
    '''

    img = SimpleITK.GetImageFromArray(f)
    img_sm = SimpleITK.CurvatureFlow(img, numberOfIterations = numberOfIterations)
    
    f_sm = SimpleITK.GetArrayFromImage(img_sm)

    return f_sm
    

######################################################################
######################################################################


def remap_field(lon, lat, f, dr = 0.5):

    '''
    Do nearest nearbor remapping on equi-distant
    local cartesian coordiante grid.
    
    Uses gradient flow filter for pixel edge
    correction.

    USAGE
    =====
    xnew, ynew, fint = remap_field(lon, lat, f, dr = 0.5)


    INPUT
    =====
    lon: longitude
    lat: latitutde
    f: 2d field
    dr: (optional) grid spacing in km


    OUTPUT
    ======
    xnew: regridded, equi-distant x-coordinate in km
    ynew: regridded, equi-distant y-coordinate in km
    fint: remapped field
    
    '''

    # convert to local coordiantes
    x, y = ll2xy(lon, lat)


    # make new grid
    xn = np.arange(x.min(), x.max(), dr)
    yn = np.arange(y.min(), y.max(), dr)
    xnew, ynew = np.meshgrid(xn, yn)

    nlon, nlat = xy2ll(xnew, ynew)
    ir, ic = create_interpolation_index(lon, lat, nlon, nlat)

    fint = curve_flow_filter(f[ir, ic], 3)

    return xnew, ynew, fint

######################################################################
######################################################################


def make_vert_cut(p1, p2, lon, lat, vg, **kwargs):


    '''
    Makes a horizontal - vertical cut through a 3d field.


    USAGE:
    =====
    s, lo, la, v = make_vert_cut(p1, p2, lon, lat, vg, **kwargs)
    

    INPUT:
    =====
    p1: 1st end of the cutting line
    p2: 2nd end of the cutting line
    lon: longitude
    lat: longitude
    vg: 3d field to be interpolated

    **kwargs: keyword arguments for the griddata routine of
              the scipy.interpolate package can be passed through  

    OUPUT:
    =====
    s: distance of km from one to the other end of the cutting line
    lo: longitude along the cut
    la: latitude along the cut
    v: field interpolated to the cutting line

    '''


    # prepare grid and start/end points ------------------------------
    # get grid in x and y
    xg, yg = ll2xy(lon,lat)

    # get points
    x1, y1 = ll2xy(p1[0],p1[1])
    x2, y2 = ll2xy(p2[0],p2[1])


    # get distance between start and end point
    d = abs((x1 - x2) + 1j * (y1 -y2))
    # ================================================================
    

    # prepare subsampling along the cutout ---------------------------

    # set smallest distance between the new grid points
    ds = 2.


    # number of new along-cutout grid points 
    N = int(d / ds)


    # subsample along the cutout
    x = np.linspace(x1,x2,N,endpoint=True)
    y = np.linspace(y1,y2,N,endpoint=True)

    pos = np.array((x,y))

    # get the along-cut distance
    x0 = x[0]
    y0 = y[0]

    lo, la = xy2ll(x, y)
    s = abs((x - x0) + 1j * (y - y0))
    # ================================================================


    
    # do interpolation -----------------------------------------------
    if not kwargs.has_key('method'):
        kwargs['method'] = 'nearest'
    v = interpolate_field_to_hor_pos(pos, xg, yg, vg, **kwargs)
    # ================================================================
    

    
    return s, lo, la, v

######################################################################
######################################################################


def simple_pixel_area(lon,lat, uncertainty = False):

    ''' 
    Approximates the grid box area of a lon - lat grid.
    
    The grid is treated as if it is an equidistant grid in local 
    Cartesian coordinates on Earth's surface (i.e. in x and y). 
    The grid point is placed at each of the for corners of the grid box
    and the corresponding average area ist calculated.

    USAGE
    =====
    a =  pixel_area(lon,lat)

    
    INPUT
    =====
    lon: longitude
    lat: latitutde

    
    OUTPUT
    ======
    a: grid box area
    
    '''

    if not np.rank(lon) == 2 and not np.rank(lat) == 2:
        print 'ERROR: please only use 2d matrices for lon and lat'
        sys.exit(0)

    # convert lon,lat into Cartesian coordinates
    x,y = ll2xy(lon,lat)

    # add some edge via linear extrapolation
    x_ext = make_add_edge(x)
    y_ext = make_add_edge(y)

    # make forward, backward difference assuming that 
    #  *) x is connected to column
    dx1 = x_ext[1:-1, 2:] - x 
    dx2 = x - x_ext[1:-1,:-2]

    #  *) y is connected to row
    dy1 = y - y_ext[2:, 1:-1]  
    dy2 = y_ext[:-2, 1:-1] - y
    
    # calculate the area assuming that the grid point is located at 
    # each of the four corners of the grid box
    ac = []
    ac.append(dx1 * dy1)
    ac.append(dx1 * dy2)
    ac.append(dx2 * dy1)
    ac.append(dx2 * dy2)

    # calculate average area
    am = np.dstack(ac).mean(axis=2)
    

    if uncertainty:
        # calculate deviation of different area 
        # (perhaps a measure of uncertainty that simple method)
        da = np.dstack(ac).std(axis=2)
        return am ,da

    else:
        return am


######################################################################
######################################################################


def i2iset(v,i):

    '''
    Returns an index set for a n-dim numpy array given an index i
    which points to a position within the same but flattened array.
    
    It is assumed that the array was flattened in C-mode.

    USAGE
    =====
    iset = i2iset(v,i)


    INPUT
    =====
    v: field which determines the rank and shape
    i: input index which points to the position in v.flatten()

    
    OUTPUT
    ======
    iset: index set which points to the position in v

    '''
   

    # get shape of input n-dim array
    s = v.shape
   
    # calculate factors which determine the jumps within 
    # the flattened array
    factors = np.cumprod(s[::-1])[:-1]
 
    # loop over different dimension of the n-dim array
    res = i
    iset = []
    for f in factors[::-1]:
        iset.append(res / f)
        res = i % f

    iset.append(res)

    return tuple(iset)



######################################################################
######################################################################


def ldiff(var, axis = -1):
    '''
    Calculates difference of a field at adjacent levels. Now, wrapper
    for numpy.diff with default on last axis.


    USAGE:
    =====
    dvar = ldiff(var, axis = -1.)

    
    INPUT
    =====
    var: field; 1d, 2d or 3d array

    
    OUTPUT
    ======
    dvar: level difference
    
    '''

 
    return np.diff(var, axis = axis)
 

######################################################################
######################################################################


def lmean(var, axis = -1):

    '''
    Calculates the layer mean of a field.

    USAGE:
    =====
    varm = lmean(var)

    
    INPUT
    =====
    var: field; 1d, 2d or 3d array

    
    OUTPUT
    ======
    varm: layer mean
    
    '''

    # select the axis for layer-average calculations
    if axis == -1:
        ndim = len(var.shape)
    else:
        ndim = axis + 1

    if ndim == 1:
        dvar = 0.5 * (var[1:] + var[:-1])

    elif ndim == 2:
        dvar = 0.5 * (var[:,1:] + var[:,:-1])

    elif ndim == 3:
        dvar = 0.5 * (var[:,:,1:] + var[:,:,:-1])

    return dvar
 


######################################################################
######################################################################



def cutout_fields(fin , slices, vaxis = 0):

    '''
    Cuts out a field or a list of similar fields given row and column slices.


    USAGE:
    =====
    fout = cutout_fields(fin , slices)

    
    INPUT
    =====
    fin: 2d or 3d field or list of fields - last two dimension are consider
    slice: tuple such that ((row1, row2), (col1, col2))

    vaxis: for 3d fields, axis which is still varying (optional, default = 0)
    

    
    OUTPUT
    ======
    fout: resulting cutout (possibly a list of cutted fields)
    
    '''

    (ir1, ir2), (ic1, ic2) = slices

    fnew = []

    if not type(fin) == type([]):
        flist = [fin]
    else:
        flist = fin

    for f in flist:
        if np.ndim(f) == 2:
            fcut = f[ir1:ir2, ic1:ic2]
        elif np.ndim(f) == 3:
            if vaxis == 0:
                fcut = f[:, ir1:ir2, ic1:ic2]

            elif vaxis == 2:
                fcut = f[ir1:ir2, ic1:ic2, :]


        fnew.append(fcut)


    if not type(fin) == type([]):
        return fnew[0]
    else:
        return fnew


######################################################################
######################################################################

def cutout_field4box(f, ind, bsize, **kwargs):

    '''
    Cuts out a field based on center pix index and box size
    '''

    # extract the vaxis keyword if there
    vaxis = kwargs.get('vaxis', 0)

    # get the fill value
    fill_value = kwargs.pop('fill_value', 0)


    # get field dimensions and prepare stacking if needed for adding edges

    if np.ndim(f) == 2:
        nrows, ncols = f.shape

    elif np.ndim(f) == 3:
        if vaxis == 0:
            ntime, nrows, ncols = f.shape

        if vaxis == 2:
            nrows, ncols, nlev = f.shape

    # get index values
    irow, icol = ind


    # set box edges
    ir1 = irow -  bsize / 2
    ir2 = irow +  bsize / 2  +  1

    ic1 = icol -  bsize / 2
    ic2 = icol +  bsize / 2  +  1

    # check for box which extend the field dimensions
    dc_left, dc_right = (0, 0)
    dr_left, dr_right = (0, 0)

    if ic1 < 0:
        dc_left =  - ic1
        ic1 = 0

    if ic2 > ncols:
        dc_right = ic2 - ncols
        ic2 = ncols  

    if ir1 < 0:
        dr_left =  - ir1
        ir1 = 0

    if ir2 > nrows:
        dr_right = ir2 - nrows
        ir2 = nrows 
   
    #set region 
    reg = ((ir1, ir2), (ic1, ic2))

    # do the cut
    fcut =  cutout_fields(f, reg, **kwargs)

    if np.ndim(f) == 2:
        fcut = np.ma.reshape(fcut, fcut.shape + (1,))

    elif np.ndim(f) == 3 and vaxis == 0:
        fcut = fcut.transpose(1,2,0)


    # and now add edges if needed
    s = fcut.shape
    if dc_left != 0:
        z = np.ma.zeros( (s[0], dc_left, s[2] ) )
        fcut = np.column_stack([z, fcut])

    if dc_right != 0:
        z = np.ma.zeros( (s[0], dc_right, s[2]) )
        fcut = np.column_stack([fcut, z])

    s = fcut.shape
    if dr_left != 0:
        z = np.ma.zeros( (dr_left,) + s[1:] )
        fcut = np.row_stack([z, fcut])

    if dr_right != 0:
        z = np.ma.zeros( (dr_right,) + s[1:] )
        fcut = np.row_stack([fcut, z])

    if np.ndim(f) == 3 and vaxis == 0:
        fcut = fcut.transpose(2,0,1)


    return fcut.squeeze().astype(f.dtype)

    
####################################################################
####################################################################


def tube_cutout4box(v3d, i0, boxsize):

    '''
    Performs tube cutout with given index set.

    INPUT
    =====
    v3d: 3d field, 1st dimension is time
    i0: (irvec, icvec) index tuple, (row index vector, co index vector)
    boxsize: boxsize

    OUTPUT
    ======
    vtube: 3d tube
    '''

    # get number of time steps
    ntime = v3d.shape[0]

    # extract index vector
    irvec, icvec = i0

    # initialize tube field
    vtube = np.ma.zeros((ntime, boxsize, boxsize))
    
    # and loop over time for cutout
    for n in range(ntime):
 
        # row and column index
        ir = irvec.flatten()[n]
        ic = icvec.flatten()[n]

        # cutout box
        vtube[n] =  cutout_field4box(v3d[n], (ir, ic), boxsize)

    return vtube


######################################################################
######################################################################


def cutout_cluster(c, nc, 
                   nedge = 30):
    
    '''
    Makes a cutout of a categorial field c for an object of class / number nc.


    INPUT
    =====
    c: 2d categorial field
    nc: category number / class which is chosen
    nedge: optional, number of edge pixels included


    OUTPUT
    =====
    ccut: cutout of categorial field c
    '''


    # get index field ------------------------------------------------
    nrow, ncol = c.shape
    irow, icol = make_index_set(nrow, ncol)
    # ================================================================


    # make mask and mask the indices ---------------------------------
    mask = (c== nc)

    ir = irow[mask]
    ic = icol[mask]
    # ================================================================

    
    # determine cutout region ----------------------------------------
    ir1 = ir.min() - nedge
    ir2 = ir.max() + nedge
    
    ic1 = ic.min() - nedge
    ic2 = ic.max() + nedge

    # check bounds
    if ir1 < 0:
        ir1 = 0

    if ir2 >= nrow:
        ir2 = nrow 

    # check bounds
    if ic1 < 0:
        ic1 = 0

    if ic2 >= ncol:
        ic2 = ncol 

    region = ((ir1, ir2), (ic1, ic2))
    # ================================================================


    # do the cutout --------------------------------------------------
    ccut = cutout_fields(c, region)
    # ================================================================

    return ccut

######################################################################
######################################################################

def region2slice(lon, lat, region):

    '''
    Given a region of ((lon1,lon2),(lat1, lat2)) the routine provides 
    a slicing ((row1, row2), (col1, col2)) than accomodates the region.


    USAGE:
    =====
    slice = region2slice(lon, lat, region)

    
    INPUT
    =====
    lon: longitude
    lat: latitude
    region:  ((lon1,lon2),(lat1, lat2)) 
    

    
    OUTPUT
    ======
    slice = ((row1, row2), (col1, col2))
    
    '''


    # get corner points
    (lon1, lon2), (lat1, lat2) = region

    # get the correspoding index of the corners 
    ir1, ic1 = get_index((lon1, lat2), lon, lat)
    ir2, ic2 = get_index((lon2, lat1), lon, lat)

    return ((ir1, ir2), (ic1, ic2))


######################################################################
######################################################################



def low2hres_region(region):

    '''
    Based on low-res. region of MSG channels an correspoding HRV region is 
    calculated.


    USAGE:
    =====
    hrv_region = low2hres_region(region)

    
    INPUT
    =====
    region : region as  ((row1, row2), (col1, col2))
    

    
    OUTPUT
    ======
    hrv_region: corresponding HRV region as ((row1, row2), (col1, col2))
    
    '''

    # calculate HRV region
    # ====================
    # It is assumed that the hrv region is a dependent variable which 
    # depends on the low res region attribute
    # That means: Given the low res region cutout hrv region is determined 
    # as corresponding cutout which 
    # (i) exactly fits with the low res region and 
    # (ii) refers to the artifical hrv full disk of 11136 x 11136
    #
    # low res pixel with the index (I,J) = (0,0)
    # has a high res pixel with index (i_m, j_m) = (2,2) in the middle
    # and (i_u, j_u) = (1,1) in the upper corner
    # see doc "MSG Level 1.5 Image Data Format Description", Fig.8
    # 
    #  which then leads to the relation between the low res pixel (I,J)
    # and its corresponding upper corner high res pixel (i_u, j_u):
    #  (i_u, j_u) = 3 * (I, J ) + 1
        
        
    (r1, r2), (c1, c2) = region 
    
    hrv_region = ((3*r1 + 1, 3*r2 + 1), (3*c1 + 1, 3*c2 + 1))
    
    return hrv_region
    

######################################################################
######################################################################


if __name__ == '__main__':
    # test cutout_field4box
 
    test = np.arange(25).reshape(( 5, 5)).astype(np.int)

    c = cutout_field4box(test, (2,2), 15)
    print c
    print c.shape


    test = np.arange(25*2).reshape((5, 5, 2)).astype(np.int)

    c = cutout_field4box(test, (2,2), 7, vaxis = 2)
    print c[:,:,0]
    print
    print c[:,:,1]
    print c.shape
