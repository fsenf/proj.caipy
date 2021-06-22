#!/usr/bin/env python

# Load libraries ---------------------------------------------------------
import glob
import xarray as xr
import sys
import netCDF4 as nc
import numpy as np
import pandas as pd
from netCDF4 import num2date

##########################################################################
##########################################################################

#def read_icon_data_3d_coarse_slice(fname, var_list, itime = 0, ilev=0):
#
#    '''
#    Reads in one slice of ICON 3d_coarse data for one particular time-step
#
#    USAGE
#    ====
#    data = read_icon_data_3d_coarse(file,[var1, var2],nt,nlev)
#
#    INPUT
#    =====
#    file: path and filename of the data
#    var_list: list of variable names
#    itime: time-step to be read in, default=first time-step
#
#    OUTPUT
#    ======
#    dset: python dictionary, single variables are accessible via
#          dset['var1']
#          example: dset['var1'].shape would give (nlat,nlon) 
#    '''  
#
#    dset = {}
#
#    # input dimensions
#    f = nc.Dataset(fname, 'r')
#
#    # input data
#    for vname in var_list:
#        v = f.variables[vname]
#        dset[vname] = np.array( v[itime,ilev,:,:] )
#
#    f.close()
#
#    return dset

######################################################################
######################################################################

#def read_icon_3d_coarse_grid(fname, var_list):

#    '''
#    Reads in ICON 3d_coarse grid information (lat,lon or height)
#
#    USAGE
#    ====
#    data = read_icon_data_3d_coarse_hor_grid(file,[var1, var2])
#
#    INPUT
#    =====
#    file: path and filename of the data
#    var_list: list of variable names
#
#    OUTPUT
#    ======
#    dset: python dictionary, single variables are accessible via
#          dset['var1']
#          example: dset['lon'].shape would give (nlon) 

#    CAVEAT
#    ======
#    function is flexible to read in horizontal and vertical grid information,
#    but both types cannot be stored in one dictionary as there would be a 
#    mismatch of dimensions
#    '''  
#
#    dset = {}
#
#    # input dimensions
#    f = nc.Dataset(fname, 'r')
#
#    # input data
#    for vname in var_list:
#        v = f.variables[vname]
#        dset[vname] = np.array( v[:] )
#
#    f.close()
#
#    return dset

######################################################################
######################################################################

def convert_frac_time_to_readable_time(time_frac):

    '''
    Converts ICON's time unit day as %Y%m%d.%f to readable times

    USAGE
    ====
    dates, frac_day, datetime_obj= convert_frac_time_to_readable_time(time_frac)

    INPUT
    =====
    time_frac: ICON's time-array

    OUTPUT
    ======
    dates: Array of time_frac.shape() containing the date as string ('20140704')
    frac_day: Array time_frac.shape() containing the seconds of the day as integer (25200)
    datetime_obj: a pandas date-time-object with time information ('2014-07-04 07:00:00')
    '''  

    # Splitting time in fractional day and dates
    frac_day, dates = np.modf(time_frac)

    # Converting frac_day to seconds
    frac_day = frac_day * 24 * 60 * 60
 
    # Rounding frac_day to nearest second
    frac_day = np.round(frac_day)

    # Convert frac_day to integer
    frac_day = frac_day.astype(int)

    # Convert frac_day to timedelta-object
    frac_day_obj = pd.to_timedelta(frac_day, unit='s')

    # Convert dates to string
    dates = dates.astype(int).astype(str)    

    # Convert dates to datetime-object
    dates_obj = pd.to_datetime(dates)

    # Combining both dates_obj and frac_day_obj
    datetime_obj = dates_obj + frac_day_obj

    return dates, frac_day, datetime_obj

#####################################################################
####################################################################

def convert_gregorian_time_to_frac_time(time, date):

    '''
    Converts YYYYMMDD HH:MM:SS to unit day time as %Y%m%d.%f 

    USAGE
    ====
    frac_time= convert_gregorian_time_to_frac_time(time, date)


    INPUT
    =====
    time: RADOLANs time-array
    date: YYYYMMDD 


    OUTPUT
    ======
    frac_time = time array in the form %Y%m%d.%f

    '''
    

    frac_day = []
    frac_time  = []

    #frac_time= np.empty([len(time)])

    for t in range(len(time)):

       TIME= num2date(time[t], units = 'seconds since 1970-01-01 00:00:00', calendar= 'gregorian')

       # converting hours+minutes+seconds to fractions of day
       frac_day = TIME.hour/(24.) + TIME.minute/(24*60.) + TIME.second/(24*60*60.)


       # time in the form %Y%m%d.%f
       frac_time0 = float(date) + frac_day
       frac_time.append(frac_time0)
   

    return frac_time 


######################################################################
#####################################################################

def convert_gregorian_time_to_readable_time(time, date):

    '''
    Converts RADOLAN time unit day from seconds since 1970-01-01 00:00:00  to readable times

    USAGE
    ====
    dates, frac_day, datetime_obj= convert_gregorian_time_to_readable_time(time, date)

    INPUT
    =====
    time: RADOLANs time-array
    date: YYYYMMDD 

    OUTPUT
    ======
    dates0: Array of time_frac.shape() containing the date as string ('20140704')
    frac_day: Array time_frac.shape() containing the seconds of the day as integer (25200)
    datetime_obj: a pandas date-time-object with time information ('2014-07-04 07:00:00')
    '''

    # convert time to a readable format 

    frac_day = []
    datetime_obj = []    

    for t in range(len(time)):

       TIME= num2date(time[t], units = 'seconds since 1970-01-01 00:00:00', calendar= 'gregorian') 
         
       # Converting hours+seconds to seconds
       frac_day0 = TIME.hour*60*60 + TIME.minute*60

       # Rounding frac_day to nearest second
       frac_day0 = np.round(frac_day0)

       # Convert frac_day to integer
       frac_day0 = frac_day0.astype(int)
       frac_day.append(frac_day0)

       # Convert frac_day to timedelta-object
       frac_day_obj0 = pd.to_timedelta(frac_day0, unit='s')

       # Convert dates to string
       dates0 = str(date)

       # Convert dates to datetime-object
       dates_obj0 = pd.to_datetime(dates0)

       # Combining both dates_obj and frac_day_obj
       datetime_obj0 = dates_obj0 + frac_day_obj0
       datetime_obj.append(datetime_obj0)

    return dates0, frac_day, datetime_obj

######################################################################
######################################################################


def read_field(**kwargs):
    
    '''
    Reads the 2d field for cluster analysis.
    
    INPUT
    =====
    kwargs: dict. of input parameters (might depent on choice)
    
    OUTPUT
    ======
    time: time vector
    lon: longitude (matrix)
    lat: latitude (matrix)
    lonv: longitude (vector)
    latv: latitude (vector)
    n_latlon: total number of valid grid points
    pos_grid: locations with valid data
    z: height above ground of the field
    f: field
    rain: radolan rainfall rates (mm/h)
    '''
    
    var_name = kwargs.get('seg_field', None)

    # --------------------------------------------------------------------------------------
    # read the fields of radolan, radvop, meteosat, synsat and wa 
    # --------------------------------------------------------------------------------------    
    if var_name == 'radolan':
        time, lon, lat, lonv, latv, n_latlon, pos_grid, z, f, rain = read_radolan(**kwargs)

    elif var_name == 'radvop':
        time, lon, lat, lonv, latv, n_latlon, pos_grid, z, f, rain = read_radvop(**kwargs)

    elif var_name == 'synsat':
        time, lon, lat, lonv, latv, n_latlon, pos_grid, z, f = read_synsat(**kwargs)

    elif var_name == 'meteosat':
        time, lon, lat, lonv, latv, n_latlon, pos_grid, z, f = read_meteosat(**kwargs) 

    #elif var_name == 'wa':
    #    time, lon, lat, z, f = read_vertical_velocity(**kwargs)



    if var_name in ['radolan', 'radvop']:
       return time, lon, lat, lonv, latv, n_latlon, pos_grid, z, f, rain
    else:
       return time, lon, lat, lonv, latv, n_latlon, pos_grid, z, f


##########################################################################
##########################################################################
 
def read_radolan(**kwargs):

    '''
    Reads RADOLAN reflectivity factor and rain rate 2d fields.
    
    INPUT
    =====
    kwargs: dict. of input parameters (might depent on choice)
    
    OUTPUT
    ======
    time: time vector
    long_mg: longitude (matrix)
    lat_mg: latitude (matrix)
    lon: longitude (vector)
    lat: latitude (vector)
    n_latlon: total number of valid grid points
    pos_grid: locations with valid data
    z: height above ground of the field
    dbz: radolan reflectivity factor
    rain: radolan rainfall rates (mm/h)
    '''

    # get parameters
    date = kwargs['date']
    ilev = kwargs['ilev']
    if ilev == 0:
       z=0    

    # set input file name
    input_dir = "/work/bm0974/cai_S5/data/radolan"
    input_file = '%s/hdfd_miub_drnet00_l3_dbz_rr_v00_%s-dwd_ICON-grid_628x530.nc' % (input_dir, date)

    print input_dir
    print input_file 

    # Open all files of the day in one dataset (ds)
    ds = nc.Dataset(input_file)


    # Reading in dimensions
    time  = ds['time'][:]    # time as seconds since 1970-01-01 00:00:00
    lat   = ds['lat'][:]     # lat
    lon   = ds['lon'][:]     # lon
    lon_mg, lat_mg = np.meshgrid(lon,lat) # grid


    # Reading in data
    dbz = ds['dbz'][:,:,:]  # RX
    rain= ds['rr'][:,:,:]   # RY

    ds.close()

    

    # #########################################################################
    # RADOLAN mask
    # #########################################################################
    mfile='/radolan_mask_3D.nc'
    mask_file=nc.Dataset(input_dir+mfile)
    mask=mask_file.variables['mask'][:,:,:]

    # mask the fields. Data outside the radolan area is replaced by NAN
    dbz[np.where(mask == 0)]=-999.9
    dbz[np.where(dbz > 92.)]=-999.9
    rain[np.where(dbz == -999.9)]=-999.9
    rain[np.where(rain == -999.0)]=-999.9   

    pos_grid=np.where(mask[0,:,:] ==1)
    n_latlon = len(pos_grid[0])


    return time, lon_mg, lat_mg, lon, lat, n_latlon, pos_grid, z, dbz, rain 



#############################################################################
#############################################################################
def read_rain_gsp_rate(date,dom):

    '''
    Reads rain rate from 2d ICON simulations.
    
    INPUT
    =====
    kwargs: dict. of input parameters (might depent on choice)
    
    OUTPUT
    ======
    rr: ICON 2d rain rate simulation
    '''

    # get parameters
    #date = kwargs['date']
    #dom  = kwargs['dom']
    #ilev = kwargs['ilev']
    #region= kwargs['region_mode']
    #identifier= kwargs['file_identifier']

    #print region
    #print identifier


    #if ilev == 0:
    #   z=0

    # set input file name
    input_dir = "/work/bm0974/cai_S5/tools/ieda/extract_rain_gsp_rate/output"
    input_file = '%s/2d_rain_gsp_rate_%s_ML_%s_genycon_ll.nc' % (input_dir,dom,date)


    #print input_dir
    print input_file

    # Open all files of the day in one dataset (ds)
    ds = nc.Dataset(input_file)

    # Reading in dimensions
    #time  = ds['time'][:]    # time as seconds since 1970-01-01 00:00:00
    #lat   = ds['lat'][:]     # lat
    #lon   = ds['lon'][:]     # lon
    #lon_mg, lat_mg = np.meshgrid(lon,lat) # grid


   # Reading in data
    rr= ds['rain_gsp_rate'][:,:,:]


    ds.close()


    return rr





###############################################################################
##############################################################################  

def read_radvop(**kwargs):

    '''
    Reads synthetic radar reflectivity factor 2d fields.
    
    INPUT
    =====
    kwargs: dict. of input parameters (might depent on choice)
    
    OUTPUT
    ======
    time: time vector
    long_mg: longitude (matrix)
    lat_mg: latitude (matrix)
    lon: longitude (vector)
    lat: latitude (vector)
    n_latlon: total number of valid grid points
    pos_grid: locations with valid data
    z: height above ground of the field
    dbz: radolan reflectivity factor
    '''

    # get parameters
    date = kwargs['date']
    dom  = kwargs['dom']
    ilev = kwargs['ilev']
 
    if ilev == 0:
       z=0

 
    # set input file name
    input_dir = "/work/bm0974/cai_S5/data/radvop"
    input_file = '%s/dbz-mie_%s_%s_ICON-grid_628x530.nc' % (input_dir, dom, date)


    #print input_dir
    print input_file

    # Open all files of the day in one dataset (ds)
    ds = nc.Dataset(input_file)

    # Reading in dimensions
    time  = ds['time'][:]    # time as seconds since 1970-01-01 00:00:00
    lat   = ds['lat'][:]     # lat
    lon   = ds['lon'][:]     # lon
    lon_mg, lat_mg = np.meshgrid(lon,lat) # grid


   # Reading in data
    dbz = ds['dbz'][:,:,:]
    ds.close()
    rain= read_rain_gsp_rate(date,dom)



    # #########################################################################
    # RADOLAN mask
    # #########################################################################
    mfile='/radolan_mask_3D.nc'
    mdir ='/work/bm0974/cai_S5/data/radolan'
    mask_file=nc.Dataset(mdir+mfile)
    mask=mask_file.variables['mask'][0:48,:,:] # radvop daily files contain 30 minutes data

    # mask the fields. Data outside the radolan area is replaced by NAN
    dbz[np.where(mask == 0)]=-999.9
    dbz[np.where(dbz == -99.99)]=-999.9
    rain[np.where(dbz == -999.9)]=-999.9
    rain[np.where(rain == -999.0)]=-999.9


    pos_grid=np.where(mask[0,:,:] ==1)
    n_latlon = len(pos_grid[0])


    return time, lon_mg, lat_mg, lon, lat, n_latlon, pos_grid, z, dbz, rain


##################################################################################
##################################################################################


def read_meteosat(**kwargs):

    '''
    Reads meteosat brightness temperature 2d fields.
    
    INPUT
    =====
    kwargs: dict. of input parameters (might depent on choice)
    
    OUTPUT
    ======
    time: time vector
    long_mg: longitude (matrix)
    lat_mg: latitude (matrix)
    lon: longitude (vector)
    lat: latitude (vector)
    n_latlon: total number of valid grid points
    pos_grid: locations with valid data
    z: height above ground of the field
    bt108: meteosat brightness temperature of the 10.8 micron channel
    '''

    # get parameters
    date = kwargs['date']
    
    # FS comment: Level makes no sense for Meteosat
    # ilev = kwargs['ilev']
    # if ilev == 0:
    z = 0

    # set input file name
    input_dir = "/work/bm0974/cai_S5/data/meteosat"
    input_file = '%s/msevi-bt108-icon_de_dom-%s.nc' % (input_dir, date)

    print input_dir
    print input_file

    # Open all files of the day in one dataset (ds)
    ds = nc.Dataset(input_file)


    # Reading in dimensions
    time  = ds['time'][:]       # day as %Y%m%d.%f
    lat   = ds['lat'][:, 0]     # lat[nlat]
    lon   = ds['lon'][0, :]     # lon[nlon]
    lon_mg, lat_mg = np.meshgrid(lon,lat) # grid


    # Reading in data
    bt108 = -ds['bt108'][:] # negative of Bt108

    ds.close()

    # #########################################################################
    # RADOLAN mask
    # #########################################################################
    mfile='/radolan_mask_3D.nc'
    mdir= '/work/bm0974/cai_S5/data/radolan'
    mask_file=nc.Dataset(mdir+mfile)
    mask=mask_file.variables['mask'][0:48,:,:] # meteosat daily files contain 30 minutes data

    # mask the fields. Data outside the radolan area is replaced by NAN
    bt108[np.where(mask == 0)]=-999.9
    bt108[np.where(bt108 == 0.0)]=-999.9  # Meteosat FillValue = 0.0 
    
    pos_grid=np.where(mask[0,:,:] ==1)
    n_latlon = len(pos_grid[0])


    return time, lon_mg, lat_mg, lon, lat, n_latlon, pos_grid, z, bt108

################################################################################
###############################################################################

def read_synsat(**kwargs):

    '''
    Reads synthetic brightness temperature 2d fields.
    
    INPUT
    =====
    kwargs: dict. of input parameters (might depent on choice)
    
    OUTPUT
    ======
    time: time vector
    long_mg: longitude (matrix)
    lat_mg: latitude (matrix)
    lon: longitude (vector)
    lat: latitude (vector)
    n_latlon: total number of valid grid points
    pos_grid: locations with valid data
    z: height above ground of the field
    bt108: synthetic brightness temperature 
    '''

    # get parameters
    date = kwargs['date']
    dom  = kwargs['dom']
    identifier = kwargs.get('file_identifier', '')
    region_mode = kwargs.get('region_mode', 'default')


    expname = '%s%s' % (region_mode, identifier)

    # FS comment: Level makes no sense for Synsat
    #ilev = kwargs['ilev']
    #if ilev == 0:
    z = 0

    # set input file name
    input_dir = "/work/bm0974/cai_S5/data/synsat"
    input_file = '%s/synsat?_3d_coarse_icon-lem-tstack_%s_%s-%s.nc' % (input_dir, dom, date, expname)


    
    print input_dir
    print input_file

    input_file = glob.glob( input_file )[0]

    # Open all files of the day in one dataset (ds)
    ds = nc.Dataset(input_file)


    # Reading in dimensions
    time_full  = ds['time'][:]       # day as %Y%m%d.%f
    frac_day_full = time_full - np.int( time_full[0] )

    lat   = ds['lat'][:, 0]     # lat[nlat]
    lon   = ds['lon'][0, :]     # lon[nlon]
    lon_mg, lat_mg = np.meshgrid(lon,lat) # grid


    # now time is much finer than half-hourly - we select only the target times
    target_frac_day = np.arange(48.) / 48.

    time_index = []
    for tfrac in target_frac_day:

        # time difference between actual time vector and target time
        dt = np.abs( frac_day_full - tfrac )

        # index of min difference
        ind = dt.argmin()
        
        time_index.append( ind )        
    time_index = np.array( time_index )

    
    print time_index

    # Reading in data
    time = ds['time'][:][time_index]
    bt108 = -ds['bt108'][:][time_index] # negative of Bt108

    ds.close()

    # #########################################################################
    # RADOLAN mask
    # #########################################################################
    mfile='/radolan_mask_3D.nc'
    mdir= '/work/bm0974/cai_S5/data/radolan'
    mask_file=nc.Dataset(mdir+mfile)
    mask=mask_file.variables['mask'][0:48,:,:] # synsat daily files contain 30 minutes data

    # mask the fields. Data outside the radolan area is replaced by NAN
    bt108[np.where(mask == 0)]=-999.9
    bt108[np.where(bt108 == 0.0)]=-999.9  # Meteosat FillValue = 0.0 

    pos_grid=np.where(mask[0,:,:] ==1)
    n_latlon = len(pos_grid[0])


    return time, lon_mg, lat_mg, lon, lat, n_latlon, pos_grid, z, bt108


###########################################################################
###########################################################################

#def read_vertical_velocity(**kwargs):

    # get parameters
#    date = kwargs['date']
#    ilev = kwargs['ilev']
#    dom  = kwargs['dom']

#    region_mode     = kwargs.get('region_mode', 'default')
#    file_identifier = kwargs.get('file_identifier', '')
    
    # Location of ICON output
#    output_dir = "/work/bm0834/k203095/OUTPUT/%s-%s%s/DATA" % (date, region_mode, file_identifier)

    # Location of ICON GRID
#    grid_dir   = "/work/bm0834/k203095/OUTPUT/GRIDS"

    # Grid file
#    grid_file  = "GRID_%s_3d_coarse_ll_%s_ML.nc" % (region_mode, dom)
    
    # Gather all 3d_coarse files for one day
#    pattern_night = "3d_coarse_night_ll_"+dom+"_ML_*"
#    pattern_day   = "3d_coarse_day_ll_"+dom+"_ML_*"
#    night_files   = sorted(glob.glob(output_dir+"/"+pattern_night))
#    day_files     = sorted(glob.glob(output_dir+"/"+pattern_day))
#    files         = night_files + day_files

    # Open all files of the day in one dataset (ds)
#    ds = xr.open_mfdataset(files)

    # Open grid file in standard way
#    grid = nc.Dataset(grid_dir+"/"+grid_file)

    # select only one slice in the vertical
    # for developing purposed read in only some time-steps slice(start,end,increment)
#    ds_sel = ds.isel(height=ilev, height_2=ilev) #time=slice(0,78,1)

    # Reading in dimensions
#    time  = ds_sel['time']            # day as %Y%m%d.%f
#    z_ifc = grid.variables['z_ifc']   # z_ifc[nhalf_levels,nlat,nlon]
#    z_mc  = grid.variables['z_mc']    # z_mc[nfull_levels,nlat,nlon]
#    lat   = grid.variables['lat']     # lat[nlat]
#    lon   = grid.variables['lon']     # lon[nlon]

    # Reading in data
#    wa = ds_sel['wa'] # vertical velocity

    # Convert to numpy array
#    wa = wa.values

    # Mean chosen model level height
#    z_ifc_mean = np.mean(z_ifc[ilev,:,:], (0,1))
#    print("Reading in {} at {} m height".format(files[0], np.round(z_ifc_mean)))

    # #########################################################################
    # Mask for the data to make analsis with RADOLAN and msevi comparable
    # Analyis should be made on the overlap region of 
    # 5.503-14.496 lon to 47.599-54.496 lat
    # #########################################################################
#    lon_mask = (lon[:] >= 5.503) & (lon[:] <= 14.496)
#    lat_mask = (lat[:] >= 47.599) & (lat[:] <= 54.496)

#    mask_2d = lon_mask[np.newaxis, :] & lat_mask[:, np.newaxis]

    # Add a third dimension to mask
#    mask_3d = np.asarray([mask_2d]*np.shape(time)[0])

#    wa  = np.ma.masked_array(wa,mask=~mask_3d) # ~ inverts the array (True->False and vice versa)
#    lon = np.ma.masked_array(lon,~lon_mask)
#    lat = np.ma.masked_array(lat,~lat_mask)
    
    
#    return time, lon, lat, z_ifc_mean, wa


##########################################################################
##########################################################################
 
if __name__ == '__main__':


    # meteosat test
    kws = dict(seg_field = 'meteosat', date = '20150704')

    # synsat test
    kws = dict(seg_field = 'synsat', date = '20150704', 
               dom = 'DOM03', 
               region_mode = 'shifted', 
               file_identifier = '-double_day')

    d = read_field(**kws)

