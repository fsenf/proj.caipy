#!/usr/bin/env python


# Authors: 
# ======
#-------------------------------------------------------------------------
# Rieke Heinze  <rieke.heinze@mpimet.mpg.de>

# MPI
# Bundesstrasse 53
# 20146 Hamburg, Germany

#-------------------------------------------------------------------------
# Ieda Pscheidt <pscheidt@uni-bonn.de>

# Meteorological Institute, University of Bonn
# Auf dem Huegel 20
# 53121 Bonn

#------------------------------------------------------------------------
# Fabian Senf   <senf@tropos.de>

# TROPOS,
# Permoserstr. 15
# 04318 Leipzig, Germany. 

#------------------------------------------------------------------------


# #########################################################################
# Calculation of convective aggregation indices based on Fabian Senf's 
# clustering algorithm and Ieda Pscheidt's analysis routines 
# applied to radar and satellite observations and ICON-LEM simulations
# #########################################################################

# Load libraries ---------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy
import math
import glob
import xarray as xr
import sys
import os
import time as ti
import segmentation as seg # Fabians segmentation routines in lib/tropy
from pyproj import Proj
from skimage import measure
from mpl_toolkits.basemap import Basemap
from matplotlib import cm  #colormaps
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages

# Own libs
from utility_functions import stats, nne_size, boot, boot_obj_size
from io_functions import *
from indices import *
from segmentation_config import basis_setup

#rc('text', usetex=True) # need to load module texlive
rc('font', size=8)
rc('xtick', direction='out')
rc('ytick', direction='out')
rc('figure', figsize=(8.27,11.69)) #Set figure size to A4 as default



#------------------------------------------------------------------------
#------------------------------------------------------------------------
#     INFORMATION TO BE PROVIDED BY THE USER
#------------------------------------------------------------------------
# calc_cai.py:
#-------------
# INPUT FIELD: radolan, radvop, meteosat, synsat or wa
# ICON DOMAIN: DOM01, DOM02 or DOM03
# DATES: 
#    radar & satellite : available from 20140101 to 20151231
#    ICON simulations: 20150704,20140815, 20150705, 20140729 


# segmentation_config.py:
#-----------------------
# thres: threshold for the segmentation algorithm
# msize: minium number of pixels for a cluster to be a cluster
# cl_method: one of the available segmentation methods: watershed, watershed_merge, connect 
#            (connect is the same as used in Tompkins and Tobin et al.)

#------------------------------------------------------------------------
#------------------------------------------------------------------------


start = ti.time()


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
#                INPUTS PROVIDED BY THE USER
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------

# Enter the days for which the program should be run:
# date_set = ["20150701", "20150702", "20150703", "20150704", "20150705", "20150706", "20150707", "20150708", "20150709",
#              "20150710", "20150711", "20150712", "20150713", "20150714", "20150715", "20150716", "20150717", "20150718",
#              "20150720", "20150721", "20150722", "20150723", "20150724", "20150725", "20150726", "20150727",
#              "20150728", "20150729", "20150730", "20150731",
#              "20150801", "20150802", "20150803", "20150804", "20150805", "20150806", "20150807", "20150808", "20150809",
#              "20150810", "20150811", "20150812", "20150813", "20150814", "20150815", "20150816", "20150817", 
#              "20150820", "20150821", "20150822", "20150823", "20150824", "20150825", "20150826", "20150827",
#              "20150828", "20150829", "20150830", "20150831",
#              "20150901", "20150902", "20150903", "20150904", "20150905", "20150906", "20150907", "20150908", "20150909",
#              "20150910", "20150911", "20150912", "20150913", "20150915", "20150916", "20150917", "20150918",
#              "20150919", "20150920", "20150921", "20150922", "20150923", "20150924", "20150925", "20150926", "20150927",
#              "20150928", "20150929", "20150930"]:


#date_set = ['20140729', '20140815', '20150704', '20150705']
date_set = ['20140729']
file_identifier = '-redone_v1'


# Enter the field to be segmented: radolan, radvop, meteosat, synsat or wa                 
#seg_field = 'meteosat'
seg_field = 'synsat'


# Enter the domain of the ICON simulation. For meteosat and radolan enter dom = 'OBS'
dom  = 'DOM03'      


# Specify the height above ground of the field (ilev = 0 for obs. and synthetic fields) 
ilev = 0        


for date in date_set:

   print date


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------





# #########################################################################
#  Fixed input parameters
# #########################################################################

   # step for the time loop: 
   if seg_field in ['radolan']:
      step = 6
   else:
      step = 1

   # segmentation parameters
   thres, msize, cl_method = basis_setup(seg_field) 
   segmentation_para = dict(thres = thres, msize = msize, cl_method = cl_method)


   
   # file identifiers of ICON simulations: double_day, redone_v2
   # region mode of ICON simulations: default, shifted 
  
   if seg_field in ['radolan', 'radvop', 'meteosat']:
      region_mode= ''
      file_identifier = ''
   else:
      if date in ['20150705']:
         region_mode     = 'shifted'
         file_identifier = ''
      if date in ['20140815']:
         region_mode     = 'default'
         file_identifier = ''
      if date in ['20150704']:
         region_mode     = 'shifted'
         file_identifier = '-double_day'
      if date in ['20140729']:
         region_mode     = 'default'
         if file_identifier is None:
            file_identifier = '-redone_v2'
      if date in ['20160603']:
         region_mode     = 'default'
         file_identifier = ''

   
   input_para = dict(ilev = ilev, dom = dom, seg_field = seg_field, \
                  date = date, region_mode = region_mode, \
                  file_identifier = file_identifier )

   print input_para


# #########################################################################
# Enable plotting
# #########################################################################
   plot_bol = False

# #########################################################################
# Save Indices
# True: returns Ior, dist, UP, LW, Fo, Fs, Fm, Fmh, Fml
# False: returns only Ior
# #########################################################################
   indices_full_output= True

# #########################################################################
# Openening a projection object - needed for transforming lon/lat to 
# Cartesian grid
# #########################################################################

# UTM-Projection of zone 32 (most of Germany) for WGS84 which is a geodetic
# reference system for positioning on the Earth
   projection = Proj(proj='utm', zone=32, ellps='WGS84')


# #########################################################################
# Read in data-sets
# #########################################################################
   print('READING IN DATA-SETS')

   if seg_field in ['radolan', 'radvop']:
      times, lon, lat, lonv, latv, n_latlon, pos_grid, z_ifc_mean, fd, rain = read_field(**input_para) 
   else:
      times, lon, lat, lonv, latv, n_latlon, pos_grid, z_ifc_mean, fd = read_field(**input_para)


   # time dimension
   dim_time = np.shape(times)


   # time set up for runs
   time_loop= range(0, dim_time[0], step) # init, end, step 



# #########################################################################
# Initialize organization indices
# ########################################################################
   # Initializing organisation indices
   cop= []
   I_org = []
   I_sh = []
   N_cl = []
   D0 = []
   D1 = []
   SCAI = []
   UPPER_org=[]
   LOWER_org=[]
   Fobs=np.empty([len(time_loop), 1000])
   Fobs[:,:]=np.nan
   Fsim=np.empty([len(time_loop),100,1000])
   Fsim[:,:,:]=np.nan
   segmfield=np.empty([len(time_loop), lat.shape[0], lon.shape[1]]) # segmented field (wa, radolan, meteosat, synsat, radvop)
   segmfield[:,:,:]=np.nan
   labelfield=np.empty([len(time_loop), lat.shape[0], lon.shape[1]]) # labeled segmented field
   labelfield[:,:,:]=np.nan
   rainfall=np.empty([len(time_loop), lat.shape[0], lon.shape[1]])  # radolan objects-rainfall 
   rainfall[:,:,:]=np.nan
   mrain = np.empty(len(time_loop))
   mrain[:]= np.nan
   Tim = np.zeros(len(time_loop))



   # Convert time-axis to a human readable version
   if seg_field in ['radolan', 'radvop']:
      dates, frac_day, datetime_obj= convert_gregorian_time_to_readable_time(times, date)
      times = convert_gregorian_time_to_frac_time(times, date) 
   else:
      dates, frac_day, datetime_obj= convert_frac_time_to_readable_time(times)


   # either time lists or still xarray object ...
   try:
      times = times.values
      dates = dates.values
   except:
       pass

   # meshgrid lat and lon:  needed for contour-plots
   lon_mg=lon
   lat_mg=lat

   # Convert lon/lat to Cartesian coordinates based on projection 
   x_lon, y_lat = projection(lon_mg,lat_mg)  # in meters


   # Create the grid to work on
   external_grid = (x_lon, y_lat) 

   #Mean delta_x and delta_y of the grid
   dx_mean = np.mean(abs(np.diff(x_lon, axis=-1)))
   dy_mean = np.mean(abs(np.diff(y_lat, axis=0)))


   # Size of a grid cell (pixel) in km - approximated, as the actual size depends
   # on the latitude
   area_grid_pixel = dx_mean*dy_mean # in m**2

   print dx_mean
   print dy_mean

# #########################################################################
# Preparation for calculation of SCAI
# #########################################################################

   # N_max: potential maximum number of objects in the domain, take the total 
   # number of grid points in the analysis domain
   N_max = n_latlon

   print 'N_max= ', N_max

   # Determine the lenght L of the characteristic domain size in km as 
   # diagonal of the domain size

   aux_x=np.empty([lat.shape[0], lon.shape[1]])
   aux_y=np.empty([lat.shape[0], lon.shape[1]])
   aux_x[:,:]=np.nan
   aux_y[:,:]=np.nan
   aux_x[pos_grid]=x_lon[pos_grid] 
   aux_y[pos_grid]=y_lat[pos_grid]
   x_min= np.nanmin(aux_x)
   y_min= np.nanmin(aux_y)
   x_max= np.nanmax(aux_x)
   y_max= np.nanmax(aux_y)

   L = ((x_max-x_min)**2+(y_max-y_min)**2)**(0.5)/1000.    # in km

   print 'L = {} km \nN_max = {}'.format(L,N_max)

   print 'x_min, x_max =', x_min, x_max
   print 'y_min, y_max =', y_min, y_max 


# #########################################################################
# Preparation of Output Filenames
# #########################################################################


   variable_name =  input_para['seg_field']
   out_directory = '../../data/indices/%s' % variable_name

   out_basename = 'cai_{}_{}_{}{}_{}m'.format(variable_name, dom,\
           date, file_identifier,int(z_ifc_mean))


   file_name_out = '%s/plot_%s' % (out_directory, out_basename)
   ds_name_out =  '%s/%s' % (out_directory, out_basename)


# #########################################################################
# Segmentation and basic properties of the objects
# #########################################################################

   with PdfPages(file_name_out+'.pdf') as pdf: # open pdf file outside the loop

     count=-1
     for nt in time_loop:
        print("Timestep ",nt, " ",datetime_obj[nt] )
        count=count+1
        Tim[count]=times[nt]
        print('   SEGMENTATION')


     # Segmentation :
     # Object labels are ordered according to size, the object with the largest
     # amount of pixels gets label 1, ...
        bin_field = fd[nt,:,:]
        bin_field[np.where(np.isnan(bin_field))]=-999.9
        bin_field[bin_field < thres] = 0

        if seg_field in ['meteosat', 'synsat']:
           bin_field = -bin_field
           bin_field[np.where(bin_field == -0.)] = 0. 

        fd_clusters = seg.clustering(bin_field, 0, cluster_method = cl_method, min_size=msize) 
       

      # Total number of clusters/objects
        N_c = np.max(fd_clusters)
        N_cl.append(N_c)
        print('   Total number of objects: {}'.format(N_c))

        if  N_c > 0:

         # save the segmented field
           sfield=bin_field
           sfield[np.where(sfield == 0.)]=np.nan
           segmfield[count,:,:]=sfield

         # save the labeled field
           fd_lab=fd_clusters.astype('float')
           fd_lab[np.where(fd_lab == 0.)]=np.nan
           labelfield[count,:,:] = fd_lab  

         # RADOLAN : Save the spatial mean rainfall of the objects 
           if seg_field in ['radolan', 'radvop']:
              rfll = rain[nt,:,:]
              rfll[np.where(np.isnan(fd_lab))]=np.nan
              rfll[np.where(rfll < 0.)]=np.nan
              rainfall[count,:,:]=rfll
              mrain[count] = np.nanmean(rfll)   # mean of field  
              if(len(np.where(rfll <0.)[0]) > 0):
                 sys.exit('Negative rainfall values')
  
 
         # Deduce the coordinates of the object's centers of mass (geometric centroid)
         # np.unique gives sorted unique elements of input: [0,1,2,...,N_c]
         # lon_clm is a vector containing the mean longitude of object 1,2,3, ... N_c:
         # lon_clm =[lon(obj1),lon(obj2), ..., lon(objN_c)]
         # index=0 is not considered since this value represents the background and not an object
           index   = np.unique(fd_clusters)
           lon_clm = scipy.ndimage.measurements.mean(lon_mg, labels = fd_clusters, index = index[1:])
           lat_clm = scipy.ndimage.measurements.mean(lat_mg, labels = fd_clusters, index = index[1:])

         # Number of pixels (grid points) in each object (uses stats-function from functions_ieda.py)
         # N_clm is a vector containing the number of pixels of object with label 1:
         # N_clm = [N_c(obj1), N_c(obj2), ..., N_c(objN_c)]
         # Gives the same result as measure.regionprops(fd_clusters)
           N_clm = stats(lon_mg, labels = fd_clusters, index = index[1:])[0] 


         # Area of the objects in m**2
           AreaObj= N_clm*area_grid_pixel

         # Equivalent radius (r_equ) of the objects in km
         # Assume the area is equal to that of a circle with radius=r_equ
           r_equ = (math.pi**(-0.5))*((N_clm*area_grid_pixel)**(0.5))/1000. # in km

         # Convert coordinates of centers of mass of the objects to a Cartesian
         # grid by means of projection p
         # Only on Cartesian coordinates the distance can be obtained with 
         # scipy.spatial.distance.cdist
         # Calling a Proj class instance with the arguments lon, lat will convert 
         # lon/lat (in degrees) to x/y native map projection coordinates (in meters).
           x_clm, y_clm = projection(lon_clm,lat_clm)


         # #########################################################################
         # SHAPE INDEX - I.shape - defined as sqrt(A)/(0.282*P) with A=area of the
         # object and P=perimeter. 0<= I_shape <= 1.0. I_shape=0 is the shape of a 
         # line which has zero area and I_shape=1 is the shape of a perfect circle 
         # => 0.5 is something like an ellipsoid
         # #########################################################################
          
           print '   SHAPE INDEX'
          
           if plot_bol:
               I_shape_mean, I_shape = calculate_Ishape(fd_clusters, full_output = True)
           else:
               I_shape_mean = calculate_Ishape(fd_clusters)

           print I_shape_mean
           I_sh.append(I_shape_mean)   


         # #########################################################################
         #  I.org 
         #  Nearest Neighbor Cumulative Distribution Function (NNCDF) 
         # #########################################################################

         # Compute nearest neighbor distance between the centers of mass for each
         # pair of objects - use x/y representation in m
           if N_c > 1: # there need to be at least two objects

              print('computing NNCDF and I.org ')

              if indices_full_output:
                 Ior, UP, LW, dist, Fobs[count,:], Fsim[count,:,:], Fm, Fmh, Fml = calculate_Iorg(x_clm, y_clm, r_equ, \
                                             external_grid, max_dist = 1000, \
                                             dist_interval = 1. ,number_of_repetitions = 100, \
                                             full_output = True)

              else:
                 Ior = calculate_Iorg(x_clm, y_clm, r_equ, \
                                             external_grid, max_dist = 1000, \
                                             dist_interval = 1. ,number_of_repetitions = 100, \
                                             full_output = False)
         # ##########################################################################       
         #  SCAI
         # ##########################################################################
            
              print('computing SCAI ')
              indx, d0_clm, d1_clm = calculate_SCAI(x_clm, y_clm, N_max= N_max, L = L)


              print('computing COP ')
              cop_id = calculate_COP_fast(x_clm/1000., y_clm/1000., AreaObj/1000000.) # in km


           else: # if N_c=1 

              d0_clm = np.nan
              d1_clm = np.nan
              indx = np.nan
              Ior = np.nan
              UP = np.nan
              LW = np.nan
              cop_id = np.nan

           # Collect all the stuff  
           D0.append(d0_clm)   
           D1.append(d1_clm)   
           I_org.append(Ior)
           UPPER_org.append(UP)
           LOWER_org.append(LW)
           SCAI.append(indx)
           cop.append(cop_id)

         # #########################################################################
         #  Plotting - only when at least two objects were identified
         # #########################################################################
          
           if (plot_bol and N_c > 1):
              print('PLOTTING')

              #plt.close('all')

            # Values for x and y-axis labels
              loncoords = np.arange( 4., 16., 2.)
              latcoords = np.arange(45., 60., 2.)

            # Contourplots on a map (llm = lat-lon-map)
              lon_min = np.min(lon)#lon.min()
              lon_max = np.max(lon)
              lat_max = np.max(lat)
              lat_min = np.min(lat)

            # Setting up the map: cyl=equidistand cylindrical projection
            # Possible setting for resolution: c, l, i, h, f
            #llm = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, \
            #              urcrnrlon=lon_max, urcrnrlat=lat_max, resolution='l')

            # Setting up the map: lcc=Lambert Conformal projection as Ieda uses it
              llm = Basemap(projection='lcc',llcrnrlon=lon_min, llcrnrlat=lat_min, \
                            urcrnrlon=lon_max, urcrnrlat=lat_max, resolution='l',  \
                            lon_0=9.4838500000000003, lat_0=51.093499999999999)

            # Project of lon_mg and lat_mg to the projection selected in basemap
              lon_llm,lat_llm = llm(lon_mg,lat_mg)

            # #########################################################################
            # Contourplot of fd
            # #########################################################################

            # Choose plotting range
              range_mintomax = np.linspace(np.min(fd),np.max(fd),21)
              range_centered = np.linspace(-2.0,2.,21)

              plt_range = range_centered

              plt.figure() 
              plt.subplot(4,2,1)

              c_plot_fd = llm.contourf(lon_llm,lat_llm,fd[nt,:,:],plt_range,cmap=cm.PiYG_r) #coolwarm

            # Appearance of map
              llm.drawcountries(color='black', linewidth=0.5)
              llm.drawcoastlines(color='black', linewidth=0.5)
              llm.drawmeridians(loncoords, labels=[0,0,0,1], linewidth=0.1)
              llm.drawparallels(latcoords, labels=[1,0,0,0], linewidth=0.1)

            # Colorbar
              cb_fd=plt.colorbar(c_plot_fd,orientation='horizontal',shrink=0.95)
              cb_fd.set_label('vertical velocity in m/s')

            # Axis labels
              #plt.xlabel('lon', labelpad=10)
              #plt.ylabel('lat', labelpad=20)


            # #########################################################################
            # Contourplot of clusters from fd (fd_clm)
            # #########################################################################
              plt.subplot(4,2,2)
              c_plot_fd_clm = llm.contourf(lon_llm,lat_llm,fd_clm[:,:],cmap=cm.rainbow_r)

            # Appearance of map
              llm.drawcountries(color='black', linewidth=0.5)
              llm.drawcoastlines(color='black', linewidth=0.5)
              llm.drawmeridians(loncoords, labels=[0,0,0,1], linewidth=0.1)
              llm.drawparallels(latcoords, labels=[1,0,0,0], linewidth=0.1)

            # Colorbar
              cb_clm=plt.colorbar(c_plot_fd_clm,orientation='horizontal',shrink=0.95)
              cb_clm.set_label('Object number')

            # Axis labels
              #plt.xlabel('lon', labelpad=10)
              #plt.ylabel('lat', labelpad=20)

            # #########################################################################
            # Plot of the CDF of dist_size (NN distances with respect to actual size
            # of an object
            # #########################################################################
              ax = plt.subplot(4,2,3) # axis object controlling the axes
              l_plot_Fo = ax.plot(dist, Fo)

            # Setting line properties
              plt.setp(l_plot_Fo, color='black', linewidth=1.5)

            # Limit of x-axis
              ax.set_xlim(0.,800.)

            # Axis labels
              ax.set_xlabel('NN-distance in km')
              ax.set_ylabel('actual NNCDF')

            # Make right and top axis invisible
              ax.spines['right'].set_visible(False)
              ax.get_yaxis().tick_left()
              ax.spines['top'].set_visible(False)
              ax.get_xaxis().tick_bottom()

            # Place title of figure in the center where is more space
              title_string_for_plot = '{}, {}, z={} m'.format(dom,datetime_obj[nt],np.round(z_ifc_mean)) 
              ax.annotate(title_string_for_plot, xy=(-0.5,1.1), xytext=(-0.5,1.1),\
                          xycoords='axes fraction', textcoords='axes fraction',fontsize=13)

            # #########################################################################
            # Plot of Fo vs. Fm, so observed NNCDF vs. inhibition NNCDF
            # #########################################################################
              ax = plt.subplot(4,2,4) # axis object controlling the axes

            # diagonal line
              l_plot_diag = ax.plot([0.,1.], [0.,1.], ':')
              plt.setp(l_plot_diag, color='black', linewidth=1.5)

            # NNCDF
              l_plot_Fo_Fm = ax.plot(Fm, Fo,'-')
              plt.setp(l_plot_Fo_Fm, color='black', linewidth=1.5, label='actual NNCDF')

            # confidence interval
              l_plot_Fo_Fm_conf = ax.plot(Fmh, Fo,'--' ,Fml, Fo, '--')
              plt.setp(l_plot_Fo_Fm_conf, color='red', linewidth=0.5)

            # Limit of x-axis
              ax.set_xlim(0.,1.0)
              ax.set_ylim(0.,1.0)

            # Axis labels
              ax.set_xlabel('inhibition NNCDF')
              ax.set_ylabel('actual NNCDF')

            # Make  top axis invisible
              ax.spines['top'].set_visible(False)
              ax.get_xaxis().tick_bottom()

            #########################################################################
            # Plot on right y-axis
              ax2 = ax.twinx()

            # distance dist(r)
              l_plot_dist = ax2.plot(Fmh, dist,'-', label='NN-distance' )
              plt.setp(l_plot_dist, color='blue', linewidth=1.5)

            # Limit of x-axis
              ax2.set_ylim(0.,800.0)
              ax2.set_xlim(0.,1.0)

            # Axis labels
              ax2.set_ylabel('distance in km')

            #########################################################################
            # Legend
              lines = l_plot_Fo_Fm+l_plot_dist 
              labels = [l.get_label() for l in lines]
              ax.legend(lines, labels, fontsize=8, loc=0)


            # #########################################################################
            # Plot of histogram of number and size of the objects
            # #########################################################################
              ax = plt.subplot(4,2,5)

            # Determine number of bins in dependence on a fixed bin size
              bin_size = 0.5 # Bin width in km
              min_edge = math.floor(np.min(r_equ))
              max_edge = math.ceil(np.max(r_equ))
              Nbins = (max_edge-min_edge)/bin_size
              bin_list = np.linspace(min_edge,max_edge, Nbins+1)

              hist_plot_N = ax.hist(r_equ,bins=bin_list,color='purple')

            # Axis labels
              ax.set_xlabel('equivalent radius in km')
              ax.set_ylabel('frequency')

            # Plot titles
              plt.title('Object size')

            # Make right and top axis invisible
              ax.spines['right'].set_visible(False)
              ax.get_yaxis().tick_left()
              ax.spines['top'].set_visible(False)
              ax.get_xaxis().tick_bottom()

            # #########################################################################
            # Plot of histogram of I_shape
            # #########################################################################
              ax = plt.subplot(4,2,6)

            # Determine number of bins in dependence on a fixed bin size
              bin_size = 0.1 # Bin width
              min_edge = math.floor(np.min(I_shape))
              max_edge = math.ceil(np.max(I_shape))
              Nbins = (max_edge-min_edge)/bin_size
              bin_list = np.linspace(min_edge,max_edge, Nbins+1)

              hist_plot_N = ax.hist(I_shape,bins=bin_list,color='purple')

            # Axis labels
              ax.set_xlabel('I.shape')
              ax.set_ylabel('frequency')

            # Plot titles
              plt.title('Shape index')

            # Make right and top axis invisible
              ax.spines['right'].set_visible(False)
              ax.get_yaxis().tick_left()
              ax.spines['top'].set_visible(False)
              ax.get_xaxis().tick_bottom()

            # #########################################################################
            # Plot of I_shape vs. r_equivalent
            # #########################################################################
              ax = plt.subplot(4,2,7) # axis object controlling the axes

            # NNCDF
              l_plot_reff_Ishape = ax.plot(r_equ, I_shape,'*')
              plt.setp(l_plot_reff_Ishape, color='black',)

            # Limit of x-axis
              #ax.set_xlim(0.,1.0)
              #ax.set_ylim(0.,1.0)

            # Axis labels
              ax.set_xlabel('equivalent radius in km')
              ax.set_ylabel('shape index I.shape')

            # Make right and top axis invisible
              ax.spines['right'].set_visible(False)
              ax.get_yaxis().tick_left()
              ax.spines['top'].set_visible(False)
              ax.get_xaxis().tick_bottom()

            # Put some results in the graphs
              tx = 'N={:4d} \nDO={:6.2f} km '\
                   '\nL={:5.1f} km \nN\_max={:6d} \nSCAI={:6.4f} \nr\_equ={:6.2f} km'\
                   '\nI.org={:6.4f} \nI.shape={:6.4f}'\
                   '\nmsize = {} px \nthres={:4.2f} m/s'\
                   '\nmethod = {}'.format(N_c, d0_clm, L, N_max,indx, np.mean(r_equ),\
                                       Ior, I_shape_mean, msize, thres, cl_method )


              ax.annotate(tx, xy=(1.2,0.4), xytext=(1.2,0.4),xycoords='axes fraction', textcoords='axes fraction')

            # #########################################################################
            # Create final plot 
            # #########################################################################
              #plt.subplot_tool()

              #plt.suptitle(title_string_for_plot, \
              #             fontsize=10) 

              #adjusts subplot params to fit into figure size
              plt.tight_layout(w_pad=1.0) # as fraction of fontsize

              #plt.show()

              # saves the current figure into a pdf page
              pdf.savefig()  

              # saves the current figure as png
              #plt.savefig('{}_{}_{}{}s_{}m.png'.format(dom,dates[nt].values,file_identifier,\
              #                  frac_day[nt].values,int(z_ifc_mean)),bbox_inches='tight')

              plt.close()
        else:
          d0_clm = np.nan
          D0.append(d0_clm)
          d1_clm = np.nan
          D1.append(d1_clm)
          indx = np.nan
          SCAI.append(indx)
          I_shape_mean = np.nan
          I_sh.append(I_shape_mean)
          Ior = np.nan
          I_org.append(Ior)
          UP= np.nan
          UPPER_org.append(UP)
          LW= np.nan           
          LOWER_org.append(LW)
          cop_id=np.nan
          cop.append(cop_id)


   print 'I_org= ', I_org
   print 'I_sh =', I_sh
   print 'UPPER_org =', UPPER_org
   print 'LOWER_org =', LOWER_org
   print 'SCAI =', SCAI
   print 'COP=', cop

# ###############################################################################################
# SAVE THE RESULTS IN A NETCDF FILE
# ###############################################################################################
   print('WRITING DATA TO NETCDF FILE')

# Replacing NaN with FillValue
   fi_val = -999.9


   D0_fv = np.where(~np.isnan(D0),D0,fi_val)
   D1_fv = np.where(~np.isnan(D1),D1,fi_val)
   N_cl_fv = np.where(~np.isnan(N_cl),N_cl,fi_val)
   SCAI_fv = np.where(~np.isnan(SCAI),SCAI,fi_val)
   cop_fv = np.where(~np.isnan(cop),cop,fi_val)
   I_sh_fv = np.where(~np.isnan(I_sh),I_sh,fi_val)
   I_org_fv = np.where(~np.isnan(I_org),I_org,fi_val)
   UPPER_org_fv = np.where(~np.isnan(UPPER_org),UPPER_org, fi_val)
   LOWER_org_fv = np.where(~np.isnan(LOWER_org),LOWER_org, fi_val)
   Fobs_fv = np.where(~np.isnan(Fobs),Fobs, fi_val)
   Fsim_fv = np.where(~np.isnan(Fsim),Fsim, fi_val)
   segmfield_fv = np.where(~np.isnan(segmfield),segmfield, fi_val)
   labelfield_fv = np.where(~np.isnan(labelfield),labelfield,fi_val)
   mrain_fv = np.where(~np.isnan(mrain),mrain,fi_val)
   rainfall_fv = np.where(~np.isnan(rainfall),rainfall,fi_val)
   times_fv = np.where(~np.isnan(Tim),Tim,fi_val) # time in frac days %Y%M%D.f


   print I_org_fv



# Global attributes for the data set:
# -----------------------------------
   if seg_field in ['wa']:

      var= 'wa'

      att_glob={'author': 'Rieke Heinze (rieke.heinze@mpimet.mpg.de)',\
          'institution': 'Max-Planck-Institute for Meteorology',\
          'description': 'Convective aggregation metrics from ICON-LEM output',\
          'segmentation': 'Method {} based on Fabian Senfs Python implementation'.format(cl_method),\
          'msize': 'Minimum size of an object: {} px'.format(msize),\
          'thres': 'Threshold in original field: {} m/s'.format(thres),\
          'segmented_field': seg_field, \
          'ICON_model_level': '{} m'.format(z_ifc_mean),\
          'ICON_model_domain': dom }

   if seg_field in ['radolan']:

      var = 'Z'

      att_glob={'author': 'Ieda Pscheidt (pscheidt@uni-bonn.de)',\
                'institution': 'Meteorological Institute, University of Bonn',\
                'description': 'Convective aggregation metrics for radar reflectivity factor',\
                'segmentation': 'Method {} based on Fabian Senfs Python implementation'.format(cl_method),\
                'msize': 'Minimum size of an object: {} px'.format(msize),\
                'thres': 'Threshold in original field: {} dBZ'.format(thres),\
                'segmented_field': seg_field }                

      att_field = {'units': 'dBZ',\
             'long_name': 'radar reflectivity factor - segmented field',
             '_FillValue': fi_val}

      att_rainfall = {'units': 'mm/h',\
                'long_name': 'objects radar rainfall rate',
                '_FillValue': fi_val}


      att_mean_rain = {'units': 'mm/h',\
                 'long_name': 'spatial mean rainfall',\
                 '_FillValue': fi_val}


   if seg_field in ['radvop']:

      var = 'Z'

      att_glob={'author': 'Ieda Pscheidt (pscheidt@uni-bonn.de)',\
                'institution': 'Meteorological Institute, University of Bonn',\
                'description': 'Convective aggregation metrics for synthetic radar reflectivity factor',\
                'segmentation': 'Method {} based on Fabian Senfs Python implementation'.format(cl_method),\
                'msize': 'Minimum size of an object: {} px'.format(msize),\
                'thres': 'Threshold in original field: {} dBZ'.format(thres),\
                'segmented_field': seg_field, \
                'ICON_model_domain': dom }

      att_field = {'units': 'dBZ',\
             'long_name': 'synthetic radar reflectivity factor - segmented field',
             '_FillValue': fi_val}


      att_rainfall = {'units': 'mm/h',\
                'long_name': 'simulations objects rainfall rate',
                '_FillValue': fi_val}


      att_mean_rain = {'units': 'mm/h',\
                 'long_name': 'spatial mean rainfall',\
                 '_FillValue': fi_val}



   if seg_field in ['meteosat']:

      var = 'BT'

      att_glob={'author': 'Fabian Senf (senf@tropos.de)',\
                'institution': 'Leibniz Institute for Tropospheric Research',\
                'description': 'Convective aggregation metrics for Meteosat brightness temperature',\
                'segmentation': 'Method {} based on Fabian Senfs Python implementation'.format(cl_method),\
                'msize': 'Minimum size of an object: {} px'.format(msize),\
                'thres': 'Threshold in original field: {} K'.format(thres),\
                'segmented_field': seg_field }

      att_field = {'units': 'K',\
             'long_name': 'satellite brightness temperature - segmented field',
             '_FillValue': fi_val}


   if seg_field in ['synsat']:

      var = 'BT'

      att_glob={'author': 'Fabian Senf (senf@tropos.de)',\
                'institution': 'Leibniz Institute for Tropospheric Research',\
                'description': 'Convective aggregation metrics for synthetic brightness temperature',\
                'segmentation': 'Method {} based on Fabian Senfs Python implementation'.format(cl_method),\
                'msize': 'Minimum size of an object: {} px'.format(msize),\
                'thres': 'Threshold in original field: {} K'.format(thres),\
                'segmented_field': seg_field, \
                'ICON_model_domain': dom } 

      att_field = {'units': 'K',\
             'long_name': 'synthetic brightness temperature - segmented field',
             '_FillValue': fi_val}




# Attributes for the single variables:
#-------------------------------------


   att_time = {'units': 'day as %Y%m%d.%f',\
               'long_name': 'Time',\
               'calendar': 'proleptic_gregorian'}


   #att_time = {'units': 'seconds since 1970-01-01 00:00:00',\
   #         'long_name': 'time',\
   #          '_FillValue': fi_val}

   att_lon = {'units': 'degrees_east',\
           'long_name': 'longitude',\
           'axis': 'X',\
           '_FillValue': fi_val}
   att_lat = {'units': 'degrees_north',\
           'long_name': 'latitude',\
           'axis': 'Y',\
           '_FillValue': fi_val}
   att_dist = {'units': 'km',\
            'long_name': 'distance',\
            '_FillValue': fi_val}
   att_nsimul = {'units': '-',\
              'long_name': 'number of simulations',\
              '_FillValue': fi_val}
   att_field_labeled = {'units': '',\
                     'long_name': 'labeled field',
                     '_FillValue': fi_val}
   att_D0 = {'units': 'km',\
         'long_name': 'geometrical mean distance',\
         '_FillValue': fi_val}
   att_D1 = {'units': 'km',\
         'long_name': 'arithmetical mean distance',\
         '_FillValue': fi_val}
   att_N_cl = {'units': '-',\
         'long_name': 'Number of objects',\
         '_FillValue': fi_val}
   att_SCAI = {'units': '-',\
          'long_name': 'simple convective aggregation index', \
          'reference': 'Tobin et. al.(2012)',\
          'L': L,\
          'Nmax': N_max,\
          '_FillValue': fi_val}
   att_cop = {'units': '-',\
          'long_name': 'convective organization potential', \
          'reference': 'white et al.(2018)',\
          '_FillValue': fi_val}

   att_I_sh = {'units': '-',\
           'long_name': 'shape index',\
           'reference': 'Maceachren (1995); Xia (1996)',\
           '_FillValue': fi_val}
   att_I_org = {'units': '-',\
           'long_name': 'organisation index',\
            'reference': 'Nair et. al.(1998), Tompkins and Semie (2017), Weger et. al. (1992)',\
           '_FillValue': fi_val}
   att_UPPER_org = {'units': '-',\
                'long_name': 'organisation index + 1.96STD ',\
                 'reference': 'Nair et. al.(1998), Tompkins and Semie (2017), Weger et. al. (1992)',\
                 '_FillValue': fi_val}
   att_LOWER_org = {'units': '-',\
                 'long_name': 'organisation index - 1.96STD',\
                 'reference': 'Nair et. al.(1998), Tompkins and Semie (2017), Weger et. al. (1992)',\
                 '_FillValue': fi_val}

   att_Fobs = {'units': '-',\
       	    'long_name': 'observed nearest neighbor cumulative distribution function',\
	    '_FillValue': fi_val}  
   att_Fsim = {'units': '-',\
	    'long_name': 'simulated nearest neighbor cumulative distribution function',\
	    '_FillValue': fi_val}


# Create the data set:
#---------------------

   if seg_field in ['radolan', 'radvop']:
      ds_out = xr.Dataset({'D0': (['time'], D0_fv, att_D0),\
                     '{}'.format(var): (['time','lat','lon'], segmfield_fv[:,:,:], att_field),\
                     '{}_labeled'.format(var): (['time','lat','lon'], labelfield_fv[:,:,:], att_field_labeled),\
	             'rainfall': (['time','lat','lon'], rainfall_fv[:,:,:], att_rainfall),\
                     'mean_rain': (['time'], mrain_fv, att_mean_rain),\
                     'D0': (['time'], D0_fv, att_D0),\
                     'D1': (['time'], D1_fv, att_D1),\
                     'N_cl': (['time'], N_cl_fv, att_N_cl),\
                     'SCAI': (['time'], SCAI_fv, att_SCAI),\
                     'COP': (['time'], cop_fv, att_cop),\
                     'I_sh': (['time'], I_sh_fv, att_I_sh),\
                     'I_org': (['time'], I_org_fv, att_I_org),\
                     'UPPER_org': (['time'], UPPER_org_fv, att_UPPER_org),\
	             'LOWER_org': (['time'], LOWER_org_fv, att_LOWER_org),\
                     'Fobs': (['time', 'dist'], Fobs_fv[:,:], att_Fobs),\
	             'Fsim': (['time', 'nsimul', 'dist'], Fsim_fv[:,:,:], att_Fsim)},\
      coords={'time': (['time'], times_fv, att_time),
      'latitude': (['lat'], latv[:], att_lat),
      'longitude': (['lon'], lonv[:], att_lon)},\
      attrs=att_glob)

   else:

  # if (seg_field == 'radvop'):
      ds_out = xr.Dataset({'D0': (['time'], D0_fv, att_D0),\
                     '{}'.format(var): (['time','lat','lon'], segmfield_fv[:,:,:], att_field),\
                     '{}_labeled'.format(var): (['time','lat','lon'], labelfield_fv[:,:,:], att_field_labeled),\
                     'D0': (['time'], D0_fv, att_D0),\
                     'D1': (['time'], D1_fv, att_D1),\
                     'N_cl': (['time'], N_cl_fv, att_N_cl),\
                     'SCAI': (['time'], SCAI_fv, att_SCAI),\
                     'COP': (['time'], cop_fv, att_cop),\
                     'I_sh': (['time'], I_sh_fv, att_I_sh),\
                     'I_org': (['time'], I_org_fv, att_I_org),\
                     'UPPER_org': (['time'], UPPER_org_fv, att_UPPER_org),\
                     'LOWER_org': (['time'], LOWER_org_fv, att_LOWER_org),\
                     'Fobs': (['time', 'dist'], Fobs_fv[:,:], att_Fobs),\
                     'Fsim': (['time', 'nsimul', 'dist'], Fsim_fv[:,:,:], att_Fsim)},\
      coords={'time': (['time'], times_fv, att_time),
      'latitude': (['lat'], latv[:], att_lat),
      'longitude': (['lon'], lonv[:], att_lon)},\
      attrs=att_glob)



  # if seg_field in ['meteosat', 'synsat']:
  #    ds_out = xr.Dataset({'D0': (['time'], D0_fv, att_D0),\
  #                   'BT': (['time','lat','lon'], segmfield_fv[:,:,:], att_field),\
  #                   'BT_labeled': (['time','lat','lon'], labelfield_fv[:,:,:], att_field_labeled),\
  ##                   'D0': (['time'], D0_fv, att_D0),\
  #                   'D1': (['time'], D1_fv, att_D1),\
  #                   'N_cl': (['time'], N_cl_fv, att_N_cl),\
  #                   'SCAI': (['time'], SCAI_fv, att_SCAI),\
  #                   'I_sh': (['time'], I_sh_fv, att_I_sh),\
  #                   'I_org': (['time'], I_org_fv, att_I_org),\
  #                   'UPPER_org': (['time'], UPPER_org_fv, att_UPPER_org),\
  #                   'LOWER_org': (['time'], LOWER_org_fv, att_LOWER_org),\
  #                   'Fobs': (['time', 'dist'], Fobs_fv[:,:], att_Fobs),\
  #                   'Fsim': (['time', 'nsimul', 'dist'], Fsim_fv[:,:,:], att_Fsim)},\
  #    coords={'time': (['time'], times_fv, att_time),
  #    'latitude': (['lat'], latv[:], att_lat),
  #    'longitude': (['lon'], lonv[:], att_lon)},\
  #    attrs=att_glob)

  # if seg_field in ['wa']:
  #    ds_out = xr.Dataset({'D0': (['time'], D0_fv, att_D0),\
  #                   'wa': (['time','lat','lon'], segmfield_fv[:,:,:], att_field),\
  #                   'wa_labeled': (['time','lat','lon'], labelfield_fv[:,:,:], att_field_labeled),\
  #                   'D0': (['time'], D0_fv, att_D0),\
  #                   'D1': (['time'], D1_fv, att_D1),\
  #                   'N_cl': (['time'], N_cl_fv, att_N_cl),\
  #                   'SCAI': (['time'], SCAI_fv, att_SCAI),\
  #                   'I_sh': (['time'], I_sh_fv, att_I_sh),\
  #                   'I_org': (['time'], I_org_fv, att_I_org),\
  #                   'UPPER_org': (['time'], UPPER_org_fv, att_UPPER_org),\
  #                   'LOWER_org': (['time'], LOWER_org_fv, att_LOWER_org),\
  #                   'Fobs': (['time', 'dist'], Fobs_fv[:,:], att_Fobs),\
  #                   'Fsim': (['time', 'nsimul', 'dist'], Fsim_fv[:,:,:], att_Fsim)},\
  #    coords={'time': (['time'], times_fv, att_time),
  #    'latitude': (['lat'], latv[:], att_lat),
  #    'longitude': (['lon'], lonv[:], att_lon)},\
  #    attrs=att_glob)






# Write into a netcdf file
   ds_out.to_netcdf(ds_name_out+'.nc')

   end = ti.time()

   total_elapsed_time = (end - start) /60. # in minutes
   print("Elapsed total time for execution: {} min".format(total_elapsed_time)) 
 
