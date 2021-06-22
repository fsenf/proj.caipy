#!/usr/bin/env python


def basis_setup(variable_name):


    if variable_name == 'wa':
        # segmentation parameters
        thres     = 1.0  # Threshold in m/s, only data is considered where wa>thres m/s
        msize     = 20   # Minium number of pixels for a cluster to be a cluster
        cl_method = 'watershed_merge' # Available segmentation methods: watershed, watershed_merge, connect

    elif variable_name in ['meteosat', 'synsat']:
        thres     = -240.0  # Threshold in K, but negative values
        msize     = 40   # Minium number of pixels for a cluster to be a cluster
        cl_method = 'watershed_merge' # Available segmentation methods: watershed, watershed_merge, connect

    elif variable_name in ['radolan', 'radvop']:
        thres     = 30.0  # Threshold in dBZ
        msize     = 30   # Minium number of pixels for a cluster to be a cluster
        cl_method = 'watershed_merge' # Available segmentation methods: watershed, watershed_merge, connect    


    return thres, msize, cl_method



