# -*- coding: utf-8 -*-
import sys
import os
import mlml
import numpy as np
import xarray as xr

sys.path.append('..')
import mlml_data_path

# get paths of data files
# paths to data should be modified in mlml_data_path.py
sw_dir = mlml_data_path.moss_landing()
nc_dir = mlml_data_path.netcdf()

# paths to intermediate files
nc_dir = os.path.join(nc_dir,'moss_landing')
nc_prefix = 'moss_landing_'

run_download = False
run_makencfile = True

if run_download:
    
    # download csv files and convert to NetCDF format
    mlml.download_station_data(sw_dir,'seawater',overwrite=True)
    
if run_makencfile:
    
    if not os.path.exists(nc_dir):
        os.mkdir(nc_dir)

    mlml.make_netcdf(sw_dir,nc_dir,nc_prefix,'seawater',download=False)

        