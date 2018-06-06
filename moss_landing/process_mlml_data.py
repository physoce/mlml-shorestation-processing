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

# path to final netcdf files
sw_csv_dir = os.path.join(sw_dir,'csv/')

# path to final netcdf files
sw_nc_dir = os.path.join(sw_dir,'netcdf/')
nc_prefix = 'moss_landing_'

run_download = False
run_makencfile = True

if run_download:
    
    # download csv files and convert to NetCDF format
    mlml.download_station_data(sw_csv_dir,'seawater',overwrite=True)
    
if run_makencfile:
    
    if not os.path.exists(sw_nc_dir):
        os.mkdir(sw_nc_dir)

    mlml.make_netcdf(sw_csv_dir,sw_nc_dir,nc_prefix,'seawater',download=False)

        