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
weather_dir = mlml_data_path.moss_landing_weather()

# path to final netcdf files
sw_csv_dir = os.path.join(sw_dir,'csv/')
weather_csv_dir = os.path.join(weather_dir,'csv/')

# path to final netcdf files
sw_nc_dir = os.path.join(sw_dir,'netcdf/')
nc_prefix = 'moss_landing_'
weather_nc_dir = os.path.join(weather_dir,'netcdf/')
weather_nc_prefix = 'moss_landing_weather_'

run_download = False
run_makencfile = True
run_download_weather = False
run_makencfile_weather = False

if run_download:

    # download csv files and convert to NetCDF format
    mlml.download_station_data(sw_csv_dir,'seawater',overwrite=True)

if run_makencfile:

    if not os.path.exists(sw_nc_dir):
        os.mkdir(sw_nc_dir)

    mlml.make_netcdf(sw_csv_dir,sw_nc_dir,nc_prefix,'seawater',download=False)
    
if run_download_weather:

    # download csv files and convert to NetCDF format
    mlml.download_station_data(weather_csv_dir,'weather',overwrite=False)

if run_makencfile_weather:

    if not os.path.exists(weather_nc_dir):
        os.mkdir(weather_nc_dir)

    mlml.make_netcdf(weather_csv_dir,weather_nc_dir,weather_nc_prefix,
                     'weather',download=False)
