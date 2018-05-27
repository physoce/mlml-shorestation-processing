# -*- coding: utf-8 -*-
import sys
import mlml

sys.path.append('..')
import mlml_data_path

# get paths of data files
# paths to data should be modified in mlml_data_path.py
sw_dir = mlml_data_path.seawater()
w_dir = mlml_data_path.weather()
nc_dir = mlml_data_path.netcdf()

sw_nc = nc_dir+'mlml_seawater.nc'
w_nc = nc_dir+'mlml_weather.nc'

ds = mlml.read_csv_data(sw_dir,format='dataset')

#mlml.make_netcdf(sw_dir,sw_nc,'seawater',download=True,overwrite=True)
#mlml.make_netcdf(w_dir,w_nc,'weather',download=True,overwrite=True)