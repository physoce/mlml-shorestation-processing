# -*- coding: utf-8 -*-
import sys
import os
import mlml

sys.path.append('..')
import mlml_data_path

# get paths of data files
# paths to data should be modified in mlml_data_path.py
sw_dir = mlml_data_path.moss_landing()
nc_dir = mlml_data_path.netcdf()

nc_tempdir = os.path.join(nc_dir,'temp')

if not os.path.exists(nc_tempdir):
    os.mkdir(nc_tempdir)

sw_nc = os.path.join(nc_tempdir,'moss_landing_fromcsv.nc')

mlml.make_netcdf(sw_dir,sw_nc,'seawater',download=False)