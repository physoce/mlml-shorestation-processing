# -*- coding: utf-8 -*-
import sys
import mlml

sys.path.append('..')
import mlml_data_path

# get paths of data files
# paths to data should be modified in mlml_data_path.py
sw_dir = mlml_data_path.seawater()
w_dir = mlml_data_path.weather()

# download csv files and convert to NetCDF format
mlml.download_station_data(sw_dir,'seawater',overwrite=True)
mlml.download_station_data(sw_dir,'weather',overwrite=True)



