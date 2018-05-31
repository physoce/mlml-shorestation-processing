# -*- coding: utf-8 -*-

"""
edit this file to define paths to MLML data

example:
-------

import mlml_data_path

data_dir = mlml_data_path.seawater()
"""

import os

# top-level data directory (default: ~/work/Data/MLML_shore_stations/)
def mlml_data():
    homedir = os.path.expanduser('~')
    datadir = os.path.join(homedir,'work/Data/MLML_shore_stations/')
    return datadir
    
# csv files 
def moss_landing():
    datadir = os.path.join(mlml_data(),'csv/moss_landing/')
    return datadir

def moss_landing_weather():
    datadir = os.path.join(mlml_data(),'csv/moss_landing_weather/')
    return datadir

# NetCDF files (created from csv files in get_mlml_data.py)
def netcdf():
    datadir = os.path.join(mlml_data(),'netcdf/')
    return datadir