# -*- coding: utf-8 -*-

"""
edit this file to define paths to MLML data

example:
-------

import mlml_data_path

data_dir = mlml_data_path.seawater()
"""

import os

# top-level data directory
def mlml_data():
    datadir = '~/work/Data/MLML/'
    return datadir
    
# csv files
def seawater():
    datadir = os.path.join(mlml_data(),'seawater/')
    return datadir

def weather():
    datadir = os.path.join(mlml_data(),'weather/')
    return datadir

# NetCDF files (created from csv files in get_mlml_data.py)
def netcdf():
    datadir = os.path.join(mlml_data(),'netcdf/')
    return datadir