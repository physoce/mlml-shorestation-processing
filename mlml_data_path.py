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
    datadir = os.path.join(mlml_data(),'moss_landing/')
    return datadir

def moss_landing_weather():
    datadir = os.path.join(mlml_data(),'moss_landing_weather/')
    return datadir

def monterey_wharf():
    datadir = os.path.join(mlml_data(),'monterey_wharf/')
    return datadir