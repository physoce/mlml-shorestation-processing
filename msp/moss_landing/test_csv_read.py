try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
import re
import os
import sys
import numpy as np
from glob import glob
from datetime import datetime, timedelta
from manual_flag import load_manual_flags
try:
    import pandas as pd
except ImportError:
    pass
try:
    import xarray as xr
except ImportError:
    pass
from mlml import list2date

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

data_dir = sw_csv_dir
# def read_csv_data(data_dir,format='dict'):
#     '''
# Read historical text data (.csv files) from the MLML seawater intake or weather station. The data must be stored locally, and can be downloaded automatically with the download_station_data() function.

# Inputs:
# data_dir - Specifies the directory where the data files are located. All files with the format yyyy-mm.csv in this directory will be read.

# Options:
#     format: output format
#         format = 'dict' (default): dictionary
#         format = 'dataframe': pandas DataFrame
#         format = 'dataset': xarray DataSet

# Output: dictionary, pandas DataFrame or xarray DataSet with keys/variable names taken from column headers
#     '''

file_list = sorted(glob(data_dir+'*.csv'))

for fi,file_name in enumerate(file_list):
    print('reading ' + file_name)

    # get list of variable names from header of first file
    f = open(file_name,'r')
    header = f.readline()
    f.close()
    header = header.strip('\r\n')
    varnames = header.split(',')

    #if first file, initialize dictionary with key and empty list for each variable
    if fi == 0:
        d = dict()
        for var in varnames:
            d[var] = []

    # specify which columns contain numeric data
    if varnames[0] == '': # extra first column in some files
        floatcols = range(3,len(varnames))
        allcols = range(1,len(varnames))
        strcols = list(set(allcols)-set(floatcols))
    else:
        floatcols = range(2,len(varnames))
        allcols = range(0,len(varnames))
        strcols = list(set(allcols)-set(floatcols))

    # get numeric data, with missing values as NaN
    try:
        # datamasked = np.genfromtxt(file_name,
        #                         skip_header=1,
        #                         delimiter=',',
        #                         missing_values='-99999',
        #                         usemask=True)
        # data = datamasked.filled(np.nan)
        datatmp = pd.read_csv(file_name,
                    na_values=-99999).values
        # filter lines where newline character in suna string puts garbage on following row
        badi, = np.where(datatmp[:, 0] == '"') 
        data = pd.read_csv(file_name,
                    na_values=-99999,
                    skiprows=badi).values

        # get string data
        datastrtmp = np.genfromtxt(file_name,
                            skip_header=1,
                            delimiter=',',
                            usecols=tuple(strcols),
                            dtype=str)
        badi2, = np.where(datastrtmp[:, 0] == '"') 
        datastr = pd.read_csv(file_name,
                    na_values=-99999,
                    usecols=tuple(strcols),
                    skiprows=badi).values

        # append data variables
        if data.size != 0:
            for si,col in enumerate(strcols):
                vname = varnames[col]
                d[vname] = np.append(d[vname],datastr[:,si])
            for col in floatcols:
                vname = varnames[col]
                if vname not in d.keys():
                    d[vname] = np.nan*np.zeros(np.shape(d[varnames[col-1]])) # create new variable (all NaN's)
                else:
                    d[vname] = np.append(d[vname],data[:,col]) # append data
            for key,value in d.items():
                if key not in varnames:
                    d[key] = np.append(d[key],np.nan*np.zeros(np.shape(data[:,2]))) # append NaN's

    except:
        print('could not read '+file_name)
        pass

# create date variables
# put in a numpy array for easy indexing
# new variable for datetime
dtime = np.array(list2date(d['utc_time'],'%Y-%m-%dT%H:%M:%SZ'))

# remove duplicate times
ii = np.where(np.diff(dtime) > timedelta(0.))[0]
dtime = dtime[ii]
for key,value in d.items():
    d[key] = d[key][ii]

# Try loading in pandas or xarray format if specified, default to dictionary format
if format == 'dataset':
    if 'xarray' not in sys.modules:
        format = 'dataframe'
        print("Warning: xarray not installed, loading MLML data in pandas dataframe format instead")
if format == 'dataframe':
    if 'pandas' not in sys.modules:
        format = 'dict'
        print("Warning: pandas not installed, loading MLML data in dictionary format instead")

if format == 'dataframe':
    # turn dictionary into pandas dataframe
    d = pd.DataFrame(d,index=dtime)
    d.index.name = 'time'
elif format == 'dataset':
    # turn dictionary in xarray dataset, using dataframe as intermediate format
    d = pd.DataFrame(d,index=dtime)
    d.index.name = 'time'
    d = xr.Dataset(d)
    d.attrs['history'] = 'dataset created using mlml.read_csv_data: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', '
else:
    # default format: dictionary containing numpy arrays
    d['dtime'] = []
    d['dtime'] = dtime

