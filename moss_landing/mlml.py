# -*- coding: utf-8 -*-

'''
Tools for working with data from MLML public data portal:
http://pubdata.mlml.calstate.edu
'''

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
try:
    import pandas as pd
except ImportError:
    pass
try:
    import xarray as xr
except ImportError:
    pass

def make_netcdf(station_dir,netcdf_dir,netcdf_prefix,station,download=False,overwrite=False):
    """
Create a netcdf file containing MLML historical seawater or weather data. The file will be created from csv and readme files already on disk, or they can be downloaded.

INPUT:
station_dir - string specifying the location of csv files (e.g. '/home/username/data/')
netcdf_dir - string specifying the location of netcdf files to be created (e.g. '/home/username/data/')
netcdf_prefix - string specifying filename pattern for netcdf files 
                the year will be appended this prefix
                (e.g. 'moss_landing_' for moss_landing_2015.nc, moss_landing_2016.nc, etc.)
station     - either 'seawater' or 'weather' (default: 'seawater')
download    - boolean specifying whether to download new files
              (default: False)
overwrite   - boolean specifying whether to overwrite the existing files, only used if downloading new data (default: False)
    """
    
    # download new data, if specified    
    if download == True:    
        download_station_data(station_dir,station,overwrite)
    
    # read data in csv files to xarray dataset
    d = read_csv_data(station_dir,format='dataset')
    
    # specify location of readme file and add metadata to dataset
    readme_file = station_dir + '1_README.TXT'
    d = add_metadata(d,station,readme_file)
    
    # Additional processing
    d = cleanup_raw(d)  
    d = add_flags(d)  
    
    d.attrs['history'] = d.attrs['history'] + 'netcdf file created using mlml.make_netcdf(station_dir'+station_dir+',netcdf_dir='+netcdf_dir+',netcdf_prefix='+netcdf_prefix+',station='+station+'download='+str(download)+',overwrite='+str(overwrite)+'): ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ')'
    
    # create netcdf files
    prefix = os.path.join(netcdf_dir,netcdf_prefix)
    years, datasets = zip(*d.groupby('time.year'))
    paths = [prefix+'%s.nc' % y for y in years]
    xr.save_mfdataset(datasets, paths)
        
def download_station_data(station_dir,station='seawater',overwrite=True):
    '''
Download all historical csv files for the MLML seawater intake or weather station. A latest version of the readme file is also downloaded. It is highly recommended to use different directories for seawater and weather, since the readme files have the same name. By default, new files are downloaded and existing files are overwritten.

INPUT:
station_dir - string specifying the local directory where you want to put 
              the data files
station     - either 'seawater' or 'weather' (default: 'seawater')
overwrite   - boolean specifying whether to overwrite the existing files 
              (default: False)
    
    '''
    # remote directories
    base_url = 'http://pubdata.mlml.calstate.edu/mlml_last/'    
    station_url = base_url + '/' + station + '/'

    # local directory
    station_dir = station_dir + '/'    
    
    # check whether a directory exists for this station
    if os.path.isdir(station_dir) == False:
        os.makedirs(station_dir)
    
    # find names of csv files that exist on the web and create a list
    # the csv filenames are in format yyyy-mm.csv
    urlpath =urlopen(station_url)
    html_string = urlpath.read().decode()
    urlpath.close()    
    file_pattern = '[0-9][0-9][0-9][0-9]-[0-9][0-9].csv'
    csv_list = re.findall(file_pattern,html_string)
    
    # get updated readme file
    urlr = urlopen(station_url + '1_README.TXT')
    fr = open(station_dir + '1_README.TXT','w')
    fr.write(str(urlr.read()))
    fr.close()
    urlr.close()
    
    # loop through each remote csv file and download if:
    # the file does not exist, the overwrite option is True or it is the last
    # file in the list (there may be new data in that file)
    for csv_name in csv_list:
        print('downloading ' + station + ': ' + csv_name)
        remote_file = station_url + csv_name
        local_file = station_dir + csv_name
        write_conditions = [os.path.exists(local_file) == False,
                            overwrite == True,
                            csv_name == csv_list[-1]]
        if any(write_conditions):
            urlfile = urlopen(remote_file)
            f = open(local_file,'w')
            filebytes = urlfile.read()
            f.write(filebytes.decode('utf8'))
            f.close()
            urlfile.close()
            
def read_csv_data(data_dir,format='dict'):
    '''
Read historical text data (.csv files) from the MLML seawater intake or weather station. The data must be stored locally, and can be downloaded automatically with the download_station_data() function.

Inputs:
data_dir - Specifies the directory where the data files are located. All files with the format yyyy-mm.csv in this directory will be read.

Options:
    format: output format
        format = 'dict' (default): dictionary
        format = 'dataframe': pandas DataFrame
        format = 'dataset': xarray DataSet

Output: dictionary, pandas DataFrame or xarray DataSet with keys/variable names taken from column headers
    '''
    
    file_list = glob(data_dir+'*.csv')

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
        floatcols = range(2,len(varnames))
        allcols = range(0,len(varnames))
        strcols = list(set(allcols)-set(floatcols))        
        
        # get numeric data, with missing values as NaN
        try:
            datamasked = np.genfromtxt(file_name,
                                 skip_header=1,
                                 delimiter=',',
                                 missing_values='-99999',
                                 usemask=True)
            data = datamasked.filled(np.nan)
 
            # get string data
            datastr = np.genfromtxt(file_name,
                             skip_header=1,
                             delimiter=',',
                             usecols=tuple(strcols),
                             dtype=str)

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
    
    if format is 'dataframe':
        # turn dictionary into pandas dataframe
        d = pd.DataFrame(d,index=dtime)
        d.index.name = 'time'
    elif format is 'dataset':
        # turn dictionary in xarray dataset, using dataframe as intermediate format
        d = pd.DataFrame(d,index=dtime)
        d.index.name = 'time'        
        d = xr.Dataset(d)
        d.attrs['history'] = 'dataset created using mlml.read_csv_data: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', ' 
    else:
        # default format: dictionary containing numpy arrays
        d['dtime'] = []    
        d['dtime'] = dtime
    
    return d
   
def list2date(datestr_list,fmt='%a %b %d %H:%M:%S %Y'):
    '''Convert a list of date strings to datetime format.

    INPUT:
    datestr_list: a list of strings that represent dates
    fmt: format of the date string, as would be input to strftime() or strptime()

    see https://docs.python.org/library/datetime.html#strftime-and-strptime-behavior

    OUTPUT:
    list of datetimes
    '''
    datetime_list = [datetime.strptime(datestr, fmt) for datestr in datestr_list]
    return datetime_list
    
def add_metadata(ds,station,readme_file):
    """
Add metadata to xarray dataset. Currently this adds lat and lon coordinates and puts the contents of the readme in an attribute. For the weather data, the anemometer height is also added as a coordinate.
    """    
    
    if station is 'seawater':
        ds.coords['lon'] = -121.7915
        ds.coords['lat'] = 36.8025
        
        ds['lon'].attrs['standard_name'] = 'longitude'
        ds['lon'].attrs['units'] = 'degrees_east'       
        
        ds['lat'].attrs['standard_name'] = 'latitude'
        ds['lat'].attrs['units'] = 'degrees_north'    
        
        ds['time'].attrs['standard_name'] = 'time'
        ds['time'].attrs['long_name'] = 'UTC time'
        
        ds.attrs['Conventions'] = 'CF-1.6'
        ds.attrs['title'] = 'Historical Seawater Data - Moss Landing Marine Labs'
        ds.attrs['institution'] = 'Moss Landing Marine Labs'
        ds.attrs['disclaimer'] = 'Moss Landing Marine Laboratories (MLML) provides these data "as is", with no warranty, expressed or implied, of the data quality or consistency. It is provided without support and without obligation on the part of MLML to assist in its use, correction, modification, or enhancement. For use in publication, authors should obtain written permission from the director of MLML, and acknowledge MLML as the data source in those publications.'
        ds.attrs['contacts'] = 'Data Manager: Jason Adelaars (jadelaars@mlml.calstate.edu), MLML Director: Jim Harvey (jharvey@mlml.calstate.edu)'
        ds.attrs['citation'] = 'Moss Landing Marine Labs. Historical Seawater Data (insert date range) [Internet]. [cited (insert date of download)]. Available from: pubdata.mlml.calstate.edu'
        ds.attrs['comment'] = 'Seawater data observations are collected from raw seawater drawn through an intake pipe. Prior to October 31, 2013, the intake opening was set at a depth of ~18.3 m (60ft) below mean lower low water (MLLW). Due to sediment build-up around the pipe, the intake opening was raised to ~16.6m (54.4ft) below MLLW on October 31, 2013. The intake opening is located at 36.8025N and 121.7915W\\r\\nThe seawater sensors are cleaned of biofouling agents on a weekly/twice-weekly interval, though some data drift can be observed in transmission, beam attenution, and fluorescence. Full maintenance Log available at: https://docs.google.com/a/mlml.calstate.edu/spreadsheet/ccc?key=0AnH2QAt-nlAldEFwY2ZyNUV6amJtLXllMDljdmM1SFE&pli=1#gid=0'
        
        flag_values = [1., 2., 3., 4., 9.]
        flag_meanings = 'pass not_evaluated suspect_or_of_high_interest fail missing_data'
        
        ds['ba'].attrs['long_name'] = 'beam attenuation'
        ds['ba'].attrs['standard_name'] = 'volume_beam_attenuation_coefficient_of_radiative_flux_in_sea_water'
        ds['ba'].attrs['units'] = 'm-1'
        ds['ba'].attrs['ancillary_variables'] = 'ba_flg'
        ds['ba'].attrs['comment'] = 'sensor: C-Star Transmissometer (10cm)'
        
        ds['ba_flg'].attrs['long_name'] = 'beam attenuation flag'
        ds['ba_flg'].attrs['standard_name'] = 'status_flag'
        ds['ba_flg'].attrs['flag_values'] = flag_values
        ds['ba_flg'].attrs['flag_meanings'] = flag_meanings
        ds['ba_flg'].attrs['comment'] =  'Sensor range test 0<ba<100, User range test 0<ba<50'
        
        ds['co2'].attrs['long_name'] = 'partial pressure of carbon dioxide'
        ds['co2'].attrs['standard_name'] = 'partial_pressure_of_carbon_dioxide_in_sea_water'
        ds['co2'].attrs['units'] = 'uatm'
        ds['co2'].attrs['ancillary_variables'] = 'co2_flg'
        ds['co2'].attrs['comment'] = 'sensor: Turner C-Sense'     
        
        ds['co2_flg'].attrs['long_name'] = 'partial pressure of carbon dioxide flag'
        ds['co2_flg'].attrs['standard_name'] = 'status_flag'
        ds['co2_flg'].attrs['flag_values'] = flag_values
        ds['co2_flg'].attrs['flag_meanings'] = flag_meanings      
        
        ds['cond'].attrs['long_name'] = 'conductivity'
        ds['cond'].attrs['standard_name'] = 'sea_water_electrical_conductivity'
        ds['cond'].attrs['units'] = 'S m-1'
        ds['cond'].attrs['ancillary_variables'] = 'cond_flg'
        ds['cond'].attrs['comment'] = 'sensor: Seabird SBE19 CTD'      
        
        ds['cond_flg'].attrs['long_name'] = 'conductivity flag'
        ds['cond_flg'].attrs['standard_name'] = 'status_flag'
        ds['cond_flg'].attrs['flag_values'] = flag_values
        ds['cond_flg'].attrs['flag_meanings'] = flag_meanings
        ds['cond_flg'].attrs['comment'] =  'Sensor range test 0<cond<7, User range test 0.75<cond<1.1, Spike test mean(cond_1.42hours)+/-2.4*stdev'
        
        ds['do2'].attrs['long_name'] = 'dissolved oxygen'
        ds['do2'].attrs['standard_name'] = 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
        ds['do2'].attrs['units'] = 'umol L-1'
        ds['do2'].attrs['ancillary_variables'] = 'do2_flg'
        ds['do2'].attrs['comment'] = 'Sensor: AADI Oxygen Optode'   
        
        ds['do2_flg'].attrs['long_name'] = 'dissolved oxygen flag'
        ds['do2_flg'].attrs['standard_name'] = 'status_flag'
        ds['do2_flg'].attrs['flag_values'] = flag_values
        ds['do2_flg'].attrs['flag_meanings'] = flag_meanings
        ds['do2_flg'].attrs['comment'] =  'Sensor range test 0<do2<500, User range test 0<do2<400, Spike test mean(do2_1.42hours)+/-2.4*stdev'
        
        ds['fluor'].attrs['long_name'] = 'chlorophyll'
        ds['fluor'].attrs['standard_name'] = 'chlorophyll_concentration_in_sea_water'
        ds['fluor'].attrs['units'] = 'ug L-1'
        ds['fluor'].attrs['ancillary_variables'] = 'fluor_flg'
        ds['fluor'].attrs['comment'] = 'Sensor: WetStar Fluorometer' 
        
        ds['fluor_flg'].attrs['long_name'] = 'chlorophyll flag'
        ds['fluor_flg'].attrs['standard_name'] = 'status_flag'
        ds['fluor_flg'].attrs['flag_values'] = flag_values
        ds['fluor_flg'].attrs['flag_meanings'] = flag_meanings
        ds['fluor_flg'].attrs['comment'] =  'Sensor range test 0.03<fluor<75, User range test 0<fluor<20, Spike test mean(fluor_1.42hours)+/-2.4*stdev'
        
        ds['osat'].attrs['long_name'] = 'dissolved oxygen percent saturation'
        ds['osat'].attrs['standard_name'] = 'fractional_saturation_of_oxygen_in_sea_water'
        ds['osat'].attrs['ancillary_variables'] = 'osat_flg'
        ds['osat'].attrs['units'] = 'percent'
        
        ds['osat_flg'].attrs['long_name'] = 'dissolved oxygen percent saturation flag'
        ds['osat_flg'].attrs['standard_name'] = 'status_flag'
        ds['osat_flg'].attrs['flag_values'] = flag_values
        ds['osat_flg'].attrs['flag_meanings'] = flag_meanings
        ds['osat_flg'].attrs['comment'] =  'Sensor range test 0<osat<100, User range test 0<osat<100, Spike test mean(osat_1.42hours)+/-2.4*stdev'
        
        ds['otemp'].attrs['long_name'] = 'optode temperature'
        ds['otemp'].attrs['standard_name'] = 'sea_water_temperature'
        ds['otemp'].attrs['units'] = 'degree_C'
        ds['osat'].attrs['ancillary_variables'] = 'otemp_flg'
        ds['otemp'].attrs['comment'] = 'Sensor: AADI Oxygen Optode'
        
        ds['otemp_flg'].attrs['long_name'] = 'optode temperature flag'
        ds['otemp_flg'].attrs['standard_name'] = 'status_flag'
        ds['otemp_flg'].attrs['flag_values'] = flag_values
        ds['otemp_flg'].attrs['flag_meanings'] = flag_meanings
        ds['otemp_flg'].attrs['comment'] =  'Sensor range test -5<temp<35, User range test 4<temp<22, Spike test mean(temp_1.42hours)+/-2.4*stdev'
        
        ds['ph'].attrs['long_name'] = 'pH total scale'
        ds['ph'].attrs['standard_name'] = 'sea_water_ph_reported_on_total_scale'
        ds['ph'].attrs['units'] = '1'
        ds['ph'].attrs['ancillary_variables'] = 'ph_flg'
        ds['ph'].attrs['comment'] = 'Sensor: Honeywell Durafet III'
        
        ds['ph_flg'].attrs['long_name'] = 'pH total scale flag'
        ds['ph_flg'].attrs['standard_name'] = 'status_flag'
        ds['ph_flg'].attrs['flag_values'] = flag_values
        ds['ph_flg'].attrs['flag_meanings'] = flag_meanings
        ds['ph_flg'].attrs['comment'] =  'Sensor range test 0<ph<14, User range test 7<ph<10, Spike test mean(ph_1.42hours)+/-2.4*stdev'

        ds['sal'].attrs['long_name'] = 'practical salinity'    
        ds['sal'].attrs['standard_name'] = 'sea_water_practical_salinity'
        ds['sal'].attrs['units'] = '1'
        ds['sal'].attrs['ancillary_variables'] = 'sal_flg'
        ds['sal'].attrs['comment'] = 'sensor: Seabird SBE19 CTD'
        
        ds['sal_flg'].attrs['long_name'] = 'practical salinity flag'
        ds['sal_flg'].attrs['standard_name'] = 'status_flag'
        ds['sal_flg'].attrs['flag_values'] = flag_values
        ds['sal_flg'].attrs['flag_meanings'] = flag_meanings
        ds['sal_flg'].attrs['comment'] =  'Sensor range test 20<sal<50, User range test 32<sal<35, Spike test mean(sal_1.42hours)+/-2.4*stdev'
        
        ds['temp'].attrs['long_name'] = 'temperature'
        ds['temp'].attrs['standard_name'] = 'sea_water_temperature'
        ds['temp'].attrs['units'] = 'degree_C'
        ds['temp'].attrs['ancillary_variables'] = 'temp_flg'
        ds['temp'].attrs['comment'] = 'sensor: Seabird SBE19 CTD'
        
        ds['temp_flg'].attrs['long_name'] = 'temperature flag'
        ds['temp_flg'].attrs['standard_name'] = 'status_flag'
        ds['temp_flg'].attrs['flag_values'] = flag_values
        ds['temp_flg'].attrs['flag_meanings'] = flag_meanings
        ds['temp_flg'].attrs['comment'] =  'Sensor range test -5<temp<35, User range test 4<temp<22, Spike test mean(temp_1.42hours)+/-2.4*stdev'
        
        ds['trans'].attrs['long_name'] = 'beam transmission'
        ds['trans'].attrs['units'] = 'percent'
        ds['trans'].attrs['ancillary_variables'] = 'trans_flg'
        ds['trans'].attrs['comment'] = 'sensor: C-Star Transmissometer (10cm)'     
        
        ds['trans_flg'].attrs['long_name'] = 'beam transmission flag'
        ds['trans_flg'].attrs['standard_name'] = 'status_flag'
        ds['trans_flg'].attrs['flag_values'] = flag_values
        ds['trans_flg'].attrs['flag_meanings'] = flag_meanings
        ds['trans_flg'].attrs['comment'] =  'Sensor range test 0<trans<100, User range test 0<trans<100'        

    elif station is 'weather':
        ds.coords['lon'] = -121.78842
        ds.coords['lat'] = 36.80040
        ds.coords['z'] = 3.3
        ds.coords['z'].attrs['name'] = 'anemometer height'
        ds.coords['z'].attrs['units'] = 'meters'
        
    with open(readme_file) as f:
        contents = f.read()
        ds.attrs['readme'] = contents
        
    ds.attrs['history'] = ds.attrs['history'] + 'attributes added to dataset using mlml.add_metadata: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', ' 
    
    return ds
    
def cleanup_raw(ds):
    '''
    Clean up variables by removing "raw" data
    '''

    # Loop through variables
    for varname, values in ds.data_vars.items():
        # If flag is NaN, data is missing
        if '_flg' in varname:
            ii, = np.where(np.isnan(values))
            ds[varname][ii] = 9.
        # Remove "raw" variables
        if '_raw' in varname:
            ds = ds.drop(varname)
    if 'phrawv' in ds.keys():
        ds = ds.drop('phrawv')
    # Remove tide variables
    if 'tide' in ds.keys():
        ds = ds.drop('tide')
    if 'tide_flg' in ds.keys():
        ds = ds.drop('tide_flg')
    # Remove extra time variables
    if 'utc_time' in ds.keys():
        ds = ds.drop('utc_time')
    if 'pst_time' in ds.keys():
        ds = ds.drop('pst_time')
    if 'unix_time' in ds.keys():
        ds = ds.drop('unix_time')
        
    ds.attrs['history'] = ds.attrs['history'] + 'raw data variables dropped using mlml.cleanup_raw: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', ' 
    
    return ds

def add_flags(ds):
    '''
    Flag questionable data based on users' subjective examination and automated criteria beyond initial QARTOD.
    '''
    
    # pH
    ph_ques_dates = [['2017-07-17 21:00','2017-07-19 19:30'],
                     ['2017-03-17 17:00','2017-05-17 19:30']]

    for drange in ph_ques_dates:
        ii, = np.where((ds['time']>=np.datetime64(drange[0]))
                       & (ds['time']<=np.datetime64(drange[1])))
        ds['ph_flg'][ii] = 3.

    # All variables
    var_list = ['ba','cond','do2','fluor','osat','otemp','ph','sal','temp','trans','co2']
    
    all_ques_dates = [['2017-05-03 18:13','2017-05-03 22:00'],
                      ['2017-06-07 19:00','2017-06-07 22:00'],
                      ['2017-07-03 06:20','2017-07-03 06:30'],
                      ['2017-07-18 18:00','2017-07-18 21:00']
                     ]
    for drange in all_ques_dates:
        ii, = np.where((ds['time']>=np.datetime64(drange[0]))
                       & (ds['time']<=np.datetime64(drange[1])))
        for var in var_list:
            ds[var+'_flg'][ii] = 3.

    # CO2
    co2_ques_dates = [['2015-01-01 00:00','2017-03-08 00:00'], # range set to 1000 ppm on first deployment (too low)
                      ['2017-07-18 20:01','2017-07-18 20:36'],
                      ['2017-07-12 21:00','2017-07-17 19:25'],
                      ['2017-07-05 23:00','2017-07-08 01:29'],
                      ['2017-03-22 17:06','2017-03-23 15:24']]

    for drange in co2_ques_dates:
        ii, = np.where((ds['time']>=np.datetime64(drange[0]))
                       & (ds['time']<=np.datetime64(drange[1])))
        ds['co2_flg'][ii] = 3.
        
    # salinity spikes    
    ii, = np.where(np.abs(np.diff(ds['sal'])) > 0.2)
    ii = ii
    ds['sal_flg'][ii] = 3.
    
    # isolated salinity points
    ii, = np.where(np.isnan(ds['sal'][1:-1]+ds['sal'][2:]) & np.isnan(ds['sal'][1:-1]+ds['sal'][:-2]))
    ii = ii + 1
    ds['sal_flg'][ii] = 3.
    
    # oxygen spikes
    ii, = np.where(np.abs(np.diff(ds['do2'])) > 100)
    ii = ii
    ds['do2_flg'][ii] = 3.
    
    # isolated oxygen points(or pairs of points)
    ii, = np.where(np.isnan(ds['do2'][1:-1]+ds['do2'][2:]) & np.isnan(ds['do2'][1:-1]+ds['do2'][:-2]))
    ii = ii + 1
    ds['do2_flg'][ii] = 3.
    ii, = np.where(np.isnan(ds['do2'][1:-1]+ds['do2'][2:]) & np.isnan(ds['do2'][1:-1]+ds['do2'][:-2]))
    ii = ii + 1
    ds['do2_flg'][ii] = 3.
    
    ds.attrs['history'] = ds.attrs['history'] + 'Flags added using mlml.add_flags: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ', ' 
    
    return ds