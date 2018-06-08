# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

sys.path.append('..')
import mlml_data_path
import qartod

def csv_to_dataset(sw_csv):
    '''
    Load Monterey Wharf csv file into an xarray dataset.
    '''

    # load csv into pandas dataframe
    df = pd.read_csv(sw_csv,parse_dates=['utc_time'])
    df = df.set_index('utc_time')

    # some missing values are reported as -99.99
    df = df.replace(-99.99, np.nan)
    
    # convert to xarray dataset
    ds = xr.Dataset.from_dataframe(df)
    
    # rename utc_time to standard name
    ds = ds.rename({'utc_time':'time'})

    ds = ds.drop('Index')

    # make sure variables do not have spaces in names
    ds = ds.rename({'BGA PE_cells':'BGA_PE_cells','ODO Conc':'ODO_Conc'})
    
    ds.attrs['history'] = 'Data read from ' + sw_csv.split('/')[-1] + ' using wharf.csv_to_dataset: ' + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    return ds

def run_qartod(ds):
    '''
    Run selected QARTOD functions on Wharf dataset for automated QC
    '''

    
    # specify variables and qartod parameters
    
    var_list = ['Temp','Salinity','Cond','ODO_Conc','pH','Turbidity','Chlorophyll_rfu','Chlorophyll_ug',
               'BGA_PE_cells','PE_rfu']
    
    sensor_range = dict()
    sensor_range['Temp'] = [-5,50]
    sensor_range['Salinity'] = [0,70]
    sensor_range['Cond'] = [0,100]      # mS/cm
    sensor_range['ODO_Conc'] = [0,50]   # mg/L
    sensor_range['pH'] = [0,14]   
    sensor_range['Turbidity'] = [0,1000]   # NTU
    sensor_range['Chlorophyll_rfu'] = [0,100]
    sensor_range['Chlorophyll_ug'] = [0,400]
    sensor_range['BGA_PE_cells'] = [0,200000]
    sensor_range['PE_rfu'] = [0,200000]
    
    user_range = dict()
    user_range['Temp'] = [4,25]
    user_range['Salinity'] = [30,35]
    user_range['Cond'] = [20,60]      # mS/cm
    user_range['ODO_Conc'] = [0,25]   # mg/L
    user_range['pH'] = [7,10]   
    
    # vectorize qartod functions
    userRange = np.vectorize(qartod.userRange,excluded=[1,2])
    sensorRange = np.vectorize(qartod.sensorRange,excluded=[1,2])
    
    # for the metadata
    flag_values = [1., 2., 3., 4., 9.]
    flag_meanings = 'pass not_evaluated suspect_or_of_high_interest fail missing_data'    
    
    for var in var_list:
        sensor_range_flg = sensorRange(ds[var],sensor_range[var][0],sensor_range[var][1])
        comment = 'sensor range test: ['+str(sensor_range[var][0])+','+str(sensor_range[var][1])+']'
        
        try:
            user_range_flg = userRange(ds[var],sensor_range[var][0],sensor_range[var][1])
            comment = comment + ', user range test: ['+str(user_range[var][0])+','+str(user_range[var][1])+']'
        except:
            user_range_flg = np.nan*np.ones(np.shape(ds[var]))
            
        all_flg = np.stack([sensor_range_flg,user_range_flg])
        max_flg = np.nanmax(all_flg,axis=0).astype(float)
        
        ds[var+'_flg'] = xr.DataArray(max_flg,coords=[ds['time']])
        ds[var+'_flg'].attrs['standard_name'] = 'status_flag'
        ds[var+'_flg'].attrs['flag_values'] = flag_values
        ds[var+'_flg'].attrs['flag_meanings'] = flag_meanings
        ds[var+'_flg'].attrs['comment'] = comment
        
    ds.attrs['history'] = ds.attrs['history'] + ', automated QC using run_qartod:' + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    return ds

def manual_flags(ds):
    
    do_ques_dates = [['2017-04-17 00:00','2018-06-08 00:00']] # extremely low values
    
    for drange in do_ques_dates:
        ii, = np.where((ds['time']>=np.datetime64(drange[0]))
                       & (ds['time']<=np.datetime64(drange[1]))
                       & np.isfinite(ds['ODO_Conc']))
        ds['ODO_Conc_flg'][ii] = 3.
    
    return ds
    
def add_metadata(ds):
    # add metadata
    
    ds.attrs['Conventions'] = 'CF-1.6'
    ds.attrs['title'] = 'Historical Shore Station Data - Monterey Wharf'
    ds.attrs['institution'] = 'Moss Landing Marine Labs'
    ds.attrs['disclaimer'] = 'Moss Landing Marine Laboratories (MLML) provides these data "as is", with no warranty, expressed or implied, of the data quality or consistency. It is provided without support and without obligation on the part of MLML to assist in its use, correction, modification, or enhancement. For use in publication, authors should contact the principal investigators, and acknowledge MLML as the data source in those publications.'
    ds.attrs['contacts'] = 'Principal Investigators: Jason Smith (jsmith@mlml.calstate.edu), Tom Connolly (tconnolly@mlml.calstate.edu), Data Manager: Jason Adelaars (jadelaars@mlml.calstate.edu), MLML Director: Jim Harvey (jharvey@mlml.calstate.edu)'
    ds.attrs['citation'] = 'Moss Landing Marine Labs. Monterey Wharf Shore Station Data (insert date range) [Internet]. [cited (insert date of download)]. Available from: pubdata.mlml.calstate.edu'
    ds.attrs['comment'] = 'Observations are collected with a YSI 6600V2-4 sonde, fixed to a pier.'

    ds.coords['lon'] = -121.8894
    ds.coords['lat'] = 36.6050

    ds['lon'].attrs['standard_name'] = 'longitude'
    ds['lon'].attrs['units'] = 'degrees_east'       

    ds['lat'].attrs['standard_name'] = 'latitude'
    ds['lat'].attrs['units'] = 'degrees_north'      
    
    ds['time'].attrs['standard_name'] = 'time'
    ds['time'].attrs['long_name'] = 'UTC time'
    
    ds = ds.drop(['pst_time','unix_time'])
    
    ds['BGA_PE_cells'].attrs['long_name'] = 'number concentration of blue-green algae cells estimated from Phycoerythrin fluorescence'
    ds['BGA_PE_cells'].attrs['units'] = 'mL-1'
    ds['BGA_PE_cells'].attrs['ancillary_variables'] = 'BGA_PE_cells_flg'    
    
    ds['PE_rfu'].attrs['long_name'] = 'Phycoerythrin relative fluorescence'
    ds['PE_rfu'].attrs['units'] = 'percent'
    ds['PE_rfu'].attrs['ancillary_variables'] = 'PE_rfu_flg'        
    
    ds['Battery'].attrs['long_name'] = 'battery voltage'
    ds['Battery'].attrs['units'] = 'volts'
    
    ds['Battery2'].attrs['long_name'] = 'battery voltage'
    ds['Battery2'].attrs['units'] = 'volts'    
    
    ds['Chlorophyll_rfu'].attrs['long_name'] = 'chlorophyll, relative fluorescence units'
    ds['Chlorophyll_rfu'].attrs['units'] = 'percent'
    ds['Chlorophyll_rfu'].attrs['ancillary_variables'] = 'Chlorophyll_rfu_flg'     
    
    ds['Chlorophyll_ug'].attrs['long_name'] = 'chlorophyll'
    ds['Chlorophyll_ug'].attrs['standard_name'] = 'chlorophyll_concentration_in_sea_water'
    ds['Chlorophyll_ug'].attrs['units'] = 'ug L-1'
    ds['Chlorophyll_ug'].attrs['ancillary_variables'] = 'Chlorophyll_ug_flg'
    
    ds['Cond'].attrs['long_name'] = 'conductivity'
    ds['Cond'].attrs['standard_name'] = 'sea_water_electrical_conductivity'
    ds['Cond'].attrs['units'] = 'mS cm-1'
    ds['Cond'].attrs['ancillary_variables'] = 'Cond_flg'    
    
    ds['Depth'].attrs['long_name'] = 'depth'
    ds['Depth'].attrs['standard_name'] = 'depth'
    ds['Depth'].attrs['units'] = 'm'     
    
    ds['ODO_Conc'].attrs['standard_name'] = 'mass_concentration_of_oxygen_in_sea_water' 
    ds['ODO_Conc'].attrs['long_name'] = 'dissolved oxygen'
    ds['ODO_Conc'].attrs['units'] = 'mg L-1'
    ds['ODO_Conc'].attrs['ancillary_variables'] = 'ODO_Conc_flg'    

    ds['Salinity'].attrs['long_name'] = 'practical salinity'    
    ds['Salinity'].attrs['standard_name'] = 'sea_water_practical_salinity'
    ds['Salinity'].attrs['units'] = '1'
    ds['Salinity'].attrs['ancillary_variables'] = 'Salinity_flg'
    
    ds['Temp'].attrs['long_name'] = 'temperature'
    ds['Temp'].attrs['standard_name'] = 'sea_water_temperature'
    ds['Temp'].attrs['units'] = 'degree_C'
    ds['Temp'].attrs['ancillary_variables'] = 'Temp_flg'
    
    ds['Turbidity'].attrs['long_name'] = 'turbidity'
    ds['Turbidity'].attrs['standard_name'] = 'sea_water_turbidity'
    ds['Turbidity'].attrs['units'] = '1'
    ds['Turbidity'].attrs['ancillary_variables'] = 'Turbidity_flg' 
    ds['Turbidity'].attrs['comment'] = 'Dimensionless quantity expressed in NTU (Nephelometric Turbidity Units).'
    
    
    ds['pH'].attrs['long_name'] = 'pH total scale'
    ds['pH'].attrs['standard_name'] = 'sea_water_ph_reported_on_total_scale'
    ds['pH'].attrs['units'] = '1'
    ds['pH'].attrs['ancillary_variables'] = 'pH_flg'
    
    ds['pHmV'].attrs['long_name'] = 'pH sensor voltage'
    ds['pHmV'].attrs['units'] = 'mV'
    
    ds.attrs['history'] = ds.attrs['history'] + ', metadata attributes added using add_metadata:' + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    return ds

def write_netcdf(ds,netcdf_dir,netcdf_prefix):
    
    ds.attrs['history'] = ds.attrs['history'] + ', written to NetCDF files using write_netcdf:' + datetime.now().strftime("%Y-%m-%d %H:%M:%S")
   
    prefix = os.path.join(netcdf_dir,netcdf_prefix)
    years, datasets = zip(*ds.groupby('time.year'))
    paths = [prefix+'%s.nc' % y for y in years]
    xr.save_mfdataset(datasets, paths)
    
if __name__ == '__main__':
    
    # get paths of data files
    # paths to data should be modified in mlml_data_path.py
    sw_dir = mlml_data_path.monterey_wharf()

    # path to csv file
    sw_csv = os.path.join(sw_dir,'csv/WharfData_all.csv')

    # path to final netcdf files
    sw_nc_dir = os.path.join(sw_dir,'netcdf/')
    nc_prefix = 'monterey_wharf_'
    
    ds = csv_to_dataset(sw_csv)
    ds = run_qartod(ds)
    ds = manual_flags(ds)
    ds = add_metadata(ds)
    write_netcdf(ds,sw_nc_dir,nc_prefix)
    
    print('done creating Monterey Wharf dataset')