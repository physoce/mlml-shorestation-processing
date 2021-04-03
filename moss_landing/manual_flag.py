import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor
import pandas as pd
import xarray as xr
import sys
import os

def manual_flag(ds):
    '''
Specify time ranges to flag as suspect with a GUI.

The suspect date ranges will be stored in a text file in manual_flag_output/
    '''

    var_list = list(ds.variables)
    print(var_list)

    var = 'do2'
    var_in = input('Choose a variable to plot (default: do2):')
    if var_in != '':
        var = var_in

    good = ds[var+'_flg'] == 1

    unix_time = (np.array(ds['time'])-np.datetime64('1970-01-01'))/np.timedelta64(1,'s')

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(unix_time[good],ds[var][good])

    now = pd.Timestamp.now()
    datetime_str = (str(now.year)+'-'+
                    str(now.month).zfill(2)+'-'
                    +str(now.day).zfill(2)+'-'
                    +str(now.hour).zfill(2)+
                    str(now.minute).zfill(2)+
                    str(now.second).zfill(2))
    print(datetime_str)

    out_dir = 'manual_flag_output/'
    out_file = var + '_flag_' + datetime_str + '.txt'
    out_path = os.path.join(out_dir,out_file)

    f = open(out_path,'w')
    f.write('# date range of questionable values, created using manual_flag() on '+datetime_str)
    f.close()

    done_flagging = False

    while not done_flagging:
        cursor = Cursor(ax, useblit=True, color='k', linewidth=1)
        zoom_ok = False

        print('\nZoom or pan to view, \npress spacebar when ready to click:\n')

        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress()

        print('Click twice to select time range ')
        val = plt.ginput(2)
        print('Selected values: ', val)

        ti1 = np.argmin(np.abs(unix_time - val[0][0]))
        ti2 = np.argmin(np.abs(unix_time - val[1][0]))

        t1 = ds['time'][ti1]
        t2 = ds['time'][ti2]

        f = open(out_path,'a')
        f.write(str(t1)+','+str(t2))
        f.close()

        response = input('Press 1 if done flagging, any other key to continue')

        if len(response) > 0:
            if response[-1] == '1':
                done_flagging = True
                break



if __name__ == '__main__':
    sys.path.append('..')
    import mlml_data_path

    # load data
    nc_dir = mlml_data_path.moss_landing() + '/netcdf/'
    ds = xr.open_mfdataset(nc_dir+'*.nc')
    ds.load()

    manual_flag(ds)
