import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.dates import date2num
import pandas as pd
import xarray as xr
import sys
import os

#def apply_manual_flags(ds):


def manual_flag(ds):
    '''
Specify time ranges to flag or unflag with a GUI.

The date ranges are stored in a text file in manual_flag_output/
The file name includes:
- the variable
- flag mode (flag or unflag)
- the date and time on which the file was created
    '''

    var_list = list(ds.variables)
    print(var_list)

    var_plt = 'do2'
    var_in = input('Choose a variable to plot (default: do2):')
    if var_in != '':
        var_plt = var_in
        if var_plt not in var_list:
            raise('variable not understood')

    var = var_plt
    var_in = input('Choose a variable to apply the flags to, or type "all" (default: plot variable):')
    if var_in != '':
        var = var_in

    flag_mode = 'flag'
    fl_in = input('Flag mode: Choose whether to flag or unflag (default: flag, type "u" to unflag):')
    if fl_in != '':
        if fl_in == "u":
            flag_mode = 'unflag'
        else:
            raise('flag mode not understood')

    good = ds[var+'_flg'] == 1

    #dnum = date2num(np.array(ds['time']))
    unix_time = (np.array(ds['time'])-np.datetime64('1970-01-01'))/np.timedelta64(1,'s')

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(1, 1, 1)
    if flag_mode == 'unflag':
        ax.plot(unix_time,ds[var],'k-',lw=0.5)
    ax.plot(unix_time[good],ds[var][good],'c-',lw=0.5)
    ax.plot(unix_time[good],ds[var][good],'bo',ms=1)
    xl = plt.xlim()
    yl = plt.ylim()

    now = pd.Timestamp.now()
    datetime_str = (str(now.year)+'-'+
                    str(now.month).zfill(2)+'-'
                    +str(now.day).zfill(2)+'-'
                    +str(now.hour).zfill(2)+
                    str(now.minute).zfill(2)+
                    str(now.second).zfill(2))
    print(datetime_str)

    out_dir = 'manual_flag_output/'
    out_file = var + '_' + flag_mode + '_' + datetime_str + '.txt'
    out_path = os.path.join(out_dir,out_file)

    f = open(out_path,'w')
    f.write('# date range of values to ' + flag_mode +
            ', created using manual_flag() on ' + datetime_str)
    f.close()

    done_flagging = False

    while not done_flagging:
        cursor = Cursor(ax, useblit=True, color='k', linewidth=1)
        zoom_ok = False

        plt.xlim(xl)
        plt.ylim(yl)

        print('\nZoom or pan to view, \npress spacebar when ready to click:\n')

        while not zoom_ok:
            zoom_ok = plt.waitforbuttonpress()

        print('Click twice to select time range ')
        val = plt.ginput(2)
        print('Selected values: ', val)

        ti1 = np.argmin(np.abs(unix_time - val[0][0]))
        ti2 = np.argmin(np.abs(unix_time - val[1][0]))

        t1 = ds['time'][ti1].values
        t2 = ds['time'][ti2].values

        f = open(out_path,'a')
        f.write(str(t1)+','+str(t2))
        f.close()

        if flag_mode == 'flag':
            ii, = np.where((ds['time'][good] >= t1) & (ds['time'][good] <= t2))
            ax.plot(unix_time[good][ii],ds[var][good][ii],'r.')
        if flag_mode == 'unflag':
            ii, = np.where((ds['time'] >= t1) & (ds['time'] <= t2))
            ax.plot(unix_time[ii],ds[var][ii],'r.')

        xl = plt.xlim()
        yl = plt.ylim()

        response = input('Press 1 if done flagging, 2 to abort, return to continue')

        if len(response) > 0:
            if response[-1] == '1':
                done_flagging = True
                break
            if response[-1] == '2':
                done_flagging = True
                os.remove(out_path)
                print(out_file+' deleted from '+out_path)
                break

if __name__ == '__main__':
    sys.path.append('..')
    import mlml_data_path

    # load data
    nc_dir = mlml_data_path.moss_landing() + '/netcdf/'
    ds = xr.open_mfdataset(nc_dir+'*.nc')
    ds.load()

    manual_flag(ds)
