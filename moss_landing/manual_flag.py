import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor
import matplotlib.dates as mdates
import pandas as pd
import xarray as xr
from glob import glob
import sys
import os

def load_manual_flags(ds,var):
    '''
Reads text files that store dates of data to be flagged or unflagged.
Returns time index of data to be flagged or unflagged in dataset.
This function does not modify the dataset or return a new dataset.

Inputs:
ds - xarray dataset
var - variable to apply flags

Outputs:
flagi - time index of data to be flagged (flag = 3)
unflagi - time index of data to be unflagged (flag = 1)
    '''

    out_dir = 'manual_flag_output/'

    flist = glob(out_dir+'/*.txt')
    flist_sorted = sorted(flist, key=os.path.getmtime)

    flagi = np.zeros(len(ds['time']),dtype=bool)
    unflagi = np.zeros(len(ds['time']),dtype=bool)

    if len(flist_sorted) > 0 :
        for fpath in flist_sorted:
            fname = fpath.rsplit('/')[-1]
            fvar = fname.split('_',2)[0]
            flag_mode = fname.split('_',2)[1]

            if var == fvar:
                df_flag = pd.read_csv(fpath,names=['t1','t2'],skiprows=1)
                for i in range(df_flag.shape[0]):
                    ti = ((ds['time'] >= np.datetime64(df_flag['t1'][i])) &
                          (ds['time'] <= np.datetime64(df_flag['t2'][i])))
                    if flag_mode == 'flag':
                        flagi[ti] = True
                    elif flag_mode == 'unflag':
                        unflagi[ti] = True

    return flagi, unflagi

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

    dnum = (np.array(ds['time'])-np.datetime64('1970-01-01'))/np.timedelta64(1,'D')

    flagi,unflagi = load_manual_flags(ds,var)

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(1, 1, 1)
    if flag_mode == 'unflag':
        ax.plot(dnum,ds[var],'k-',lw=0.5)
    ax.plot(dnum[good],ds[var][good],'c-',lw=0.5)
    ax.plot(dnum[good],ds[var][good],'bo',ms=1)
    ax.plot(dnum[flagi],ds[var][flagi],'r.')
    ax.plot(dnum[unflagi],ds[var][unflagi],'c.')
    xl = plt.xlim()
    yl = plt.ylim()

    formatter = mdates.DateFormatter("%Y-%m-%d")
    ax.xaxis.set_major_formatter(formatter)

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
            ', created using manual_flag() on '
            + datetime_str + '\n')
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

        ti1 = np.argmin(np.abs(dnum - val[0][0]))
        ti2 = np.argmin(np.abs(dnum - val[1][0]))

        t1 = ds['time'][ti1].values
        t2 = ds['time'][ti2].values

        f = open(out_path,'a')
        f.write(str(t1) + ',' + str(t2) + '\n')
        f.close()

        if flag_mode == 'flag':
            ii, = np.where((ds['time'][good] >= t1) & (ds['time'][good] <= t2))
            ax.plot(dnum[good][ii],ds[var][good][ii],'r.')
        if flag_mode == 'unflag':
            ii, = np.where((ds['time'] >= t1) & (ds['time'] <= t2))
            ax.plot(dnum[ii],ds[var][ii],'c.')

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
    #load_manual_flags(ds)
