# MLML shore station processing

Tom Connolly and Jason Adelaars, Moss Landing Marine Labs

Contact: tconnolly@mlml.calstate.edu

Data processing for CeNCOOS shore stations operated by Moss Landing Marine Labs.

This software takes original text files as input, adds metadata, and creates CF-compliant NetCDF files. Quality control flags are automatically applied where necessary, and additional questionable data are flagged manually.

### Data access

TBD

### Shore station information

Moss Landing Marine Labs Seawater Intake
* https://www.cencoos.org/data/shore/mosslanding
* http://pubdata.mlml.calstate.edu/seawater

Monterey Municipal Wharf 2
* https://www.cencoos.org/data/shore/monterey

### Software requirements

* pandas https://pandas.pydata.org
* xarray http://xarray.pydata.org

### Running the software

1. Edit [mlml_data_path.py](mlml_data_path.py) to specify the paths to input and output file directories.

The default directory structure (located in top-level directory `~/work/Data/` by default) is:
```
MLML_shore_stations/
  monterey_wharf/
    csv/          
    netcdf/       
  moss_landing/
    csv/          
    netcdf/       
```
Input files are located in `csv` directories. Output files are create in `netcdf` directories.

2. To process Monterey Wharf data, run [process_monterey_wharf.py](moss_landing/process_monterey_wharf.py). This requires the input data file `WharfData_all.csv`

3. To process Moss Landing seawater intake data, run [process_mlml_data.py](moss_landing/process_mlml_data.py). This file can be edited to automatically download the csv files from the MLML public data server at http://pubdata.mlml.calstate.edu/seawater
