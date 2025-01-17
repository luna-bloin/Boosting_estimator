import numpy as np
import xarray as xr
import glob
import sys
sys.path.append("../utils")
import preproc as pc
## paths
# source data folder
cesm_path = '/net/meso/climphys/cesm212/b.e212.B1850cmip6.f09_g17.001/archive/atm/hist/'
# processed data
output_path = '/net/xenon/climphys/lbloin/boost_proba/'

# === Preprocessing all PI control years 1-4003 ===
print("preprocessing PI control run")
#selecting all existing PI runs
h1_filenames = sorted(glob.glob(f"{cesm_path}b.e212.B1850cmip6.f09_g17.001.cam.h1.*.nc"))
#open and preprocess
ds_PI = xr.open_mfdataset(h1_filenames,combine="by_coords",parallel=True, preprocess=pc.preprocess)
print("opened")
# calculate climatology and anomalies
clim_PI = ds_PI.groupby("time.dayofyear").mean().Tx5d
ds_PI["Tx5d_anom"] = ds_PI.Tx5d.groupby("time.dayofyear")-clim_PI
#save
ds_PI.to_netcdf(f'{output_path}PNW_PI_control_simulation.nc')
print("saved")
# === Saving the short window time periods as separate files
# define reference periods
time_periods = {
    'T1':slice("1801","1850"),
    'T2':slice("1851","1900"),
    'T3':slice("1801","1900"),
}
print("saving sub-periods of PI control")
anomalies_tim = {}
for tim in time_periods:
    ds_tim = ds_PI.sel(time=time_periods[tim]) #select period
    clim_tim = ds_tim.Tx5d.rolling(time=30,center=True).mean().groupby("time.dayofyear").mean() # climatology for sub period
    anomalies_tim[tim] = clim_tim
    #anomalies
    ds_tim["Tx5d_anom"] = ds_tim.Tx5d.groupby("time.dayofyear")-clim_tim
    #save
    ds_tim.to_netcdf(f'{output_path}PNW_PI_test_slice_{tim}.nc')


    