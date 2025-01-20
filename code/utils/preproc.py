import numpy as np
import xarray as xr

## paths
# source data folder
cesm_path = '/net/meso/climphys/cesm212/b.e212.B1850cmip6.f09_g17.001/archive/atm/hist/'
boost_path = "/net/meso/climphys/cesm212/boosting/archive/"
# processed data
output_path = '/net/xenon/climphys/lbloin/boost_proba/'

#defining the area to crop
min_lon = -123+360
min_lat = 45
max_lon = -119+360
max_lat = 52
KtoC = -273.15


# preprocessing function used to crop the area and compute 5-day running means over the area
def preprocess(ds,mean=True):
    da = ds["TREFHTMX"]
    if mean == True:
        da = weighted_mean(da)
    ds_out = (da.rolling(time=5, center=True).mean()+KtoC).to_dataset(name="Tx5d")
    return ds_out

def weighted_mean(da):
    da_slice = da.sel(lat=slice(min_lat,max_lat),lon=slice(min_lon,max_lon))
    wgt = np.cos(np.deg2rad(da.lat))
    da_slice = da_slice.weighted(wgt).mean(("lat","lon"))
    return da_slice