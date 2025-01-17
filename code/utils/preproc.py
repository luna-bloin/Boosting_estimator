import numpy as np
import xarray as xr


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
    ds_out = (da.rolling(time=5, center=True).mean()+KtoC).rename({"TREFHTMX":"Tx5d","Z500":"Z500x5d"})
    return ds_out

def weighted_mean(da):
    da_slice = da.sel(lat=slice(min_lat,max_lat),lon=slice(min_lon,max_lon))
    wgt = np.cos(np.deg2rad(da.lat))
    da_slice = da_slice.weighted(wgt).mean(("lat","lon"))
    return da_slice