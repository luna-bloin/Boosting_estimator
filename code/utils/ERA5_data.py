import xarray as xr
import preproc as pc

path = "/net/xenon/climphys/lbloin/boost_proba/"
def open_era5(var="T2M"):
    if var == "T2M_clim":
        detrend = "_detrended_clim"
    else:
        detrend = "" 
    if var[0:3] == "T2M":
        var_file = "temperature"
    elif var == "Z":
        var_file = "z500"
    ds = xr.open_dataset(f'{path}{var_file}{detrend}_ERA5_2021.nc')[var]
    #resample to daily and do rolling mean
    ds = ds.resample(time="1D").max().rolling(time=5,center=True).mean()
    return ds

def get_era5_2021_PNW_temperature_anomaly():
    # open preprocessed June-July ERA5 temperature data
    temp = open_era5()
    # open detrend June-July ERA5 temperature climatology
    temp_clim_detrend = open_era5(var="T2M_clim")
    temp_anom_detrend_PNW = pc.weighted_mean(temp-temp_clim_detrend) 
    # get max anom and date of max anom
    PNW_heatw_anomaly = temp_anom_detrend_PNW.max()
    PNW_heatw_anomaly_day = temp_anom_detrend_PNW.idxmax()
    return PNW_heatw_anomaly, PNW_heatw_anomaly_day

def get_era5_map_at_peak():
    # open preprocessed June-July ERA5 temperature data
    temp = open_era5()
    # open detrend June-July ERA5 temperature climatology
    temp_clim_detrend = open_era5(var="T2M_clim")
    # open preprocessed June-July ERA5 Z500 data
    z500 = open_era5(var="Z")
    # max anom day 
    PNW_heatw_anomaly_day = get_era5_2021_PNW_temperature_anomaly()[1]
    # get maps
    era5_map = (temp-temp_clim_detrend).sel(time= PNW_heatw_anomaly_day).to_dataset(name="Tx5d_anom")
    era5_map["Z500"] = z500.sel(time= PNW_heatw_anomaly_day).squeeze("plev")/9.81
    return era5_map
        
