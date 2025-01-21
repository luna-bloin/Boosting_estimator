import xarray as xr
import numpy as np
import glob
import sys
sys.path.append("../code/utils")
import return_calc as rc
from datetime import datetime, timedelta
import datetime as dt
from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from numpy.random import default_rng
rng = default_rng()

import string
import cartopy.crs as ccrs


pnw_le30 = xr.open_dataset(f"{path}optim_boost/TREFHTMX_PNW_2005-2035.nc").rolling(time=5,center=True).mean()
pnw_le30 = pnw_le30.groupby("time.season")["JJA"].TREFHTMX-273.15

pnw_le100 = xr.open_dataset(f"{path}boost_proba/PNW_LE_100_member.nc")
pnw_le100["member"] = [int(sorted(glob.glob(f"/net/meso/climphys/fischeer/CESM-ETH/CESM2-LE/tasmax_*_r*i1p1.1850-2014.nc"))[i][59:-17]) for i in range(100)]
pnw_le100 = pnw_le100.Tx5d
pnw_le_100_jja = pnw_le100.groupby("time.season")["JJA"]

# open files
era5_2021 = xr.open_dataset(f'{path}boost_proba/temperature_ERA5_2021.nc').resample(time="1D").max().rolling(time=5,center=True).mean().T2M
era5_2021_clim_detrend = xr.open_dataset(f'{path}boost_proba/temperature_detrended_clim_ERA5_2021.nc').resample(time="1D").max().rolling(time=5,center=True).mean().T2M_clim
era5_2021_z = xr.open_dataset(f'{path}boost_proba/tmp_2021_06_z500_ERA5.nc').resample(time="1D").max().rolling(time=5,center=True).mean().Z


# calculate maximum anomaly and time of occurrence
era5_anom = (era5_2021-era5_2021_clim_detrend).sel(lat=slice(45,52),lon=slice(237,241)).mean(("lat","lon"))
era5_anomaly = era5_anom.max()
era5_anomaly_day = era5_anom.idxmax()
print(f"maximum temp anmaly = {era5_anomaly.values:.3}, at {era5_anomaly_day.values}")


era5 = (era5_2021-era5_2021_clim_detrend).sel(time= era5_anomaly_day).to_dataset(name="Tx5d_anom")
era5["Z500"] = era5_2021_z.sel(time= era5_anomaly_day).squeeze("plev")/9.81



# find anomalies, corrected for non-stationarity of the climate using the 100-member LE
anomalies_30 = []
anomalies_100 = []
for year in range(2005,2036):
    # climatology is 100-member mean of 3 years surrounding each year in 30-member LE
    clim = get_clim(year,pnw_le100)
    # find and append anomaly time series for all members
    anomalies_30.append(pnw_le30.sel(time=str(year)).groupby("time.dayofyear") - clim)
    anomalies_100.append(pnw_le_100_jja.sel(time=str(year)).groupby("time.dayofyear") - clim)
anomalies_30 = xr.concat(anomalies_30,dim="time")
anomalies_100 = xr.concat(anomalies_100,dim="time")

# find top annual maximum values
ann_max = anomalies_30.groupby("time.year").max().stack(dim=("member","year"))
ann_max_100 = anomalies_100.groupby("time.year").max().rename({"year": "dim"}).stack(year=("dim","member"))
top = ann_max.sortby(ann_max,ascending=False)[0:13]
events = {}
Tref = 100
peaks = {}
ranks = {}
for i,t in enumerate(top):
    peak = anomalies_30.where(anomalies_30 == t,drop=True).time.dt.strftime("%Y-%m-%d").item()
    mem = int(t.member.values)
    year = int(t.year.values)
    events[f"{mem:02d}:{year}"] = peak
    print(f"{t.dim.values}, {t.values:.3}, {peak}")
    peaks[f"{mem:02d}:{year}"] = peak
    Tref = np.min([Tref,(t.values)])
    ranks[f"{mem:02d}:{year}"] = i + 1

# deleting the events that aren't boosted
for event in ["18:2016","05:2013","17:2009", "25:2012","20:2023","25:2009"]:#,"10:2007"
    del events[event]

#Find and print Tref and Ptref
P_Tref = (ann_max.where(ann_max >=Tref).count()/ann_max.count()).values
P_Tref_100 = (ann_max_100.where(ann_max_100>=Tref).count()/ann_max_100.count()).values
print(f"Tref = {Tref}, P_Tref = {P_Tref}, RE = {((P_Tref-P_Tref_100)/P_Tref_100):.4}")


# preproc boosted runs
full_boost = []
boost_around_peak = []
for event in events:
    # open, get Tx5d
    ds = xr.open_dataset(f"{path}optim_boost/TREFHTMX_PNW_boosted_{event}.nc").TREFHTMX-273.15
    ds = ds.rolling(time=5,center=True).mean()
    # get anomaly values
    year = ds.time.dt.year[0]
    clim = get_clim(year,pnw_le100)
    ds_anom = ds.groupby("time.dayofyear")-clim
    
    #change start date into lead time
    peak = datetime.strptime(events[event], "%Y-%m-%d")
    ds_anom["start_date"] = [(datetime.strptime(date, "%Y-%m-%d") - peak).days for date in ds_anom["start_date"].values]
    full_boost.append(ds_anom)
    # get max in 10-day window around parent peak
    boost_around_peak.append(ds_anom.sel(time=slice((peak - timedelta(days=5)).strftime("%Y-%m-%d"),(peak + timedelta(days=5)).strftime("%Y-%m-%d"))).max("time"))
full_boost = xr.concat(full_boost,dim="case").sel(member=slice(1,100))
full_boost["case"] = list(events.keys())
boost_simple_max = full_boost.max("time")
boost_around_peak = xr.concat(boost_around_peak,dim="case").sel(member=slice(1,100))
boost_around_peak["case"] = list(events.keys())

boost_only_21 = []
for case in full_boost.case:
    b_here_only_21 = []
    b_here = full_boost.sel(case=case).dropna(dim="start_date",how="all")
    for start_date in b_here.start_date:
        b_here_only_21.append(b_here.sel(start_date=start_date).dropna(dim="time",how="all").isel(time=slice(3,21)))
    b_here_only_21 = xr.concat(b_here_only_21,dim= "start_date")
    b_here_only_21["start_date"] = b_here.start_date
    boost_only_21.append(b_here_only_21)
boost_only_21 = xr.concat(boost_only_21,dim="case").max("time")
boost_only_21["case"] = full_boost.case