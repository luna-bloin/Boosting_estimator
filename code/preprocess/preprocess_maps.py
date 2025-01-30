import xarray as xr
import pandas as pd
from tqdm import tqdm
import glob

# to run this script, you need to have run preprocess_CESM2_LE.py and ../plots/infor_boosted_LE.py

## paths
# source data folder
LE100_path = '/net/meso/climphys/fischeer/CESM-ETH/CESM2-LE/'
LE30_path = "/net/meso/climphys/cesm212/"
# processed data
output_path = '/net/xenon/climphys/lbloin/boost_proba/'

# open LE100 full Tx5d clim
le100_full = xr.open_dataset(f"{output_path}PNW_LE_100_member_full.nc").Tx5d

for LE_mem in [30,100]:
    maps_of_tops = []
    # read in the top events in this LE and their peak
    peaks = pd.read_csv(f"../../inputs/top_le{LE_mem}.csv",delimiter=";") 
    # loop over each event to get the map of Tx5d and Z500
    for peak in tqdm(peaks["peak"]):
        case = peaks[peaks["peak"] == peak]["case"].item()
        mem = int(case[0:2])
        # file handling for LE30
        if LE_mem == 30:
            if int(peak[0:4]) < 2015:
                period = "HIST"
            else:
                period = "SSP370"
            path = glob.glob(f"{LE30_path}b.e212.B{period}cmip6.f09_g17.001.2005.ens{mem:03d}/archive/atm/hist/b.e212.B*cmip6.f09_g17.001.2005.ens{mem:03d}.cam.h1.{peak[0:4]}-01-01-00000.nc")
            ds = xr.open_dataset(path[0])[["TREFHTMX","Z500"]].rolling(time=5,center=True).mean().sel(time=peak).max("time")
            ds["Tx5d"] = ds.TREFHTMX-273.15
        # file handling for LE100
        elif LE_mem == 100:
            if int(peak[0:4]) < 2015:
                period = "hist"
            else:
                period = "ssp370"
            for var in ["tasmax", "z500"]: # open temperature and z500 (they are in different files)
                file = glob.glob(f"{LE100_path}{var}_{period}_r{mem}i1p1.*.nc")[0]
                ds_tmp = xr.open_dataset(file)
                if var =="tasmax":
                    ds = ds_tmp["TREFHTMX"].to_dataset(name="Tx5d")-273.16
                elif var == "z500":
                    ds["Z500"] = ds_tmp["Z500"]
            ds = ds.rolling(time=5,center=True).mean().sel(time=peak).max("time")
        #temperature climatology of peak day
        clim_here = le100_full.sel(time=peak).max("time")
        ds["lat"] = clim_here["lat"]
        #get temperature anomaly
        ds["Tx5d_anom"] = ds["Tx5d"]-clim_here
        maps_of_tops.append(ds)
    maps_of_tops = xr.concat(maps_of_tops,dim="case")
    maps_of_tops["case"] = list(peaks.case)
    maps_of_tops.to_netcdf(f"{output_path}maps_of_tops_LE{LE_mem}.nc")
