import xarray as xr
import glob
import csv
import datetime as dt
import sys
sys.path.append("../utils")
import preproc as pc
import utils as ut

cases_boosted = { 
    #"10:2017": [dt.date(2017, 7, 10), dt.date(2017, 7, 28)],
    #"10:2007": [dt.date(2007, 7, 13), dt.date(2007, 7, 28)],
    "29:2033": [dt.date(2033, 6, 23), dt.date(2033, 7, 4)],
    "06:2013": [dt.date(2013, 6, 5), dt.date(2013, 6, 30)],
    "12:2031": [dt.date(2031, 7, 1), dt.date(2031, 8, 24)],
    "23:2033": [dt.date(2033, 5, 25), dt.date(2033, 6, 11)],
    "18:2034": [dt.date(2034, 6, 23), dt.date(2034, 7, 5)],
}
        
pnw_le30 = xr.open_dataset(f'{pc.output_path}PNW_LE_30_member.nc').Tx5d
pnw_le100 = xr.open_dataset(f'{pc.output_path}PNW_LE_100_member.nc').Tx5d
pnw_le100["member"] = [int(sorted(glob.glob(f"/net/meso/climphys/fischeer/CESM-ETH/CESM2-LE/tasmax_*_r*i1p1.1850-2014.nc"))[i][59:-17]) for i in range(100)]

for case in cases_boosted:
    print(case)
    # find parent peak
    parent = pnw_le30.sel(member=int(case[0:2]),time=case[3:7])
    peak = parent.idxmax().item().dayofyr
    # find climatology (corrected for non-stationarity)
    clim = ut.get_clim(int(case[3:7]),pnw_le100)
    # loop over all dates for each case
    mem = case[0:-5]
    start_date = cases_boosted[case][0]
    end_date = cases_boosted[case][1]
    current_date = start_date
    #file naming conventions depending on year
    if int(case[3:]) < 2015:
        fi_len = [77,80] 
    else:
        fi_len = [79,82] 
    # open and preprocess for all lead times
    boost = []
    dates = []
    while current_date <= end_date:
        h1_filenames_boost = sorted(glob.glob(f"{pc.boost_path}B*cmip6.000*{mem}.{current_date}.ens*/atm/hist/B*cmip6.*.*.ens*.cam.h1.*-00000.nc"))
        files = [fi for fi in h1_filenames_boost if "old" not in fi] #don't include old files
        print(current_date,len(files))
        if files != []:
            dates.append(current_date.timetuple().tm_yday)
            ds = xr.open_mfdataset(files, preprocess=pc.preprocess,concat_dim="member", combine="nested",parallel=True)
            print("opened")
            members =  [int(fi[fi_len[0]:fi_len[1]]) for fi in files] # find the exact member number from the file (in case one is missing)
            print(members)
            ds["member"] = members
            ds = ds.set_coords('member')
            # add anomalies
            ds["Tx5d_anom"] = ds.Tx5d.groupby("time.dayofyear")-clim
            boost.append(ds.isel(time=slice(0,21)))
        current_date += dt.timedelta(days=1)
    # gathering all lead times into one file
    boost=xr.concat(boost,dim="start_date")
    boost["start_date"] = [int(dt - peak) for dt in dates]
    boost = boost.set_coords('start_date')
    print("processed")
    boost.to_netcdf(f"{pc.output_path}PNW_LE_boosted_{case}.nc")
    print("saved")
