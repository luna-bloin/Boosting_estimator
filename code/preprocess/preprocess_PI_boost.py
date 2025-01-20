import xarray as xr
import glob
import datetime as dt
import sys
sys.path.append("../utils")
import preproc as pc


# === Preprocess boosted runs, resulting from PiControl_select.ipynb, in same way ===
# chosen cases
cases_boosted = { 
    'T1':{
        "1808-06-06": ["1808-05-19", "1808-05-30"],
        "1832-08-24": ["1832-08-06", "1832-08-17"],
        "1828-08-09": ["1828-07-22", "1828-08-02"],
        "1835-06-22": ["1835-06-04", "1835-06-15"],
        "1829-06-17": ["1829-05-30", "1829-06-10"],
    },
    'T2':{
        "1857-07-15": ["1857-06-27", "1857-07-08"],
        "1878-08-30": ["1878-08-12", "1878-08-23"],
        "1895-07-17": ["1895-06-27", "1895-07-10"],
        "1870-06-25": ["1870-06-07", "1870-06-18"],
        "1856-08-07": ["1856-07-20", "1856-07-31"],
    },
}

for tim in cases_boosted:
    print(tim)
    clim_tim = xr.open_dataset(f'{pc.output_path}PI_test_slice_{tim}.nc').Tx5d.rolling(time=30,center=True).mean().groupby("time.dayofyear").mean()
    for case in cases_boosted[tim]:
        print(case)
        # loop over all dates for each case
        date_range = cases_boosted[tim][case]
        start_date = dt.date(int(date_range[0][0:4]),int(date_range[0][5:7]),int(date_range[0][8:10]))
        end_date = dt.date(int(date_range[1][0:4]),int(date_range[1][5:7]),int(date_range[1][8:10]))
        peak_date = dt.date(int(case[0:4]),int(case[5:7]),int(case[8:10]))
        current_date = start_date
        boost = []
        dates = []
        while current_date <= end_date:
            print(int((current_date - peak_date).days))
            # open all 100 boosts for each date
            h1_filenames_boost = sorted(glob.glob(f"/net/meso/climphys/cesm212/boosting_piControl/archive/B1850cmip6.1000001.{current_date}*/atm/hist/B1850cmip6.1000001.{current_date}*.cam.h1.*-00000.nc"))
            members = [int(ls[87:90]) for ls in h1_filenames_boost] # find the exact member number from the file (in case one is missing)
            if h1_filenames_boost != []:
                ds = xr.open_mfdataset(h1_filenames_boost,combine="nested", concat_dim="member", coords='minimal',preprocess=pc.preprocess).isel(time=slice(0,21))
                print("opened")
                ds["member"] = members
                # add anomalies
                ds["Tx5d_anom"] = ds.Tx5d.groupby("time.dayofyear")-clim_tim
                boost.append(ds)
                # next date
                dates.append(current_date)
            current_date += dt.timedelta(days=1)
        boost = xr.concat(boost,dim="start_date")
        boost["start_date"] = [int((dt - peak_date).days) for dt in dates]
        boost.to_netcdf(f'{pc.output_path}PNW_PI_boosted_{case}.nc')
        print("saved")
