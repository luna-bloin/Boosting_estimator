import xarray as xr
from tqdm import tqdm
import glob
import sys
sys.path.append("../utils")
import preproc as pc

## paths
# source data folder
LE100_path = '/net/meso/climphys/fischeer/CESM-ETH/CESM2-LE/'
LE30_path = "/net/meso/climphys/cesm212/"
# processed data
output_path = '/net/xenon/climphys/lbloin/boost_proba/'

def preproc(ds):
    return pc.preprocess(ds).sel(time=slice("2003","2052"))

# === Preprocessing all 100-member LE ===
print("preprocessing 100-member LE")
#selecting all existing PI runs
pnw_le100 = []
for period in tqdm(["1850-2014","2015-2100"]):
    h1_filenames = sorted(glob.glob(f"{LE100_path}tasmax_*_r*i1p1.{period}.nc"))
    #open and preprocess
    pnw_le100.append(xr.open_mfdataset(h1_filenames,preprocess=preproc,combine="nested",concat_dim="member",parallel=True))
pnw_le100["member"] = [int(h1_filenames[i][59:-17]) for i in range(100)]
print("opened")
pnw_le100 = xr.concat(pnw_le100,dim="time")
print("concatenated")
pnw_le100.to_netcdf(f"{output_path}PNW_LE_100_member.nc")
print("saved PNW files")

# === Preprocessing all 30-member LE ===
print("preprocessing 30-member LE")
#selecting all existing PI runs
pnw_le30 = []
for mem in tqdm(range(1,31)): #TODO: what about extra members 31-35?
    file = sorted(glob.glob(LE30_path+f"b.e212.B*cmip6.f09_g17.001.2005.ens{mem:03d}/archive/atm/hist/b.e212.B*cmip6.f09_g17.001.2005.ens{mem:03d}.cam.h1.*-01-01-00000.nc"))
    with xr.open_mfdataset(file,preprocess = preproc) as ds:
        pnw_le30.append(ds.sel(time=slice("2005","2035"))) # there is a file containing only the first day of january 2036
pnw_le30 = xr.concat(pnw_le30,dim = "member")
pnw_le30["member"] = range(1,31)
pnw_le30 = pnw_le30.set_coords('member')
print("opened")
pnw_le30.to_netcdf(f"{output_path}PNW_LE_30_member.nc")
print("saved PNW files")