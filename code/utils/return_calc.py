from scipy.stats import genextreme as gev
import numpy as np
import xarray as xr
from tqdm import tqdm
import matplotlib.pyplot as plt
from numpy.random import default_rng
rng = default_rng()

# ============================================================================
# === Functions for naive return time estimation (non-boosted simulations) ===
# ============================================================================

def fit_gev(data):
    # Fit the GEV distribution to your data and return the shape, location, and scale
    shape, loc, scale = gev.fit(data)
    return shape, loc, scale

def gev_return_level(shape, loc, scale, return_time):
    # Calculate return level for a given return time T
    return gev.isf(1/return_time, shape, loc=loc, scale=scale)    
    
def gev_return_lev_for_time_series(data,max_return_time):
    # Calculate return level for a  return time series
    shape, loc, scale = fit_gev(data)
    # Define return times, including beyond the dataset
    start = len(data)/(len(data)-1) # start at second point in data set
    return_times = np.logspace(np.log10(start), np.log10(max_return_time+1),num=1000)  # From 1-year to `max_return_time`-year return times
    return_levels = return_level(shape, loc, scale, return_times)
    return return_levels

def return_time_bootstrap(data,bootstrap = 1,max_return_time = 100000):
    # find the return times of a given data set of return levels data, with option for bootstrap
    new_maxes = []
    bootstrap_data = rng.choice(data, size=(bootstrap, len(data)), replace=True)
    bootstrap_data = xr.DataArray(data=bootstrap_data,coords={"bootstrap":range(bootstrap),"year":data.year})
    for i in tqdm(range(bootstrap)):
        selected_years = bootstrap_data.sel(bootstrap=i)
        new_maxes.append(return_lev_for_time_series(selected_years,max_return_time))
    
    start = len(data)/(len(data)-1) # start at second point in data set
    sorted_new_maxes = xr.DataArray(
        new_maxes, 
        coords={
            'bootstrap': list(range(bootstrap)), 
            'return_time': np.logspace(np.log10(start), np.log10(max_return_time+1),num=1000)  # From 1-year to `max_return_time`-year return times
        }, 
        dims=['bootstrap','return_time']
    )
    return sorted_new_maxes

def naive_estimator(TXx5d):
    return_levels_sorted = TXx5d.sortby(TXx5d,ascending=False)
    return_times = len(TXx5d)/np.arange(1, len(TXx5d) + 1)
    return return_times, return_levels_sorted
    
# ===============================================================================
# === Functions for boosting return time estimation (non-boosted simulations) ===
# ===============================================================================

def boosting_estimator(TXx5d, lead_time, Tref, P_Tref, bootstrap=1000):
    """
    calculates the return times from a data array TXx5d with the boosting estimator

    Parameters:
    TXx5d (xarray dataarray): Data array of TXx5d of boosted simulations
    lead_time (int, slice or None): if int (or slice), the lead time (or lead time range) at which to calculate. If None, it pools all lead times in TXx5d together
    Tref (float): reference temperature
    P_Tref (float): probability of reference temperature in non-boosted simulation
    bootstrap (int): number of times to repeat calculation
    Returns:
        return_time_da (xarray dataarray): data array of bootstrapped return times and return values
    """
    #pool all TXx5d together
    if lead_time == None:
        TXx5d = TXx5d.stack(dim=("case","member","start_date")).dropna(dim="dim")
    else:
        TXx5d = TXx5d.sel(start_date=lead_time).stack(dim=("case","member")).dropna(dim="dim")
    # find all TXx5d above Tref
    above_Tref_all = TXx5d.where(TXx5d>= Tref,drop=True)
    T_ext_all = boost_above_Tref_all.sortby(boost_above_Tref_all,ascending=False)
    
    bootstrap_TXx5d = xr.DataArray(data = rng.choice(TXx5d, size=(bootstrap, len(TXx5d)), replace=True))
    return_times_T_ext = []
    # find P(Tref|AC) for each bootstrap sample
    P_Tref_AC = bootstrap_TXx5d.where(bootstrap_TXx5d>= Tref).count(dim="dim_1").values/bootstrap_TXx5d.count(dim="dim_1").values
    # loop over all T_ext (>= Tref) values to calculate probability
    for text in tqdm(T_ext_all):
        # find P(T_ext|AC)
        P_Text_AC = bootstrap_TXx5d.where(bootstrap_TXx5d>= text).count(dim="dim_1").values / bootstrap_here.count(dim="dim_1").values
        #find P(T_ext)
        P_Text = (P_Tref * (P_Text_AC/P_Tref_AC))
        return_times_T_ext.append(1/P_Text)
    return_times_T_ext = xr.DataArray(data = return_times_T_ext, dims=["T_ext","bootstrap"],
    coords=dict(
        T_ext=T_ext_all.values,
        bootstrap=range(bootstrap),
    ),)
    return return_times_T_ext

def find_trefAC(TXx5d,lead_time,Tref,bootstrap = 1):
    """
    calculates P(Tref | AC) in a boosted data array of TXx5d

    Parameters:
    TXx5d (xarray dataarray): Data array of TXx5d of boosted simulations
    lead_time (int or None): if int, the lead time at which to calculate. If None, it pools all lead times in TXx5d together
    Tref (float): reference temperature
    bootstrap (int): number of times to repeat calculation
    Returns:
        P_Tref_AC (list): list of bootstrapped values of P(Tref | AC)
    """
    if ld == None:
        boost_here = boost_here.stack(dim=("case","member","start_date")).dropna(dim="dim")
    else:
        boost_here = boost_here.sel(start_date=ld).stack(dim=("case","member")).dropna(dim="dim")    
    if bootstrap == 1:
        boost_above_Tref = boost_here.where(boost_here>= Tref,drop=True)
        P_Tref_AC = boost_above_Tref.count().values/boost_here.count().values 
    else:
        bootstrap_here = xr.DataArray(data = rng.choice(boost_here, size=(bootstrap, len(boost_here)), replace=True))
        P_Tref_AC = []
        for i in range(bootstrap):
            boost_above_Tref = bootstrap_here[i].where(bootstrap_here[i]>= Tref,drop=True)
            P_Tref_AC.append(boost_above_Tref.count().values/boost_here.count().values )
    return P_Tref_AC

