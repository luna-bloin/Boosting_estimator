from scipy.stats import genextreme as gev
import numpy as np
import xarray as xr
from tqdm import tqdm
import matplotlib.pyplot as plt
from numpy.random import default_rng
rng = default_rng()

def fit_gev(data):
    # Fit the GEV distribution to your data and return the shape, location, and scale
    shape, loc, scale = gev.fit(data)
    return shape, loc, scale

def return_level(shape, loc, scale, return_time):
    # Calculate return level for a given return time T
    return gev.isf(1/return_time, shape, loc=loc, scale=scale)    
    
def return_lev_for_time_series(data,max_return_time):
    shape, loc, scale = fit_gev(data)
    # Define return times, including beyond the dataset
    start = len(data)/(len(data)-1) # start at second point in data set
    return_times = np.logspace(np.log10(start), np.log10(max_return_time+1),num=1000)  # From 1-year to `max_return_time`-year return times
    return_levels = return_level(shape, loc, scale, return_times)
    return return_levels

def return_time_bootstrap(data,bootstrap = 1,max_return_time = 10000):
    # find the return times of a given data set data, with option for bootstrap
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

def plot_return(data,sorted_new_maxes,label,ax,color="skyblue",msize=5):
    # plot the return period
    ax.plot(len(data)/np.arange(1, len(data) + 1), data.sortby(data,ascending=False), ".",label=label,markersize=msize,color=color)
    sorted_new_maxes.median("bootstrap").plot(color=color,ax=ax)
    # to plot 
    x = sorted_new_maxes.return_time
    quant_low = sorted_new_maxes.quantile(0.025,"bootstrap")
    quant_high = sorted_new_maxes.quantile(0.975,"bootstrap")
    ax.fill_between(x,quant_low,quant_high,alpha=0.3,color=color)
    return None

def calculate_ret(ld,boost_here,Tref, P_Tref):
    # select from boosted dataset
    if ld == None:
        boost_here = boost_here.stack(dim=("case","member","start_date")).dropna(dim="dim")
    else:
        boost_here = boost_here.sel(start_date=ld).stack(dim=("case","member")).dropna(dim="dim")
    # find parameters for equation
    boost_above_Tref = boost_here.where(boost_here>= Tref,drop=True)
    P_Tref_AC = boost_above_Tref.count().values/boost_here.count().values 
    Text = boost_above_Tref.sortby(boost_above_Tref,ascending=False)
    P_Text_AC = np.arange(1,len(boost_above_Tref)+1) / len(boost_here)
    # equation
    P_Text = P_Tref * (P_Text_AC/P_Tref_AC)
    return_times = 1/P_Text
    return Text,return_times

def calculate_ret_bootstrap(ld, boost_here, Tref,P_Tref,bootstrap=1000):
    if ld == None:
        boost_here = boost_here.stack(dim=("case","member","start_date")).dropna(dim="dim")
    else:
        boost_here = boost_here.sel(start_date=ld).stack(dim=("case","member")).dropna(dim="dim")
    boost_above_Tref_all = boost_here.where(boost_here>= Tref,drop=True)
    T_ext_all = boost_above_Tref_all.sortby(boost_above_Tref_all,ascending=False)
    
    bootstrap_here = xr.DataArray(data = rng.choice(boost_here, size=(bootstrap, len(boost_here)), replace=True))
    return_times_text = []
    P_Tref_AC = bootstrap_here.where(bootstrap_here>= Tref).count(dim="dim_1").values/bootstrap_here.count(dim="dim_1").values
    for text in tqdm(T_ext_all):
        P_Text_AC = bootstrap_here.where(bootstrap_here>= text).count(dim="dim_1").values / bootstrap_here.count(dim="dim_1").values
        P_Text = (P_Tref * (P_Text_AC/P_Tref_AC))
        return_times_text.append(1/P_Text)
    return_times_text = xr.DataArray(data = return_times_text, dims=["T_ext","bootstrap"],
    coords=dict(
        
        T_ext=T_ext_all.values,
        bootstrap=range(bootstrap),
    ),)
    return return_times_text

def find_trefAC(ld,boost_here,Tref,bootstrap = 1):
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