# === Script to open the PI_control runs =================
# === T0: the total reference climate model simulation ===
import pandas as pd
import plot_config as pco
import return_calc as rc
import xarray as xr
from tqdm import tqdm
from glob import glob


time_periods = {
    "T1": [1801,1850],
    "T2": [1851,1900],
}

plot_labels = {"T0":r"Control period: $N =4000$",
          "T1":'Test slice: $N = 50$',
          "T2":'Test slice: $N = 50$',
          "T3":'Test slice: $N = 100$',
         }

def read_parents(test_slice):
    """
    Reads selected_events.csv to find the parents in each test slice (T1, T2 or T3)

    Parameters:
    test_slice (str): Either T1, T2, T3
    Returns:
        list: A list of years (integers) from the specified test slice.
    """
    # Read the CSV file
    df = pd.read_csv('../../inputs/selected_events.csv', sep=';', dtype=str)
    if test_slice == 'T3': # reads from both T1 and and T2 and combines it into one list
        years = pd.concat([df[col].astype(int) for col in ["T1","T2"]]).tolist()
    else:
        years = df[test_slice].astype(int).tolist()
    return years


class PI_simulation():
    def __init__(self, test_slice):
        """
        Initializes the PI_simulation class.
    
        Parameters:
        test_slice (str): test slice. Either T1, T2, T3 or T0 (for entire reference simulation)
    
        Returns:
        None
        """
        #assign attributes
        self.test_slice = test_slice
        if self.test_slice == "T0":
            self.full_simulation = xr.open_dataset(f"{pco.path}PNW_PI_control_simulation.nc").Tx5d_anom.sel(time=slice(None,"4000"))
        else:
            self.full_simulation = xr.open_dataset(f"{pco.path}PNW_PI_test_slice_{self.test_slice}.nc").Tx5d_anom
        self.TXx5d = self.full_simulation.groupby("time.season")["JJA"].groupby("time.year").max("time")
        # find return time for all TXx5d
        self.bootstrapped_GEV = rc.return_time_bootstrap(self.TXx5d,bootstrap = 1000)
        # Only the TXx5d of the years that are selected to be boosted
        if test_slice != "T0":
            self.parents = self.TXx5d.sel(year=read_parents(self.test_slice))

    def pref_tref(self):
        # find Tref and its probability
        if self.test_slice != "T3":
            top = 4
            Tref = self.TXx5d.sortby(self.TXx5d ,ascending=False)[top] # 5th highest value
            P_Tref = 0.1 #by construction 
        else: # if test slice == T3, Tref is the minimum of Tref from T1, T2, and Ptref must be calculated
            year = min([PI_simulation(test_slice = t).pref_tref()[0] for t in ["T1","T2"]]).year
            Tref = self.TXx5d.where(self.TXx5d.year == year, drop=True).squeeze("year")
            P_Tref = (self.TXx5d.where(self.TXx5d>=Tref,drop=True).count()/self.TXx5d.count()).values
        return Tref, P_Tref

    def calculate_and_plot_return_period(self,ax):
        # plot non-boosted
        for j,pi_sim in enumerate([PI_simulation(test_slice="T0"),self]):
            pco.plot_return(
                pi_sim.TXx5d,
                pi_sim.bootstrapped_GEV,
                plot_labels[pi_sim.test_slice],
                ax,
                color=pco.colors[j],
                msize=3
            )
        #highlight parent events
        if self.test_slice == "T3": #go parent by parent to find the correct return time index (since the selected 10 events in T3 are not top 10)
            sorted_TXx5d = self.TXx5d.sortby(self.TXx5d,ascending=False)
            for i,time_slice in enumerate(["T1","T2"]):
                for j,parent in enumerate(self.parents):
                    if i == 0 and j == 0:
                        lab = "Parent Events"
                    else:
                        lab = "_nolab"
                    index = len(self.TXx5d.where(self.TXx5d >= parent.values,drop=True))
                    pco.plot_parents(
                        ax,
                        rc.naive_estimator(self.TXx5d)[0][index-1],
                        parent.values,
                    )
        else:
            return_times,return_levels = rc.naive_estimator(self.TXx5d)
            pco.plot_parents(
                        ax,
                        return_times[0:5],
                        return_levels[0:5]
                    )

class boosted_PI_simulation():
    def __init__(self, test_slice=None):
        """
        Initializes the boosted_PI_simulation class.
    
        Parameters:
        test_slice (str): test slice. Either T1, T2, T3 (for all boosted simulations)
    
        Returns:
        None
        """
        self.test_slice = test_slice
        #parent simulations
        self.parent_sim = PI_simulation(test_slice)
        # opening boosted runs
        files = glob(f'{pco.path}PNW_PI_boosted_*.nc')
        # open all boosted realizations
        full_boost = xr.open_mfdataset(files,concat_dim ="case",combine="nested").Tx5d_anom
        full_boost["case"] = [int(f[54:58]) for f in files]
        self.full_simulation = full_boost.sel(case = read_parents(self.test_slice)).dropna(dim="time",how="all").load()
        self.TXx5d = self.full_simulation.max("time")
        
    def plot_boosted_PI(self,lead_time,ax):
        Tref,P_Tref = self.parent_sim.pref_tref()
        # calculate boosting estimator (with bootstrap)
        ret = rc.boosting_estimator(
            self.TXx5d,
            lead_time, 
            Tref,
            P_Tref,
        )
        #plot median
        ret.median("bootstrap").plot.line(
            ".",
            color=pco.colors[2],
            label=fr"Boosted simulations, $N_b = {len(ret)}$",
            markersize=3,
            y="T_ext",
            ax=ax
        )
        #plot upper and lower bound
        ax.fill_betweenx(
            ret.T_ext, 
            ret.quantile(0.025,"bootstrap"),
            ret.quantile(0.975,"bootstrap").fillna(10000000), # so that it counts the infs (up to 10^4 which is the figure limit)
            color=pco.colors[2],
            alpha=0.5
        )
        ax.set_xscale('log')
        ax.set_xlim(.9,3e4)
        ax.set_ylim(6,14.2)
        ax.axhline(Tref,linestyle="--",color=pco.colors[2],label=r'$T_{\mathrm{ref}}$')