# === Script to open the LE runs =================
import plot_config as pco
import return_calc as rc
import utils as ut
import xarray as xr
from tqdm import tqdm
from glob import glob

path = "/net/xenon/climphys/lbloin/boost_proba/"

class LE_simulation():
    def __init__(self, tot_members, anomalies=False):
        """
        Initializes the PI_simulation class.
    
        Parameters:
        tot_members (int): number of members in LE to open. Either 30 or 100
        """
        #assign attributes
        self.tot_members = tot_members
        
        self.full_simulation = xr.open_dataset(f"{pco.path}PNW_LE_{self.tot_members}_member.nc").Tx5d
        if anomalies == True:            
            # # find Tx5d anomalies corrected for non-stationarity
            anomalies = []
            for year in range(2005,2036):
                # climatology is 100-member mean of 3 years surrounding each year in 30-member LE
                clim = ut.get_clim(year,LE_simulation(tot_members=100,anomalies=False).full_simulation)
                # find and append anomaly time series for all members
                anomalies.append(self.full_simulation.sel(time=str(year)).groupby("time.dayofyear") - clim)
            anomalies = xr.concat(anomalies,dim="time")
            self.full_simulation = anomalies
        self.TXx5d = self.full_simulation.groupby("time.season")["JJA"].groupby("time.year").max("time")
        TXx5d_1d = self.TXx5d.rename({'year':'dim'}).stack(year=("member","dim")) # stacking years and members to get a 1D dataarray
        self.TXx5d_1D_sorted = TXx5d_1d.sortby(TXx5d_1d,ascending=False)

    def find_top_TXx5d(self,top=13):
        # finds the member and year of events with highest intensity. top denotes how many (ordered) events to output
        top_events = self.TXx5d_1D_sorted[0:top]
        top_event_info = {}
        for i,t in enumerate(top_events):
            peak = self.full_simulation.where(self.full_simulation == t,drop=True).time.dt.strftime("%Y-%m-%d").item() # date of highest intensity
            mem = int(t.member.values)
            year = int(t.dim.values)
            top_event_info[f"{mem:02d}:{year}"] = peak
        return top_events, top_event_info

    def find_parents(self,remove_parent=[]):
        parent_files = glob(f"{path}/PNW_LE_boosted*.nc")
        parent_list = [(int(parent_files[i][54:56]),int(parent_files[i][57:61])) for i in range(len(parent_files))]
        if len(remove_parent) > 0:
            for parent in remove_parent:
                parent_list.remove((int(parent[0:2]),int(parent[3:])))
        return self.TXx5d_1D_sorted.sel(year=parent_list)
    
    def pref_tref(self):
        top_events = self.find_top_TXx5d()[0]
        Tref = top_events[-1].values
        P_Tref = rc.find_proba_naive(self.TXx5d_1D_sorted,Tref,bootstrap=1000) #(self.TXx5d.where(self.TXx5d >=Tref).count()/self.TXx5d.count()).values
        return Tref, P_Tref

    def bootstrapped_GEV(self,remove_a_parent=[]):
        # find return time for all TXx5d. remove_a_parent is either false, for the full return time estimation, or equal a parent to remove (format str, "MM:YYYY" M=member, Y=year)
        if len(remove_a_parent) == 0:
            return rc.return_time_bootstrap(self.TXx5d_1D_sorted,bootstrap = 1000)
        else:
            new_sorted_TXx5d_1D = self.get_new_sorted_TXx5d_1D(remove_a_parent[0])
            return rc.return_time_bootstrap(new_sorted_TXx5d_1D,bootstrap = 1000)
    
    def get_new_sorted_TXx5d_1D(self,remove_a_parent):
        return self.TXx5d_1D_sorted.where(self.TXx5d_1D_sorted != self.TXx5d_1D_sorted.sel(year=(int(remove_a_parent[0:2]),int(remove_a_parent[3:]))),drop=True)
    
    def calculate_and_plot_return_period(self,ax,remove_a_parent=[],color=pco.colors[1]):
        # plot non-boosted
        if len(remove_a_parent) > 0:
            TXx5d = self.get_new_sorted_TXx5d_1D(remove_a_parent[0])
        else:
            TXx5d = self.TXx5d_1D_sorted
        pco.plot_return(
            TXx5d,
            self.bootstrapped_GEV(remove_a_parent=remove_a_parent),
            f"{self.tot_members}-member Large Ensemble",
            ax,
            color=color,
            msize=3
        )
        #highlight highest parent events
        if self.tot_members == 100:
            # For Le 100, this means highlighting ~E1 with a star
            pco.plot_highest_parent(
                ax,
                len(TXx5d),
                TXx5d[0],
                r"$\tilde{E}1$",
                color=pco.colors[0]
            )
        else:
            # For LE 30, this means highlighting E1 with a star and E2-13 with diamonds
            if remove_a_parent != ["10:2007"]:
                pco.plot_highest_parent(
                    ax,
                    len(TXx5d),
                    TXx5d[0],
                    r"$E1$",
                    color=pco.colors[1]
                )
            #highlight rest of parent events
            return_times,return_levels = rc.naive_estimator(TXx5d)
            # remove chosen parent + 10:2007 from highlight
            if len(remove_a_parent) >0:
                parents = self.find_parents(remove_parent=list(set(["10:2007",remove_a_parent[0]])))
            else:
                parents = self.find_parents(remove_parent=["10:2007"])
            # highlight parent by parent
            for i,parent in enumerate(parents):
                if i == 0:
                    label = True
                else:
                    label = False
                levels_above = return_levels.where(return_levels >= parent,drop=True)
                level = levels_above[-1].values
                time = return_times[len(levels_above)-1]
                pco.plot_parents(
                    ax,
                    time,
                    level,
                    label = label
                )

class boosted_LE_simulation():
    def __init__(self,drop_case=[]):
        """
        Initializes the boosted_LE_simulation class.
        Parameters:
        drop_case (list): list of cases whose simulations you wouldn't take into account

        """
        
        # opening boosted runs
        files = glob(f'{pco.path}PNW_LE_boosted_*.nc')
        # open all boosted realizations
        full_boost = xr.open_mfdataset(files,concat_dim ="case",combine="nested").Tx5d_anom.sel(member=slice(1,100)) # 100 members for each lead time for each parent
        full_boost["case"] = [f[54:61] for f in files]
        # remove boosted simulations from one or several parents if drop_case isnt empty:
        if len(drop_case) > 0:
            print(f"droping cases {drop_case}")
            for case in drop_case:
                full_boost = full_boost.drop_sel(case=case)
        self.full_simulation = full_boost.load()
        self.TXx5d = self.full_simulation.max("time")
        #parent simulations
        self.parent_sim = LE_simulation(30,anomalies=True)
        
    def plot_boosted_LE(self,lead_time,ax,era5_anomaly=None):
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
            label=r"Boosted simulations, $N_b$ = " + str(round(self.TXx5d.sel(start_date=lead_time).count().item(),-2)),
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
        ax.set_xlim(.9,1e5)
        ax.set_ylim(8,None)
        ax.axhline(Tref,linestyle="--",color=pco.colors[2],label=r'$T_{\mathrm{ref}}$')
        
