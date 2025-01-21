# ========================================================================================================================
# === Fig 6: Link between boosted simulations and perturbation time in (blue) test slice 1 and (orange) test slice 2.  ===
# ========================================================================================================================

import matplotlib.pyplot as plt
import sys
sys.path.append("../utils")
import return_calc as rc
import plot_config as pco
import PI_runs_parent_boosted as pi
import utils as ut
import string
from tqdm import tqdm
import xarray as xr

# === Necessary functions ===
def get_spread_time_data(boost):
    """
    Creates a new data array from boost, with a new time axis, which is time elapsed since perturbation
    parameters:
    boost (xarray datarray): array of boosted Tx5d simulations (so with date time axis)
    Returns:
    spread_time_data (xarray datarray): array of boosted Tx5d simulations (so with time since perturbation axis)
    """
    spread_time_data = []
    for i, case in tqdm(enumerate(boost.case.values)):
        boost_case = boost.sel(case=case)
        spread_time_data_case = []
        for j,st in enumerate(boost_case.start_date.values):
            this_boost = boost_case.sel(start_date=st).dropna(dim="member",how="all").dropna(dim="time",how="all")
            spread_time_data_case.append(
                xr.DataArray(
                    data=(this_boost-this_boost.mean("member")).values,
                    dims=["member","spread_time"],
                    coords=dict(
                        member=this_boost.member,
                        spread_time=[t for t in range(len(this_boost.time.values))],
                    ),
                )
            )
        spread_time_data.append(xr.concat(spread_time_data_case,dim="start_date"))
    spread_time_data_all = xr.concat(spread_time_data,dim="case")
    spread_time_data_all["case"] = boost.case.values
    spread_time_data = ut.multi_to_single_index(spread_time_data_all,dims=("start_date","member"),new_name="ID")
    return spread_time_data

f,axs = plt.subplots(1,2,figsize=(7,3.25))

# === Panel a: spread and saturation plot ===
ax = axs[0]
# find climatological standard deviation
le_std_dev = pi.PI_simulation(test_slice="T0").full_simulation.groupby("time.season")["JJA"].groupby("time.dayofyear").std(("time")).mean()

test_slices = ["T1","T2"]
for x,test_slice in enumerate(test_slices): #plot for test slice 1 and 2
    boosted_simulation = pi.boosted_PI_simulation(test_slice=test_slice)
    # create new data array, with time since initialization (spread time) instead of time
    spread_time_data = get_spread_time_data(boosted_simulation.full_simulation)
    # plot std dev divided by climatological std dev
    to_plot = spread_time_data.std(dim="ID")/le_std_dev
    #median
    to_plot.median("case").plot(
        ax=ax,
        label=f"Test slice {x+1}",
        color=pco.colors[x+1]
    )
    #upper and lower bound
    ax.fill_between(
        to_plot.spread_time, 
        to_plot.min("case"), 
        to_plot.max("case"),
        alpha=0.3,
        color=pco.colors[x+1]
    )
# figure configs
ax.axhline(
    1,
    linestyle="dotted", 
    color=pco.colors[0],
    alpha=0.7
)
ax.set_xlabel("Days since perturbation")
ax.set_ylabel(r"Standard deviation relative to climatology")
ax.set_title(r"Spread of boosted simulations with time")
ax.set_xlim(0,None)
ax.set_ylim(0,None)
for i,ax in enumerate(axs):
    pco.set_grid(ax)
    ax.legend(loc="lower right")
    ax.text(0.025,0.94,r"$\textbf{"+string.ascii_lowercase[i]+r"}$",transform=ax.transAxes)
plt.tight_layout()


# === Panel b: Ptref AC plot ===

ax = axs[1]
bt = 1000 #bootstrap
lead_times = range(-18,-6)
for i,test_slice in enumerate(test_slices): #plot for test slice 1 and 2
    boosted_simulation = pi.boosted_PI_simulation(test_slice=test_slice)
    Tref = boosted_simulation.parent_sim.pref_tref()[0]
    # find P(Tref|AC) for all lead times and plot this (with boostrapped error bars)
    p_tref_acs = []
    for lead_time in tqdm(lead_times):
        p_tref_acs.append(
            rc.find_trefAC(
                boosted_simulation.TXx5d,
                lead_time,
                Tref,
                bootstrap = bt
            )
        )
    p_tref_acs = xr.DataArray(p_tref_acs,coords = {"lead_time":lead_times,"bootstrap":range(bt)})
    #plot errorbars
    lower = p_tref_acs.quantile(0.025,dim="bootstrap")
    upper = p_tref_acs.quantile(0.975,dim="bootstrap")
    ax.errorbar(
        p_tref_acs.lead_time,
        p_tref_acs.median("bootstrap"),
        yerr=(p_tref_acs.median("bootstrap")-lower,upper-p_tref_acs.median("bootstrap")),
        capsize=3,
        fmt="o",
        label=f"Test slice {i+1}",
        color=pco.colors[i+1]
    )
#fig configs
pco.set_grid(ax)
ax.set_xlim(-19,-6)
ax.set_ylim(0,1)
ax.set_xlabel("Lead time [days]")
ax.set_ylabel(r"$\hat{p}_{T \geq T_{\mathrm{ref}} \ | \ \mathrm{AC}_t^\epsilon}$")
ax.set_title(r"Effect of lead time on $\hat{p}_{T \geq T_{\mathrm{ref}} \ | \ \mathrm{AC}_t^\epsilon}$")
ax.set_xticks([-17,-14,-11,-8])

plt.savefig(f"{pco.out_path}Fig6_lead_time_effect.png",bbox_inches="tight")