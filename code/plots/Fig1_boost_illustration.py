# =====================================================================================================================
# === Code necessary to reproduce Figure 1: Illustration of the Ensemble Boosting algorithm for heatwaves =============
# === This code reproduces the panels in the figure separately, which are then put together in Microsoft PowerPoint ===
# =====================================================================================================================

import sys
sys.path.append("../utils")
import plot_config as pco
import PI_runs_parent_boosted as pi
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import string

test_slice = pi.PI_simulation(test_slice="T1")
## Panel a: initial climate model simulation
tref,ptref = test_slice.pref_tref()

f,ax=plt.subplots(figsize=(7,1))
##plot non-boost
test_slice.TXx5d.plot(
    color=pco.colors[0], 
    label="Parent Ensemble",
    marker=".",
    linestyle=""
)
pco.set_grid(ax)
plt.ylabel("TXx5d [$^\circ$C]")
plt.xlabel("Year")
ax.text(
    0.015,
    0.85,
    r"$\textbf{"+string.ascii_lowercase[0]+r"}$",
    transform=ax.transAxes
)
plt.savefig(
    f"{pco.out_path}illustration_annual_max.svg",
    transparent=True,
    bbox_inches="tight",
)

# Panel b-e: Boosted runs for different lead times

f,ax=plt.subplots(
    5,
    1,
    figsize=(7,6.5),
    sharey=True,
    height_ratios=[1, 0.5, 1, 1,1]
)
ax[1].axis('off') # white space

#choose parent event
parent_year = 1828
parent = test_slice.full_simulation.sel(time=str(parent_year)).groupby("time.season")["JJA"].copy().groupby("time.dayofyear").mean()
peak = parent.idxmax()
parent["lead_time"] = parent.dayofyear-peak

# open boosted parent
boost_simulation_one_case = pi.boosted_PI_simulation("T3").full_simulation.sel(case = parent_year)

#plot boosted realizations
lead_times = [-16,-12,-7]
for i,dt in enumerate(lead_times):
    boost_here = boost_simulation_one_case.groupby("time.season")["JJA"].groupby("time.dayofyear").mean().sel(start_date=dt)
    boost_here["lead_time"] = boost_here.dayofyear - peak
    boost_here.plot(
        x="lead_time",
        ax=ax[i+2],
        hue="member",
        add_legend=False,
        color=pco.colors[2],
        alpha=0.3
    )
    parent_here = parent.where(parent.lead_time >= dt,drop=True).where(parent.lead_time <= dt+5,drop=True)
    parent_here.plot(
        x="lead_time", 
        color=pco.colors[2],
        linewidth=1.2,
        ax=ax[i+2]
        )
# figure style + plot parent
for i in range(len(ax)):
    if i ==1:
        continue
    if i > 1:
        x = i
    else:
        x = i+1
    ax[i].text(
        0.015,
        0.86,
        r"$\textbf{"+string.ascii_lowercase[x]+r"}$",
        transform=ax[i].transAxes
    )
    #plot parent event everywhere
    parent.plot(
        x="lead_time",
        ax=ax[i],
        color=pco.colors[0],
        linewidth=1,
        zorder=3,
        label="Parent event"
    ) 
    #fig configs
    pco.set_grid(ax[i])
    ax[i].set_xlim(-17,15)
    ax[i].set_title("")
    if i < 4 and i !=0:
        ax[i].set_xlabel("")
    else:
        ax[i].set_xlabel("Lead time [days]")
    ax[i].set_ylabel("Tx5d anomaly [$^\circ$C]")
    ax[i].axvline(0,color=pco.colors[0],linewidth=1)
    if i != 0:
        ax[i].axvline(
            lead_times[i-2],
            color=pco.colors[0],
            linestyle="--",
            linewidth=1
        )
ax[0].legend(loc="lower left")
ax[4].legend(handles = [Patch(facecolor=pco.colors[2], label='Boosted batch')],loc="lower left")
for i in [2,3]:
    ax[i].set_xticklabels([])
plt.savefig(
    f"{pco.out_path}illustration_boost.svg",
    transparent=True,
    bbox_inches="tight",
)
