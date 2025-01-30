#================================
# === Fig 5: Hummingbird plot ===
#================================

import matplotlib.pyplot as plt
import sys
sys.path.append("../utils")
import return_calc as rc
import plot_config as pco
import PI_runs_parent_boosted as pi
import string


f,axs = plt.subplots(2,5,figsize=(7,4),sharey=True,sharex=True)
ax = axs.flatten()
#open all boosted simulations
boost_sim = pi.boosted_PI_simulation(test_slice="T3")
TXx5d = boost_sim.TXx5d 
parents = boost_sim.parent_sim.parents

#loop over all parent cases
for i,case in enumerate(TXx5d.case): 
    # plot median + 90th percentile range of boosted
    to_plot = TXx5d.sel(start_date=slice(-18,-7),case=case)
    upper =  to_plot.quantile(0.95,dim="member")
    lower =  to_plot.quantile(0.05,dim="member")    
    to_plot.median("member").plot(ax=ax[i],color=pco.colors[2],zorder=3,label="Batch median")
    ax[i].fill_between(to_plot.start_date,upper,lower,alpha=0.5,color=pco.colors[2], label=r"5$^{\mathrm{th}} - $ 95$^{\mathrm{th}}$ percentile values")
    # plot max boosted
    to_plot.max("member").plot.line("+",color=pco.colors[2],ax=ax[i],label = "Batch maximum")
    # plot parent TXx5d
    ax[i].axhline(parents.sel(year=case),linestyle="--",color=pco.colors[1],label="Parent")
    # fig configs
    if i < 5:
        ax[i].set_xlabel("")
    else:
        ax[i].set_xlabel("Lead time [days]")
    if i%5!=0:
        ax[i].set_ylabel("")
    else:
        ax[i].set_ylabel("TXx5d [$^\circ$C]")
    pco.set_grid(ax[i])
    ax[i].text(0.025,0.9,r"$\textbf{"+string.ascii_lowercase[i]+r"}$",transform=ax[i].transAxes)
    if i < 5:
        ax[i].set_title(f"Parent {i%5+1}")
    else:
        ax[i].set_title("")
handles, labels = ax[0].get_legend_handles_labels()
f.legend(handles, labels, loc='lower center',ncol=4)
f.text(0.5, 0.95, "Test slice 1", ha='center', fontsize=14) 
f.text(0.5, 0.51, "Test slice 2", ha='center', fontsize=14)
plt.tight_layout(rect=[0, 0.42, 1, 0.95])
plt.subplots_adjust(bottom=0.21)
pco.convert_ticklabels_to_strings(f)
plt.savefig(f"{pco.out_path}Fig5_hummingbird.png", bbox_inches="tight")