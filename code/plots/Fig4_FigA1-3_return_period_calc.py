# ============================================================================================================================================================
# === Fig 4 (+ appendix figure 1-3): Return periods of TXx5d in the PNW region estimated with the boosting estimator under stationary climate conditions.  ===
# ============================================================================================================================================================

import matplotlib.pyplot as plt
import sys
sys.path.append("../utils")
import return_calc as rc
import PI_runs_parent_boosted as pi
import plot_config as pco
import string

# === Fig 4 ===
fig,axs = plt.subplots(2,2,figsize = (7,5),sharex=True,sharey=True)
test_slices = ["T1","T2"]

# calculate return time for each slice of lead time
lead_times = [slice(-18,-13),slice(-12,-7)]
for i,test_slice in enumerate(test_slices):
    # get boosted simulation
    boosted_simulation = pi.boosted_PI_simulation(test_slice=test_slice)
    for x,lead_time_range in enumerate(lead_times):
        ax = axs[i][x]
        #plot return periods for parents
        boosted_simulation.parent_sim.calculate_and_plot_return_period(ax)
        #plot return periods for boosted
        boosted_simulation.plot_boosted_PI(lead_time_range,ax)
        #fig configs
        ax.set_title(f"Test slice {i+1}: Lead times from {lead_time.start} to {lead_time.stop} days")
        if x%2 == 0:
            ax.set_ylabel("TXx5d [$^\circ$C]")
        else:
            ax.set_ylabel("")
        if i == 1:
            ax.set_xlabel("Return period [year]")
        else:
            ax.set_xlabel("")
        ax.text(0.025,0.92,r"$\textbf{"+string.ascii_lowercase[2*i + x]+r"}$",transform=ax.transAxes)
        pco.set_grid(ax)
handles, labels = ax[0][0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center',ncol=4)
plt.tight_layout()
plt.subplots_adjust(bottom=0.2)
plt.savefig(f"{pco.out_path}Fig4_boosted_return_times.png",bbox_inches="tight")

# === Fig A1 and A2 ===
# lead time by lead time return time calculation

test_slices = ["T1","T2"]
for test_slice in test_slices:
    # figure
    fig,axs = plt.subplots(4,3,figsize = (7,6),sharex=True,sharey=True)
    ax = axs.flatten()
    # get return times for non-boosted
    boosted_simulation = pi.boosted_PI_simulation(test_slice=test_slice)
    #loop over all lead times
    lead_times = range(-18,-6)
    for i,lead_time in enumerate(lead_times):
        #plot return periods for parents
        boosted_simulation.parent_sim.calculate_and_plot_return_period(ax[i])
        #plot return periods for boosted
        boosted_simulation.plot_boosted_PI(lead_time,ax[i])
        #fig config
        ax[i].set_title(f"Lead time: {lead_time} days")
        if i%3 ==0:
            ax[i].set_ylabel("TXx5d [$^\circ$C]")
        else:
            ax[i].set_ylabel("")
        if i > 8: 
            ax[i].set_xlabel("Return period [year]")
        else:
            ax[i].set_xlabel("")
        ax[i].text(0.025,0.85,r"$\textbf{"+string.ascii_lowercase[i]+r"}$",transform=ax[i].transAxes)
        pco.set_grid(ax[i])
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center',ncol=3)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(f"{pco.out_path}FigA{test_slice[1]}_boosted_return_times_lead_by_lead.png",bbox_inches="tight")

# === Fig A3 ===
# Test slice 1 + 2 together (T3)
test_slice = "T3"
fig,axs = plt.subplots(1,2,figsize = (7,3.25),sharex=True,sharey=True)
# calculate return time for each lead time
lead_times = [slice(-18,-13),slice(-12,-7)]
# get return times for non-boosted
boosted_simulation = pi.boosted_PI_simulation(test_slice=test_slice)
for x,lead_time in enumerate(lead_times):
    ax = axs[x]
    #plot return periods for parents
    boosted_simulation.parent_sim.calculate_and_plot_return_period(ax)
    #plot return periods for boosted
    boosted_simulation.plot_boosted_PI(lead_time_range,ax)
    #fig config
    ax.set_title(f"Test slice 1 + 2: Lead times ({lead_time.start}, {lead_time.stop}) days")
    if x%2 == 0:
        ax.set_ylabel("TXx5d [$^\circ$C]")
    else:
        ax.set_ylabel("")
    ax.set_xlabel("Return period [year]")
    ax.text(0.025,0.94,r"$\textbf{"+string.ascii_lowercase[x]+r"}$",transform=ax.transAxes)
    pco.set_grid(ax)
handles, label = ax[0].get_legend_handles_labels()
fig.legend(handles, label, loc='lower center',ncol=3)
plt.tight_layout()
plt.subplots_adjust(bottom=0.32)
plt.savefig(f"{pco.out_path}FigA3_boosted_return_times_{test_slice}.png",bbox_inches="tight")