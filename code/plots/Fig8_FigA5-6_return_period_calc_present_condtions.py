# =========================================================================================================================================================
# === Fig * (+ appendix figure 5-6): Return periods of TXx5d in the PNW region estimated with the boosting estimator under current climate conditions.  ===
# =========================================================================================================================================================

import matplotlib.pyplot as plt
import sys
sys.path.append("../utils")
import return_calc as rc
import LE_runs_parent_boosted as le
import plot_config as pco
import string
import ERA5_data as e5
# === Open necessary data ===
#100 member LE

LE_100 = le.LE_simulation(100,anomalies=True)
PNW_2021_ERA5_anomaly = e5.get_era5_2021_PNW_temperature_anomaly()[0]



# # === Fig 8 ===
# fig,axs = plt.subplots(1,2,figsize = (7,3.25),sharex=True,sharey=True)
# # calculate return time pooled lead times
# lead_time_range = slice(-18,-13)
# tit = ["All parent events", "Most extreme parent event removed"] 

# for i,drop_case in enumerate([[],["10:2007"]]):
#     # get boosted simulation
#     boosted_simulation = le.boosted_LE_simulation(drop_case=drop_case)
#     ax = axs[i]
#     #plot return periods for parents
#     LE_100.calculate_and_plot_return_period(ax,color=pco.colors[0])
#     # plot ERA5 anomaly
#     ax.axhline(PNW_2021_ERA5_anomaly.values,linestyle="--",color=pco.colors[0], label="ERA5: 2021 PNW heatwave")
#     boosted_simulation.parent_sim.calculate_and_plot_return_period(ax, remove_a_parent=drop_case)
#     #plot return periods for boosted
#     boosted_simulation.plot_boosted_LE(lead_time_range,ax,era5_anomaly=PNW_2021_ERA5_anomaly.values)
#     #fig configs
#     ax.set_title(tit[i])
#     if i == 0:
#         ax.set_ylabel("TXx5d [$^\circ$C]")
#     else:
#         ax.set_ylabel("")
#     ax.set_xlabel("Return period [year]")
#     ax.text(0.025,0.92,r"$\textbf{"+string.ascii_lowercase[i]+r"}$",transform=ax.transAxes)
#     pco.set_grid(ax)
#     ax.set_xticks([1,10,100,1000,10000,100000])
#     ax.set_yticks([8,10,12,14,16])


# handles, labels = axs[0].get_legend_handles_labels()
# fig.legend(handles, labels, loc='lower center',ncol=3)
# pco.convert_ticklabels_to_strings(fig,scientificx=True)
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.35)
# plt.savefig(f"{pco.out_path}Fig8_boosted_return_time_LE.png", bbox_inches="tight")



# # === Fig A5 ===
# # lead time by lead time return time calculation
# fig,axs = plt.subplots(4,3,figsize = (7,6),sharex=True,sharey=True)
# ax = axs.flatten()
# # get return times for non-boosted
# boosted_simulation = le.boosted_LE_simulation()
# #loop over all lead times
# lead_times = range(-18,-6)
# for i,lead_time in enumerate(lead_times):
#     #plot return periods for parents
#     LE_100.calculate_and_plot_return_period(ax[i],color=pco.colors[0])
#     # plot ERA5 anomaly
#     ax[i].axhline(PNW_2021_ERA5_anomaly.values,linestyle="--",color=pco.colors[0], label="ERA5: 2021 PNW heatwave")
#     boosted_simulation.parent_sim.calculate_and_plot_return_period(ax[i])
#     #plot return periods for boosted
#     boosted_simulation.plot_boosted_LE(lead_time,ax[i],era5_anomaly=PNW_2021_ERA5_anomaly.values)
#     #fig config
#     ax[i].set_title(f"Lead time: {lead_time} days")
#     if i%3 ==0:
#         ax[i].set_ylabel("TXx5d [$^\circ$C]")
#     else:
#         ax[i].set_ylabel("")
#     if i > 8: 
#         ax[i].set_xlabel("Return period [year]")
#     else:
#         ax[i].set_xlabel("")
#     ax[i].text(0.025,0.85,r"$\textbf{"+string.ascii_lowercase[i]+r"}$",transform=ax[i].transAxes)
#     pco.set_grid(ax[i])
#     ax[i].set_xticks([1,10,100,1000,10000,100000])
#     ax[i].set_yticks([8,10,12,14,16])
# handles, labels = ax[0].get_legend_handles_labels()
# fig.legend(handles, labels, loc='lower center',ncol=3)

# pco.convert_ticklabels_to_strings(fig,scientificx=True)
# plt.tight_layout()
# plt.subplots_adjust(bottom=0.2)
# plt.savefig(f"{pco.out_path}FigA5_boosted_return_times_lead_by_lead_LE.png",bbox_inches="tight")



# === Fig A6 ===
print("plotting A6")
fig,axs = plt.subplots(4,2,figsize = (7,7.5),sharex=True,sharey=True)
ax = axs.flatten()
# calculate return time this lead time range
lead_time_range = slice(-18,-13)

# get a list of top parents sorted by intensity
boosted_simulation = le.boosted_LE_simulation()
parents = boosted_simulation.parent_sim.find_parents()
parents = parents.sortby(parents, ascending=False)
parent_list = [[f"{el.member.item():02d}:{el.dim.item()}"] for el in parents]
parent_list.insert(0,[])

for i,drop_case in enumerate(parent_list):
    # get boosted simulation
    boosted_simulation = le.boosted_LE_simulation(drop_case=drop_case)
    #plot return periods for parents
    LE_100.calculate_and_plot_return_period(ax[i],color=pco.colors[0])
    # plot ERA5 anomaly
    ax[i].axhline(PNW_2021_ERA5_anomaly.values,linestyle="--",color=pco.colors[0], label="ERA5: 2021 PNW heatwave")
    boosted_simulation.parent_sim.calculate_and_plot_return_period(ax[i], remove_a_parent=drop_case)
    #plot return periods for boosted
    boosted_simulation.plot_boosted_LE(lead_time_range,ax[i],era5_anomaly=PNW_2021_ERA5_anomaly.values)
    # find rank of parent - for title
    if i > 0:
        mem = int(drop_case[0][0:-5])
        year = int(drop_case[0][-4:])
        top_events = boosted_simulation.parent_sim.find_top_TXx5d()[0]
        parent_val = top_events.where(top_events.member == mem,drop=True).where(top_events.dim == year,drop=True).values
        parent_rank = len(top_events.where(top_events >= parent_val,drop=True))
    #fig configs
    if i%2 == 0:
        ax[i].set_ylabel("TXx5d [$^\circ$C]")
    else:
        ax[i].set_ylabel("")
    if i > 5:
        ax[i].set_xlabel("Return period [year]")
    else:
        ax[i].set_xlabel("")
    ax[i].text(0.025,0.9,r"$\textbf{"+string.ascii_lowercase[i]+r"}$",transform=ax[i].transAxes)
    pco.set_grid(ax[i])
    ax[i].set_xticks([1,10,100,1000,10000,100000])
    ax[i].set_yticks([8,10,12,14,16])
    if i == 0:
        tit = "All parent events"
    else:
        tit = f"E{parent_rank} removed"
    ax[i].set_title(f"{tit}")

handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center',ncol=3)
pco.convert_ticklabels_to_strings(fig,scientificx=True)
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.savefig(f"{pco.out_path}FigA6_return_time_removing_LE.png", bbox_inches="tight")


    


    

    

    








