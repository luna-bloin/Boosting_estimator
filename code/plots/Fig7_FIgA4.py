#================================================================================
# === Fig 7 (+ appendix figure A4): Relationship between parent and offspring ===
#================================================================================

import matplotlib.pyplot as plt
sys.path.append("../utils")
import return_calc as rc
import plot_config as pco
import PI_runs_parent_boosted as pi
import utils as ut
import string

# === Fig 7 ===
f,ax = plt.subplots(1,2,figsize=(7,3.25),sharey=True)
vals_90th_percentile_boost_parent = []
test_slices = ["T1","T2"]
for j, test_slice in enumerate(test_slices): #plot for test slice 1 and 2
    boosted_simulation = pi.boosted_PI_simulation(test_slice=test_slice)
    parents = boosted_simulation.parent_sim.parents
    for i,case in enumerate(parents.year):
        vals = boosted_simulation.TXx5d.sel(start_date=slice(-18,-13)).stack(dim=("start_date","member")).sel(case=int(case))
        # append 90th percentile TXx5d of boosted values, and TXx5d of parent for each parent event
        vals_90th_percentile_boost_parent.append([vals.quantile(0.9).values,case.values])
        #boxplot
        ax[j].boxplot(
            vals.values[~np.isnan(vals.values)],
            positions=[i],
            patch_artist=True, 
            whiskerprops = dict(linewidth=0.75), 
            capprops = dict(linewidth=0.75), 
            flierprops = dict(markersize=3, linewidth=0.75), 
            medianprops=dict(linewidth=0.75,color=pco.colors[0]),
            boxprops=dict(linewidth=0.75,facecolor=color2_with_alpha, edgecolor="k",color=pco.colors[0]),
        )
        # plot parent
        pco.plot_parents(
            ax[j],
            i, 
            boosted_PI_simulation.parent_sim.TXx5d.sel(year=case),
        )
        #fig configs
        if i == 0:
            lab = "Parent event"
            lab_1 = "95$^{\mathrm{th}}$ percentile of boosted simulations"
        else:
            lab = "__no"
            lab_1 = "__no"
        ax[j].set_xticks(range(len(parents[tim])),[f"Parent {i+1}" for i in range(len(parents[tim]))], rotation = 0)
        pco.set_grid(ax)
        ax[j].set_title(f"Test slice {j+1}")
    ax[j].text(0.025,0.94,r"$\textbf{"+string.ascii_lowercase[j]+r"}$",transform=ax[j].transAxes)
ax[0].set_ylabel("TXx5d [$^\circ$C]")
# legend and fig configs
legend_elements = [Patch(facecolor=pco.color2_with_alpha,edgecolor="k", label='Boosted simulations')]
handles, labels = ax[0].get_legend_handles_labels()
combined_handles = handles + legend_elements
f.legend(handles=combined_handles, loc='lower center',ncol=2)
plt.tight_layout()
plt.subplots_adjust(bottom=0.17)
plt.savefig(f"{pco.out_path}Fig7_PI_heatwave_parent_offspring_boxplot.png", bbox_inches="tight")

# == Fig A4: correlation plot ===
f,ax = plt.subplots(figsize=(3.5,3.5),sharey=True)

vals_90th_percentile_boost_parent = np.array(vals_90th_percentile_boost_parent).transpose()
#find correlation
correlation = spearmanr(vals_90th_percentile_boost_parent[0],vals_90th_percentile_boost_parent[1])[0]
ax.plot(
    vals_90th_percentile_boost_parent[1],
    vals_90th_percentile_boost_parent[0],
    ".",
    color=pco.colors[1],
    label=f"Correlation = {correlation:.4}"
)
#plot x=y line
x = np.linspace(8,12)
ax.plot(x,x,label="x=y",color="k")
#fig configs
ax.set_ylabel("$P_{90}$ boosted TXx5d [$^\circ$C]")
ax.set_xlabel("Parent TXx5d [$^\circ$C]")
ax.legend(loc="lower right")
pco.set_grid(ax)
plt.tight_layout()
ax.set_xlim(8.5,12)
ax.set_ylim(8.5,12)
ax.set_xticks(np.arange(9,12.5,0.5))
plt.savefig(f"{pco.out_path}FigA4_correlation_parent_boosted.png", bbox_inches="tight")