import cartopy.crs as ccrs
import matplotlib.patches as patches
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.append("../utils")
import plot_config as pco
import string
import ERA5_data as e5

def add_PNW_box(ax):
    PNW_rect = patches.Rectangle((237-180, 45), 4, 7, linewidth=1, edgecolor='k', facecolor='none')
    ax.add_patch(PNW_rect)

def around_PNW(ds):
    return ds.sel(lat=slice(20,75),lon=slice(150,350))

def plot_temp(to_plot,col,col_wrap,figsize):
    f = to_plot.Tx5d_anom.plot(x="lon", y="lat",col=col,figsize=figsize,transform=ccrs.PlateCarree(),subplot_kws={'projection':proj, "aspect": 1.5},col_wrap=col_wrap,vmin=-15,vmax=15,cmap="RdBu_r")
    f.cbar.set_label("Tx5d anomaly [$^\circ$C]", size=10)# Set label and font size
    return f

def plot_z500(ax,min_all,max_all,to_plot):
    g = to_plot.Z500.plot.contour(ax=ax,vmin=min_all,vmax=max_all,levels=[min_all +x*50 for x in range(int((max_all-min_all)/50)+1)],cmap="k",transform=ccrs.PlateCarree(),linewidth=0.5)
    ax.clabel(g, inline=True, fontsize=6,levels=g.levels[::3])

# projection for maps
proj = ccrs.PlateCarree(central_longitude=180)
#open ERA5 map
era5_map = e5.get_era5_map_at_peak()

title_E = {30: "E", 100:r"\tilde{E}"} # so that for LE 30 you get E1 etc and for LE100 you get tilde E1 etc
fig = {30:"9",100:"A8"} #for saving

for x,LE_mem in enumerate([30,100]):
    # open CESM parent maps
    LE_map = xr.open_dataset(f"{pco.path}maps_of_tops_LE{LE_mem}.nc")
    era5_map["lat"] = LE_map["lat"] #there are small (10^-10) differences in lat definition between cesm2 and era5
    # ========================
    # === Fig 9 and Fig A8 ===
    # ========================
    titles = [fr"${title_E[LE_mem]}1$","ERA5: 2021 PNW heatwave",fr"${title_E[LE_mem]}2-{title_E[LE_mem]}3$ mean",fr"${title_E[LE_mem]}2-{title_E[LE_mem]}13$ mean"]
    # concat maps to plot
    to_plot = around_PNW(xr.concat([
        LE_map.isel(case=0), # E1
        era5_map, #ERA5 data
        LE_map.isel(case=[1,2]).mean("case"), #E2-3 mean
        LE_map.isel(case=slice(1,13)).mean("case") #E2-13 mean
    ],dim= "events"))
    # find min and max Z500 values 
    min_all = int(round(to_plot.Z500.min().values.item(),-2))
    max_all = int(round(to_plot.Z500.max().values.item(),-2))
    #plot temperature anomalies
    f = plot_temp(to_plot,"events",2,(7,3))
    for i,a in enumerate(f.axs.flat):
        # plot Z500
        plot_z500(a,min_all,max_all,to_plot.isel(events=i))
        #plot PNW region box
        add_PNW_box(a)
        #fig configs
        a.set_title(titles[i])
        a.coastlines()
        a.text(152-180,77,r"$\textbf{"+string.ascii_lowercase[i]+r"}$")
    pco.convert_colorbar_ticks_to_strings(f.cbar)
    plt.savefig(f"{pco.out_path}Fig{fig[LE_mem]}_maps.png",bbox_inches="tight")

    # =====================
    # === Fig A7 and A9 ===
    # =====================
    to_plot = around_PNW(LE_map) #plot all top 13 maps
    #plot temperature anomalies
    f = plot_temp(to_plot,"case",3,(7,4.5)) 
    for i,a in enumerate(f.axs.flat):
        if i ==13:
            break
        # plot Z500
        plot_z500(a,min_all,max_all,to_plot.isel(case=i))
        #plot PNW region box
        add_PNW_box(a)
        a.set_title(fr"${title_E[LE_mem]}{i+1}$")
        a.coastlines()
        a.set_xlim(150-180,350-180)
        a.set_ylim(20,75)
        a.text(152-180,77,r"$\textbf{"+string.ascii_lowercase[i]+r"}$")
    pco.convert_colorbar_ticks_to_strings(f.cbar)
    plt.savefig(f"{pco.out_path}FigA{7+2*x}_maps_case_by_case.png",bbox_inches="tight")