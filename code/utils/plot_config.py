import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import return_calc as rc
import numpy as np

### Matplotlib parameters
plt.rcParams.update({
    "text.usetex": True,
    'axes.linewidth': 0.5,
    "font.size":10,
    'lines.linewidth': 0.75,
    'lines.markersize':  3,
    'patch.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    "figure.dpi": 600,
    "savefig.dpi": 600
})

colors = ["k", ( 112/256, 160/256, 205/256), (196/256, 121/256, 0/256)] #IPCC color scheme
color2_with_alpha = to_rgba(colors[2], alpha=0.5)

#in/out paths
path = '/net/xenon/climphys/lbloin/boost_proba/'
out_path = '../../figures/'

# Function to configure the grid
def set_grid(ax=None):
    """
    Apply consistent grid styling to the specified Axes object
    """
    ax.grid(linestyle='dotted', linewidth=0.3)

# Plot that highlights the parents 
def plot_parents(ax,x,y,label=True):
    if label == True:
        lab = "Parent Events"
    else:
        lab = "__nolab"
    ax.plot(
        x,
        y,
        "D",
        mec=colors[0],
        mew=0.5,
        color=colors[1],
        markersize=4,
        label=lab,
        zorder=3
    )
    return None

def plot_highest_parent(ax,x,y,label,color=colors[0]):
    ax.plot(
        x,
        y,
        "*",
        mec=colors[0],
        color=color,
        label=label,
        markersize=6,
        mew=0.5,
        zorder=3
    )

# plotting the return time of a non-boosted simulation, with GEV fit
def plot_return(TXx5d,bootstrapped_GEV,label,ax,color = colors[1],msize=3):
    # plot the naively estimated return time
    ax.plot(
        rc.naive_estimator(TXx5d)[0], 
        rc.naive_estimator(TXx5d)[1],
        ".",
        label=label,
        markersize=msize,
        color=color
    )
    #plot bootstrapped GEV return time estimation
    # median
    bootstrapped_GEV.median("bootstrap").plot(color=color,ax=ax)
    # upper and lower quantile
    x = bootstrapped_GEV.return_time
    quant_low = bootstrapped_GEV.quantile(0.025,"bootstrap")
    quant_high = bootstrapped_GEV.quantile(0.975,"bootstrap")
    ax.fill_between(x,quant_low,quant_high,alpha=0.3,color=color)
    return None

def convert_ticklabels_to_strings(f,only_y=False,scientificx=False,scientificy=False):
    """
    Convert xticklabels and yticklabels to strings for each subplot in figure f.
    
    Parameters:
    f (matplotlib.figure.Figure): The figure containing subplots.
    """
                
    
    for ax in f.get_axes():  # Iterate over all axes in the figure
        if only_y == False:
            if ax.get_xticklabels():  # Check if xticklabels exist
                ax.set_xticklabels([f"10$^{{{int(np.log10(label))}}}$" if scientificx == True and label > 0 else f"{label:.6g}" for label in ax.get_xticks()])
        if ax.get_yticklabels():  # Check if yticklabels exist
            ax.set_yticklabels([f"10$^{{{int(np.log10(label))}}}$" if scientificy == True and label > 0 else f"{label:.6g}" for label in ax.get_yticks()])

def convert_colorbar_ticks_to_strings(cbar):
    """
    Convert colorbar tick labels to to strings for each subplot in figure f.
    
    Parameters:
    cbar (matplotlib.colorbar.Colorbar): The colorbar whose ticks need formatting.
    """
    ticks = cbar.get_ticks()
    formatted_ticks = [f"{tick:.6g}" for tick in ticks]
    cbar.set_ticks(ticks)  # Ensure the ticks remain the same
    cbar.set_ticklabels(formatted_ticks)  # Apply formatted labels