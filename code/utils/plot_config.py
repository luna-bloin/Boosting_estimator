import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import return_calc as rc

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
    or the current Axes if none is provided.
    """
    if ax is None:
        ax = plt.gca()  # Get current Axes
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

# plotting the return time of a non-boosted simulation, with GEV fit
def plot_return(TXx5d,bootstrapped_GEV,label,ax,color = colors[1],msize=3):
    # plot the naively estimated return time
    ax.plot(
        rc.naive_estimator(TXx5d)[0], 
        rc.naive_estimator(TXx5d)[1],
        ".",
        label=label,
        markersize=5,
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
