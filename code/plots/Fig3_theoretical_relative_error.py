# ====================================================================================================================================================================================================================
# === Code necessary to reproduce Figure 3: Theoretical relative error of the boosting estimator with N=50 and (solid orange line) N_b = 100, (dashed orange line) N_b = 500  and (dotted orange line) N_b = 3000. ===
# ====================================================================================================================================================================================================================

import numpy  as np
import matplotlib.pyplot as plt
import string
import sys
sys.path.append("../utils")
import plot_config as pco

## Additional plot configs
plt.rcParams['xtick.top'] = False
plt.rcParams['ytick.right'] = False
plt.rcParams['text.latex.preamble']= r'\usepackage{amsfonts}'
plt.rcParams['lines.linewidth'] = 0.75
linestyles = ["-","--","dotted"]

## Relative error function
def compute_relative_error_ensemble_boosting_probability_estimator( N , p_Tref_ini , N_b , p_Tref_boost , p_Text_boost ):
    ## Description of parameters
    # N is the number of days in the initial long run
    # p_Tref_ini is P(T>T_ref) in the initial ensemble
    # N_b is the number of members in the boosted ensemble
    # p_Tref_boost is P(T>T_ref) in the boosted ensemble
    # p_Text_boost is P(T>T_ext) in the boosted ensemble (must be below p_Tref_boost)
    
    ## Elements to compute the variance
    E_p_Tref_ini_square = N * p_Tref_ini * ( 1 + ( N - 1 ) * p_Tref_ini ) / N**2
    E_N_square = N_b * p_Text_boost * ( 1 + ( N_b - 1 ) * p_Text_boost )
    E_D_square = N_b * p_Tref_boost * ( 1 + ( N_b - 1 ) * p_Tref_boost )
    E_D_fourth = N_b * p_Tref_boost  + 4 * N_b * ( N_b - 1 ) * p_Tref_boost**2  +  6 * N_b * ( N_b - 1 ) * ( N_b - 2 ) * p_Tref_boost**3  +  ( N_b**4 - 6*N_b*(N_b-1)*(N_b-2) - 4*N_b*(N_b-1) - N_b ) * p_Tref_boost**4
    V_D_square = E_D_fourth - E_D_square**2
    E_ND_square = N_b * p_Text_boost  + 2 * N_b * ( N_b - 1 ) * p_Text_boost * p_Tref_boost  +  2 * N_b * ( N_b - 1 ) * p_Text_boost**2  +  5 * N_b * ( N_b - 1 ) * ( N_b - 2 ) * p_Tref_boost * p_Text_boost**2  +  N_b * ( N_b - 1 ) * ( N_b - 2 ) * p_Tref_boost**2 * p_Text_boost  +  ( N_b**4 - 6*N_b*(N_b-1)*(N_b-2) - 4*N_b*(N_b-1) - N_b ) * p_Tref_boost**2 * p_Text_boost**2
    Cov_ND_square = E_ND_square - E_N_square * E_D_square
    E_p_Text_ini = p_Text_boost * p_Tref_ini / p_Tref_boost

    ## 
    E_p_Text_ini_square = E_p_Tref_ini_square * E_N_square / E_D_square * ( 1 - Cov_ND_square / ( E_N_square * E_D_square ) + V_D_square / E_D_square**2 )
    V_p_Text_ini = E_p_Text_ini_square - E_p_Text_ini**2
    
    return np.sqrt(V_p_Text_ini) / E_p_Text_ini

## Plotting results
f,axs=plt.subplots(1,2,figsize=(7,4),sharey=True,sharex=True)
N = 50
for i,p_Tref_boost in enumerate([0.75,0.3]): # probability to be above T_ref in the boosted ensemble (P(Tref AC) )
    ax = axs[i]
    # define P(Tref), calculate P(Text AC) and P(Text)
    p_Tref_ini =0.1
    p_Text_boost = np.linspace( 0 , p_Tref_boost , 200 ) # probability to be above T_ext in the boosted ensemble (must be below p_Tref_boost)
    p_Text_ini = p_Tref_ini * p_Text_boost / p_Tref_boost # P(T>T_ext) in the initial ensemble (ie climatological probability)

    # loop over three conifguration of N_b the number of boosted simulations
    for x,N_b in enumerate([100,500,3000]):
        N = 50 # for boosting estimator
        # calculate error for boosted estimator and plotting it
        relative_error_boost = compute_relative_error_ensemble_boosting_probability_estimator( 
            N,
            p_Tref_ini, 
            N_b, 
            p_Tref_boost, 
            p_Text_boost 
        )
        ax.plot( 
            p_Text_ini, 
            relative_error_boost * 100, 
            label=fr"Boosting, $N=$ 50, $N_b$ = {(N_b)}",
            color=pco.colors[2],
            linestyle=linestyles[x]
        )
        # calculate error for naive estimator with equivalent computational resource use
        N_equiv = (N+(N_b)/(100/21))
        relative_error_naive_estimator = np.sqrt( p_Text_ini * ( 1 - p_Text_ini ) / N_equiv ) / p_Text_ini
        ax.plot( 
            p_Text_ini, 
            relative_error_naive_estimator * 100, 
            color=pco.colors[1], 
            label=fr"Naive, $N$= {(N)}+{(int(N_b/(100/21)))}",
            linestyle=linestyles[x]
        )
        # print naive estimator with N = 50 and N= 4000 (placed here so they appear correctly on the legend)
        if x < 2:
            #full PI-control simulation
            if x == 1:
                N = 4000
            elif x ==0:
                N=50
            relative_error_naive_estimator = np.sqrt( p_Text_ini * ( 1 - p_Text_ini ) / N ) / p_Text_ini
            ax.plot( 
                p_Text_ini, 
                relative_error_naive_estimator * 100, 
                color=pco.colors[0], 
                label=fr"Naive, $N=$ {N}",
                linestyle=linestyles[x]
            )
    # figure parameters     
    if i == 0:
        ax.set_ylabel('Relative error [\%]')
    ax.set_xlabel('$\hat{p}_{T \geq T_{\mathrm{ext}}}$')    
    ax.set_title(r"$\hat{p}_{T\geq T_{\mathrm{ref}} | \mathrm{AC}_t^{\epsilon}}$="+ str(p_Tref_boost))
    ax.text(0.025,0.92,r"\textbf{"+string.ascii_lowercase[i]+r"}",transform=ax.transAxes)
    ax.set_ylim( 0 , 105 )
    ax.set_xlim( 5*0.0001 ,0.1)
    pco.set_grid(ax)
    ax.set_xscale('log')

handles, labels = axs[0].get_legend_handles_labels()
f.legend(handles, labels, loc='lower center',ncol=3)
plt.tight_layout()
plt.subplots_adjust(bottom=0.3)
pco.convert_ticklabels_to_strings(f,scientificx=True)
plt.savefig(f"{pco.out_path}Fig3_theoretical_relative_error.png",bbox_inches="tight")