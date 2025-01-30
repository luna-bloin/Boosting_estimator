import sys
sys.path.append("../utils")
import LE_runs_parent_boosted as le
import return_calc as rc
import ERA5_data as e5
import csv
import numpy as np

prefs = []
LEs = {}
for members_in_LE in [100,30]:
    print(f"LE with {members_in_LE} members")
    LE = le.LE_simulation(tot_members = members_in_LE,anomalies=True)
    LEs[members_in_LE] = LE
    top_events, top_event_info = LE.find_top_TXx5d()
    Tref,P_Tref = LE.pref_tref()
    P_Tref = np.median(P_Tref)
    print(top_event_info)
    # write top event info into file (used to preprocess maps later)
    with open(f"../../inputs/top_le{members_in_LE}.csv", mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        # Write the header
        writer.writerow(['case','peak'])
        # Write the data
        for case, peak in top_event_info.items():
            writer.writerow([case,peak])
    prefs.append(P_Tref)
    # print Tref P_Tref and its relative error in comparison with LE100 (RE will be 0 be design for LE100)
    print(f"Tref = {Tref}, P_Tref = {P_Tref}, RE = {((P_Tref-prefs[0])/prefs[0]):.4}") 

# === Calculate return period estimate of ERA5 2021 PNW heatwave ===
return_level = e5.get_era5_2021_PNW_temperature_anomaly()[0]

# regular large ensembles
for mem in [30,100]:
    print(f"{mem}-member LE")
    ret = rc.find_return_time_naive_gev(LEs[mem].TXx5d_1D_sorted,return_level)

# LE 30, dropping one parent at the time
for case in top_event_info:
    print(f"30-member LE without {case}")
    ret = rc.find_return_time_naive_gev(LEs[30].get_new_sorted_TXx5d_1D(case),return_level)

# boosted estimator
boost = le.boosted_LE_simulation().TXx5d.sel(start_date=slice(-18,-13))
boost_stacked = boost.stack(dim=("case", "start_date", "member"))
print("All boosted simulations")
ret = rc.find_return_time_boost(P_Tref, Tref, boost_stacked,return_level)
for case in boost.case.values:
    print(f"All boosted simulations except {case}")
    boost_stacked = boost.drop_sel(case=case).stack(dim=("case", "start_date", "member"))
    ret = rc.find_return_time_boost(P_Tref, Tref, boost_stacked,return_level)

    