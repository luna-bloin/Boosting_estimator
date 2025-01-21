import sys
sys.path.append("../utils")
import PI_runs_parent_boosted as pi


# ====================================================================
# === Find Tref of different time slices and compare to that of T0 ===
# ====================================================================
control_sim = pi.PI_simulation(test_slice="T0").TXx5d
tref_t0 = control_sim.sortby(control_sim,ascending=False)[int(len(control_sim)/10)-1]

for test_slice in ["T1","T2"]:
    test_sim = pi.PI_simulation(test_slice=test_slice)
    tref = test_sim.pref_tref()[0]
    ptref_t0 = control_sim.where(control_sim>=tref,drop=True).count()/control_sim.count()

    print(f"{test_slice}: Tref = {tref:.4}, RE = {((tref-tref_t0)/tref_t0).values:.4}, compared to T0 = {tref_t0.values:.4}")
    print(f"{test_slice}: P_tref = 0.1, RE = {((0.1-ptref_t0)/ptref_t0).values:.4}, compared to T0 = {ptref_t0.values:.4}")















