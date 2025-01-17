import xarray as xr

def multi_to_single_index(ds,dims=("case","lead_time"),new_name="lead_ID"):
    stacked = ds.stack(event=dims)
    coords = [str(st) for st in stacked["event"].values]
    stacked = stacked.reset_index("event")
    stacked["event"] = coords
    stacked = stacked.drop_vars(dims)
    stacked = stacked.rename({"event":new_name})
    return stacked