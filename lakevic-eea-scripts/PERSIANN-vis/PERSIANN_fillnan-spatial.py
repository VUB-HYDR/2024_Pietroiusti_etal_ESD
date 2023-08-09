"""
rosa.pietroiusti@vub.be
Pietroiusti et al 2024 ESD

Fill Persiann nan data outputting a spatial file by replacing missing days with climatology for calendar day
"""

#%%
import sys, os 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr


# set working directory
os.chdir(os.path.dirname(__file__)) #os. getcwd() to check

#set data-path
data_path = os.path.join('..', '..', "data") 

#set figures-path
fig_path = os.path.join('..', '..', "figures/figures_jun23_fillnan") 
if os.path.exists(fig_path):
    print("The directory", fig_path, "exists!")
else:
    os.makedirs(fig_path)
    print("The directory", fig_path, "was made!")


out_path = os.path.join(data_path, 'data-modified', 'PERSIANN')


#%%

ncfile = os.path.join(out_path, 'PERSIANN-CDR_v01r01_1983_2020_remapped_owngrid_reorder_timelatlon.nc')
var = 'precipitation'

# Open dataset and get key information
with xr.open_dataset(ncfile, engine='netcdf4') as ds:
    print(str(list(ds.keys())))
    print(str(list(ds.dims)))
    ds_dims = list(ds.dims) # time,lon, lat
    ds_varnames = list(ds.keys()) # get names of variables precipitation
    time = np.array(ds.get('time'))
    PERSIANN_raw = ds[str(var)] 



# %%

# find nan days 

slice = PERSIANN_raw.sel(time=slice("1983-01-01", "1983-12-31"))
test = slice.isnull()



# %%
print(test.values.sum() / (130*130))

# 10 days in 1983

# %%
# test
indices = []
for t in slice.time:
    if np.isnan(slice.sel(time=t).values).all() == True:
        indices.append(t.values)

# %%
# do for whole dataset 

indices = []
data = PERSIANN_raw
for t in data.time:
    if np.isnan(data.sel(time=t).values).all() == True:
        indices.append(t.values)
print(len(indices))

# %%

climatology_daily = PERSIANN_raw.groupby("time.dayofyear").mean('time')

# %%

PERSIANN_fillnan_climato = PERSIANN_raw.copy()
for i in indices:
    print(i)
    PERSIANN_fillnan_climato.loc[dict(time=i)] = climatology_daily[PERSIANN_fillnan_climato.sel(time=i).time.dt.dayofyear.values - 1].values



#%%
# check if there are nans left 

check = []
data = PERSIANN_fillnan_climato
for t in data.time:
    if np.isnan(data.sel(time=t).values).all() == True:
        check.append(t.values)
print(len(check))


#%%

# check yearly accumulation




# %%

# save
#PERSIANN_fillnan_climato.to_netcdf(os.path.join(out_path, 'PERSIANN-CDR_v01r01_1983_2020_remapped_owngrid_reorder_timelatlon_fillnan_climato.nc' ))


# %%
yearly_tseries_raw = PERSIANN_raw.mean(dim=['lat','lon']).resample(time='Y', label='left').sum()
yearly_tseries_fillnan = PERSIANN_fillnan_climato.mean(dim=['lat','lon']).resample(time='Y', label='left').sum()
yearly_tseries_fillnan.plot(label='fillnan')
yearly_tseries_raw.plot(label='raw')
plt.legend()
plt.ylabel("precipitation (mm/year)")

#plt.savefig(os.path.join(fig_path,'PERSIANN_spatial_fillnan_fulldomain.png'),dpi=300)

# %%

# mask over basin and check yearly accumuation in PERSIANN_anomalytseries_LVB_analysis
