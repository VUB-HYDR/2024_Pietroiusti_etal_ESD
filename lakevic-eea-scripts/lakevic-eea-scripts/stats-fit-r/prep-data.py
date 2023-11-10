# -*- coding: utf-8 -*-
"""
Created on Fri May  5 12:15:21 2023

@author: rpietroi

Clean and prep messy data for use in R stat-fits
"""


import sys, os 
import numpy as np
import statsmodels.api as sm
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr
import glob

# set working directory
os.chdir(os.path.dirname(__file__)) #os.getcwd() to check

#%%

"""
1. GMST data obs from ncfile to csv for easier use in Py/R
"""
#gmst pbs data 
gmst_path = os.path.join('..', '..', "data/input-data/gmst-obs/igiss_al_gl_a_4yrlo.nc") 

# clean
gmst_raw = xr.open_dataarray(gmst_path, decode_times=False)
gmst_raw = gmst_raw.assign_coords(time=np.arange(1879, 1879+len(gmst_raw)))
gmst_raw = gmst_raw.sel(time=slice(1879, 2022))
gmst_raw = gmst_raw.rename({'time': 'year'}).to_dataframe().rename(columns={'Ta':'gmst'})
gmst_raw['gmst'].plot(label = 'gmst anomaly wrt 1951-1980 gisstemp (deg C)')

# save
#gmst_raw.to_csv(os.path.join(gmst_path,'..', 'igiss_al_gl_a_4yrlo.csv'), index=True)

#%%

"""
2. IOD obs: detrended and not detrended (wrt time or should i do it wrt to sst or gmst??)
"""
#path
iod_path = os.path.join('..', '..', "data/data-modified/iod-obs/dmi_seasonal_1870_2021.csv") 
iod_raw = pd.read_csv(iod_path, index_col=0)
iod_ond = iod_raw.iloc[:,0]
#iod_ond.plot()

def detrend(year, data):
         fit = sm.OLS(data, sm.add_constant(year)).fit()
         return fit.resid

def scale(data):
         return (data - data.mean()) / data.std()

iod_ond_detrend = detrend(iod_ond.index, iod_ond)
iod_ond_detrend.plot()

# save
#iod_ond.to_csv(os.path.join(iod_path,'..', 'dmi_ond_1870_2021.csv'), index=True)
#iod_ond_detrend.to_csv(os.path.join(iod_path,'..', 'dmi_ond_detrended_1870_2021.csv'), index=True)

#%%

"""
3. GMST models: 4-yr smoothed
"""
#path
in_path = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\data\input-data\gmst-models'
out_path = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\data\data-modified\gmst-models'

models = ['CanESM5','CNRM-CM6-1','GFDL-ESM4','IPSL-CM6A-LR','MIROC6','MRI-ESM2-0']
scenarios = ['historical', 'ssp370', 'hist-nat']

# hist-nat
fig, ax = plt.subplots()

scenario = scenarios[2] 
for i in models:
    gmst_path = glob.glob(os.path.join(in_path, scenario, i, '*.txt'))[0]
    gmst = pd.read_csv(gmst_path, delim_whitespace=True, comment='#', header=None, names=['years','gmst'], index_col=0)
    gmst_smo = gmst.rolling(4, min_periods=1, center=True).mean()
    
    #as anomaly wrt 1951-1980
    baseline = gmst.loc[1951:1980].mean()
    gmst_ano = gmst - baseline
    gmst_ano_smo = gmst_ano.rolling(4, min_periods=1, center=True).mean()
     
    # save
    #gmst_smo.to_csv(os.path.join(out_path, scenario, 'gmst_{}_hist-nat_{}_{}_4yrlo.csv'.format(i, gmst_smo.index[0], gmst_smo.index[-1])), index=True, float_format='%.3f')
    #gmst_ano_smo.to_csv(os.path.join(out_path, scenario, 'wrt1951-80', 'gmst_{}_hist-nat_{}_{}_wrt1951-80_4yrlo.csv'.format(i, gmst_smo.index[0], gmst_smo.index[-1])), index=True, float_format='%.3f')

    # as anomaly wrt PI
    baseline = gmst.loc[1850:1900].mean()
    gmst_ano = gmst - baseline
    gmst_ano_smo = gmst_ano.rolling(4, min_periods=1, center=True).mean()
    
    gmst_ano_smo['gmst'].plot(ax=ax, label=i)
    
    # save
    #gmst_ano_smo.to_csv(os.path.join(out_path, scenario, 'wrt1850-1900', 'gmst_{}_hist-nat_{}_{}_wrt1850-1900_4yrlo.csv'.format(i, gmst_smo.index[0], gmst_smo.index[-1])), index=True, float_format='%.3f')

    
plt.legend()

#%%
# hist-rcp370 to 2100 and to 2020
fig, ax = plt.subplots()

for i in models:
    gmst_path_a = glob.glob(os.path.join(in_path, scenarios[0], i, '*.txt'))[0]
    gmst_path_b = glob.glob(os.path.join(in_path, scenarios[1], i, '*.txt'))[0]
    print(gmst_path_a, gmst_path_b)
    gmst_a = pd.read_csv(gmst_path_a, delim_whitespace=True, comment='#', header=None, names=['years','gmst'], index_col=0)
    gmst_b = pd.read_csv(gmst_path_b, delim_whitespace=True, comment='#', header=None, names=['years','gmst'], index_col=0)
    gmst = pd.concat([gmst_a, gmst_b])
    gmst_smo = gmst.rolling(4, min_periods=1, center=True).mean()
    
    # as anomaly wrt preindustrial
    baseline = gmst.loc[1850:1900].mean()
    gmst_ano = gmst - baseline
    gmst_ano_smo = gmst_ano.rolling(4, min_periods=1, center=True).mean()
    
    # plot
    #gmst.plot()
    #gmst_smo.plot()
    #gmst_ano['gmst'].plot(ax=ax, label=i)
    gmst_ano_smo['gmst'].plot(ax=ax, label=i)
    
    # save
    #gmst_smo.to_csv(os.path.join(out_path, 'hist-rcp370', 'gmst_{}_hist-rcp370_{}_{}_4yrlo.csv'.format(i, gmst_smo.index[0], gmst_smo.index[-1])), index=True, float_format='%.3f')
    #gmst_ano_smo.to_csv(os.path.join(out_path, 'hist-rcp370', 'wrt1850-1900', 'gmst_{}_hist-rcp370_{}_{}_wrt1850-1900_4yrlo.csv'.format(i, gmst_smo.index[0], gmst_smo.index[-1])), index=True, float_format='%.3f')
    #gmst_smo.loc[1850:2020].to_csv(os.path.join(out_path, 'hist-rcp370', 'gmst_{}_hist-rcp370_{}_{}_4yrlo.csv'.format(i, gmst_smo.loc[1850:2020].index[0], gmst_smo.loc[1850:2020].index[-1])), index=True, float_format='%.3f')
    
plt.legend()



#%%

"""
4. IOD models: detrended and not detrended - get from KX
"""
