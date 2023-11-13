# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:13:20 2023

@author: rpietroi

WBM evaluation 
Analysing Python WBM output : look at BIAS cfr to OBS


"""
#%%


#=============#
#== IMPORTS ==#
#=============#

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os , glob
import datetime
from datetime import timedelta
import csv
import math
import pandas as pd
from matplotlib import transforms
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdates
import seaborn as sns
import statsmodels.api as sm


#%%

#===========#
#== DATA  ==#
#===========#

workDIR = os.path.join(r"C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea" )

inDIR_Py = os.path.join(workDIR,r'WBM-git\output\2022\output\obs\PERSIANN\v01r01_out\\') 
inDIR_Obs = os.path.join(workDIR,r'lakevic-eea-scripts\data\data-modified\lakelevels' ) 
fig_path = os.path.join(workDIR,r'lakevic-eea-scripts\figures\figures_5_mar23_wbmeval') 


# Python version
#=============================================================================
filepath = os.path.join(inDIR_Py, 'WBM_run_PERSIANN_obs_v01r01_1983_2020_mm.csv') 
WBM_obs_raw = pd.read_csv(filepath)

filepath = os.path.join(workDIR,r'lakevic-eea-scripts\output\event_def_dt180\obs-WBM\PERSIANN\out\dLdt_1983_2020_dt180_PERSIANN_obs_out.csv')
WBM_obs_dLdt = pd.read_csv(filepath, parse_dates=['date'], index_col=('date')) # with semi-obs outflow


# Observational lakelevels 
#=============================================================================
filepath = os.path.join(inDIR_Obs, 'lakelevel_all_intr_DH22.csv') 
lakelevels_obs_raw = pd.read_csv(filepath)


filepath = os.path.join(workDIR,r'lakevic-eea-scripts\output\event_def_dt180\obs\dLdt_1948_2022_dt180_obs.csv')
LL_obs_dLdt = pd.read_csv(filepath, parse_dates=['date'], index_col=('date'))

# Data  
#=============================================================================
name_obs = 'WBM obs (1983-2020)'
name_obs_plot = 'WBM observational run \n(1983-2020)'

#WBM_obs
df = WBM_obs_raw[['date','L_obs', 'P_lake', 'Q_in', 'E_lake', 'Q_out']]
df['date'] = pd.to_datetime(df['date'])
WBM_obs = df.set_index(['date'])

#LL_observed
df = lakelevels_obs_raw
df['date'] = pd.to_datetime(df['date'])
LL_obs = df.set_index(['date'])

startYEAR = WBM_obs.index[0].year
endYEAR = WBM_obs.index[len(WBM_obs)-1].year

# bias
left = WBM_obs.index[0]
right = WBM_obs.index[len(WBM_obs)-1]

residual = WBM_obs.loc[str(left):str(right)]['L_obs'] - LL_obs.loc[str(left):str(right)]['water_level']
df_res = pd.DataFrame(residual)
df_res.columns = ['residual']


#%% Calculations, analysis of bias 

print(df_res.min(), df_res.max())
print(df_res[df_res['residual'] == df_res['residual'].min()])
print(df_res[df_res['residual'] == df_res['residual'].max()])

#%%
print('record-breaking 2020 \n', LL_obs[LL_obs['water_level'] == LL_obs['water_level'].max()])

print(WBM_obs[WBM_obs['L_obs'] == WBM_obs['L_obs'].max()])

print('\n observed november\n', LL_obs.loc['2019-11-19'])

print('\n model may 17 \n', WBM_obs.loc['2020-05-17'], '\n')
print('model november \n ', WBM_obs.loc['2019-11-19'], '\n')

print('180 days before =', datetime.datetime.strptime('2020-05-17', "%Y-%m-%d") - timedelta(days=180))
print('nov-may model', WBM_obs['L_obs'].loc['2020-05-17'] - WBM_obs['L_obs'].loc['2019-11-19'])
print('nov - may obs', LL_obs['water_level'].loc['2020-05-17'] - LL_obs['water_level'].loc['2019-11-19'])



print('\n365 days before =', datetime.datetime.strptime('2020-05-17', "%Y-%m-%d") - timedelta(days=365))
print('may - may model ', WBM_obs['L_obs'].loc['2020-05-17'] - WBM_obs['L_obs'].loc['2019-05-18'])
print('may - may obs', LL_obs['water_level'].loc['2020-05-17'] - LL_obs['water_level'].loc['2019-05-18'])


print('\nmay19 - nov19 model ', WBM_obs['L_obs'].loc['2019-11-19'] - WBM_obs['L_obs'].loc['2019-05-18'])
print('may19 - nov19 obs', LL_obs['water_level'].loc['2019-11-19'] - LL_obs['water_level'].loc['2019-05-18'])


#%%

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

#%%

predictions = WBM_obs.loc[str(left):str(right)]['L_obs'] 
targets = LL_obs.loc[str(left):str(right)]['water_level']

print('rmse is', rmse(predictions,targets), 'm')

print('mean bias', np.mean(residual), 'm')
# %%

# review phase 

left= 2007
right= 2007  # to get nicer labels 
print(f'mean bias {left}-{right} =', residual.loc['{}-01-01'.format(left):'{}-12-31'.format(right)].mean())





# %%
