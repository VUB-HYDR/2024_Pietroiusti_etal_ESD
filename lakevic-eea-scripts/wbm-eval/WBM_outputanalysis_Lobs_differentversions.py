"""


WBM after correcting persiann data for nan days, no longer closes, figure this out


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


#%%


#===========#
#== DATA  ==#
#===========#

workDIR = os.path.join(r"C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea" )

inDIR_Py_nan_out = os.path.join(workDIR,r'WBM-git\output\2022\output\obs\PERSIANN\v01r01_out') # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output'
inDIR_Py_nan_ac = os.path.join(workDIR,r'WBM-git\output\2022\output\obs\PERSIANN\v01r02_ac') # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output'
inDIR_Py_nan_1985 = os.path.join(workDIR,r'WBM-git\output\2023\PERSIANN\v01r02_nan1985') # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output'

inDIR_Py_fillnan_out = os.path.join(workDIR,r'WBM-git\output\2023\PERSIANN\v01r02_out')  # PERSIANN fillnan
inDIR_Py_fillnan_ac = os.path.join(workDIR,r'WBM-git\output\2023\PERSIANN\v01r02_ac')
inDIR_Py_fillnan_out_1985 = os.path.join(workDIR,r'WBM-git\output\2023\PERSIANN\v01r02_1985')

inDIR_Obs = os.path.join(workDIR,r'lakevic-eea-scripts\data\data-modified\lakelevels' ) # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/output/lakelevels/'

fig_path = os.path.join(workDIR,r'lakevic-eea-scripts\figures\figures_5_mar23_wbmeval') # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/figures/figures_9_may22_WBMoutput_py'


# Open python versions 
#=============================================================================

def open_LL(dirname):
    filename = 'WBM_run_PERSIANN_obs_v01*_2020_mm.csv'    
    filepath = glob.glob(os.path.join(dirname,filename))[0]  
    WBM_df = pd.read_csv(filepath, parse_dates=['date'], index_col=('date'))
    #WBM_df.columns = ['L_wb', 'P_lake', 'E_lake', 'Q_in', 'Q_out']
    LL_df = WBM_df.iloc[:,0].to_frame(name="L_wb")
    return LL_df

def open_LL_obs(dirname, filename):
    filepath = os.path.join(dirname, filename) 
    LL_df = pd.read_csv(filepath, parse_dates=['date'], index_col=('date')).rename(columns={ 'water_level': 'L_obs'})

    return LL_df
#%%
# open model output
LL_nan_out = open_LL(inDIR_Py_nan_out)
LL_nan_ac = open_LL(inDIR_Py_nan_ac)
LL_nan_1985 = open_LL(inDIR_Py_nan_1985)

LL_fillnan_out= open_LL(inDIR_Py_fillnan_out)
LL_fillnan_ac= open_LL(inDIR_Py_fillnan_ac)
LL_fillnan_out_1985 = open_LL(inDIR_Py_fillnan_out_1985)



# open observations
LL_obs = open_LL_obs(inDIR_Obs, 'lakelevel_all_intr_DH22.csv')
#%%

# plot all different options 
fig, ax = plt.subplots(figsize=(10,7)) 
LL_nan_out['L_wb'].plot(ax=ax, label = 'nan out')
LL_nan_ac['L_wb'].plot(ax=ax, label = 'nan AC')
LL_nan_1985['L_wb'].plot(ax=ax, label = 'nan out 1985-2020')

LL_fillnan_out['L_wb'].plot(ax=ax, label = 'fillnan out') # not ok
LL_fillnan_ac['L_wb'].plot(ax=ax, label = 'fillnan AC') # ok
LL_fillnan_out_1985['L_wb'].plot(ax=ax, label = 'fillnan out 1985-2020') # not ok

LL_obs['L_obs'].loc["1983-01-01":"2020-12-31"].plot(ax=ax, label = 'observations', c='k', ls='--')
plt.legend()

#plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_fillnanprecip_cfr_v2.png'),dpi=300)
#plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_fillnanprecip_cfr_v2.pdf'),dpi=300)


# %%
fig, ax = plt.subplots(figsize=(9,6)) 
LL_nan_out['L_wb'].plot(ax=ax, label = 'uncorrected precipitation, semi-observed outflow')
LL_nan_ac['L_wb'].plot(ax=ax, label = 'uncorrected precipitation, agreed curve outflow')
LL_fillnan_ac['L_wb'].plot(ax=ax, label = 'corrected precipitation, agreed curve outflow')

LL_obs['L_obs'].loc["1983-01-01":"2020-12-31"].plot(ax=ax, label = 'observations', c='k', ls='--')
plt.legend()

#plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_fillnanprecip_cfr_v2.png'),dpi=300)
#plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_fillnanprecip_cfr_v2.pdf'),dpi=300)

# include this as an SI ?? 

# %%
