# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 21:08:14 2023

@author: rpietroi

Look at isimip3a lake level simulations (and then dLdt representation in obsclim and counterclim)
Could also eval the precip data spatial + seasonal pattern 

"""



#=============#
#== IMPORTS ==#
#=============#

import numpy as np
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

inDIR_Obs = os.path.join(workDIR,r'lakevic-eea-scripts\data\data-modified\lakelevels' ) # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/output/lakelevels/'
fig_path = os.path.join(workDIR,r'lakevic-eea-scripts\figures\figures_16_jun23_isimip3a') # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/figures/figures_9_may22_WBMoutput_py'


# wbm-obsclim
filepath = glob.glob(os.path.join(workDIR,'WBM-git/output/2023/obsclim/20CRv3-ERA5/v01r01/*_mm.csv'))[0]
WBM_obs_raw = pd.read_csv(filepath, parse_dates=[0], index_col=0)

# wbm-counterclim
filepath = glob.glob(os.path.join(workDIR, 'WBM-git/output/2023/counterclim/20CRv3-ERA5/v01r01/*_mm.csv'))[0]
WBM_cou_raw = pd.read_csv(filepath, parse_dates=[0], index_col=0)

# observed lake levels 

filepath = os.path.join(inDIR_Obs, 'lakelevel_ext_intr_HCDH22.csv') 
lakelevels_obs_raw = pd.read_csv(filepath, parse_dates=[0], index_col=0)


# add also paths to dLdt when I calc this 

#%%

WBM_obs = WBM_obs_raw[['L_wb']]
WBM_cou = WBM_cou_raw[['L_wb']]
LL_obs = lakelevels_obs_raw

startYEAR = LL_obs.index[0].year
endYEAR = LL_obs.index[len(LL_obs)-1].year

# change to in-situ, by shifting to lake level in 1901 
#shift = WBM_obs.iloc[0].values - LL_obs.loc[WBM_obs.index[0]].values

#WBM_obs = WBM_obs - shift
#WBM_cou = WBM_cou - shift

shift = 1123.32
LL_obs = LL_obs + shift

#%% Plot 

fig, ax = plt.subplots(figsize=(7,4))
ax.plot(LL_obs['water_level'],label="observed", color='C2' ) 
ax.plot(WBM_obs['L_wb'], label="modelled obsclim" )
ax.plot(WBM_cou['L_wb'], label="modelled counterclim" )
ax.set_ylabel("lake level (m a.s.l.)")
ax.legend(frameon=False, loc='upper left')
ax.set_axisbelow(True)

#plt.savefig(os.path.join(fig_path,'LL_obsclim_counterclim_1901-2020.png'),dpi=300)
#plt.savefig(os.path.join(fig_path,'LL_obsclim_counterclim_1901-2020.pdf'),dpi=300)


#%%

WBM_obs_raw['P_lake'].resample('Y').sum().plot()
WBM_obs_raw['Q_in'].resample('Y').sum().plot()
WBM_obs_raw['Q_out'].resample('Y').sum().plot() # Qout very low !  initialise model with lake level in late 1900s!! 
#WBM_obs_raw['E_lake'].resample('Y').sum().plot() 
plt.legend()

#plt.savefig(os.path.join(fig_path,'WB_obsclim_1901-2020.png'),dpi=300)
#plt.savefig(os.path.join(fig_path,'WB_obsclim_1901-2020.pdf'),dpi=300)


# way off !!! do an eval of the ISIMIP3a precip data in the supercomputer and see if the w5e5 dataset is much better! 

# probably needs a bit of spin up time !! 