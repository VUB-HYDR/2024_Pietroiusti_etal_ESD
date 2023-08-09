#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:15:28 2022

@author: rosapietroiusti

THESIS : LAKE VICTORIA EVENT ATTRIBUTION: ATTRIBUTION

    Script 2. 

Extract climate variable (dL/dt for dt=180 days) from timeseries of lake levels
Extract blockmax
Save, for analysis in climate explorer or R 

"""

#=============#
#== IMPORTS ==#
#=============#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os , glob
import math
import pandas as pd
from matplotlib import transforms
from matplotlib.ticker import MaxNLocator
import scipy as scipy
from scipy import stats
from matplotlib.pyplot import cm

#===========#
#== DATA  ==#
#===========#

inDIR_Obs = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\data\data-modified\lakelevels' #'/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/output/lakelevels/'
inDIR_PyObs = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\WBM-git\output\2022\output\obs' # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion_HPC/local/output'
inDIR_PyGCM = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\WBM-git\output\2022\output' # '/Users/rosapietroiusti/Documents/Terminal_HPC/Thesis/WBM/output'

#fig_path = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/figures/figures_11_jun22_obsprob'
out_path = 'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\output\event_def_dt180' #'/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/output/event_def_dt180'

#%%

# Select flags 
run_type = 2    
                #0: observed lake levels 
                #1: GCM 
                #2: WBM with observational precip
               
                
if run_type == 1:
    run_sim = 1 
                #0: hist
                #1: hist-nat

if run_type ==1:
    run_model = 5
                # 0: CanESM5
                # 1: CNRM-CM6-1
                # 2: GFDL-ESM4
                # 3: IPSL-CM6A-LR
                # 4: MIROC6
                # 5: MRI-ESM2-0

if run_type == 2:
    prec_source = 0
                # 0: PERSIANN
    out_type = 1
                # 0: semi-observational outflow
                # 1: agreed curve outflow 

#%%

# Get to the correct directory

ver = 'v01r01'
# observational lake levels 
if run_type == 0:
    filepath = os.path.join(inDIR_Obs, 'lakelevel_all_intr_DH22.csv' ) #LL_observed, interpolated 
    sim_type = run_save = dir_name = 'obs'
# GCM, ISIMIP3b
elif run_type == 1: 
    GCM_list = ['CanESM5','CNRM-CM6-1','GFDL-ESM4','IPSL-CM6A-LR','MIROC6','MRI-ESM2-0']
    if run_sim == 0: 
        sim_type = 'hist'
        dir_name =  'hist-rcp70'
        if run_model == 2:
            inDIR = os.path.join(inDIR_PyGCM,dir_name , GCM_list[run_model], ver)
        else:
            inDIR = os.path.join(inDIR_PyGCM, dir_name, GCM_list[run_model], ver)
    if run_sim == 1: 
        sim_type = dir_name = 'hist-nat'
        inDIR = os.path.join(inDIR_PyGCM, sim_type , GCM_list[run_model], ver)
    filepath = os.path.join(inDIR, 'WBM_run_ISIMIP3b_{}_{}_{}_1850_2020_mm.csv'.format(sim_type, GCM_list[run_model], ver ))
    run_save = str(sim_type + '_' +GCM_list[run_model])
# WBM observational run
elif run_type == 2: 
    if prec_source == 0: 
        prec_source = 'PERSIANN'
        if out_type == 0:
            ver_name = ver
            indir = '{}_out'.format(ver_name)
            run_save = '{}_obs_out'.format(prec_source)
        elif out_type == 1:
            ver_name = 'v01r02'
            indir = '{}_ac'.format(ver_name)
            run_save = '{}_obs_ac'.format(prec_source)
    dir_name = 'obs-WBM/{}'.format(prec_source)        
    filepath = os.path.join(inDIR_PyObs, prec_source, indir, 'WBM_run_{}_obs_{}_1983_2020_mm.csv'.format(prec_source, ver_name))
    


#%%
# Read and clean data 

WBM_df = pd.read_csv(filepath,
                     index_col = 'date',
                     parse_dates=True)

LL_df = pd.DataFrame(WBM_df.iloc[:,0]) # Select only first column
LL_df.columns = LL_colname = ['LL']


#%%

# # Get years 
startYEAR = LL_df.index[0].year #1949 # cut to full years?
endYEAR = LL_df.index[len(LL_df)-1].year # 2021

#%%
# Get dLdt for each day

# set dt
dt_list = [180]

for dt in dt_list:
    LL_df['dLdt_{}'.format(dt)] = np.zeros(len(LL_df))


    # loop over each day from dt to the end, calculate dL/dt, insert right
    for i in list(LL_df.index[dt:len(LL_df)]):  
        deltaL =  LL_df[LL_colname].loc[i] - LL_df[LL_colname].loc[i - pd.Timedelta(days=dt)]
        LL_df['dLdt_{}'.format(dt)].loc[i] = deltaL[0].round(6)

#%% Save dLdt
#====================
df_save = LL_df 
#df_save = df_save.drop(columns=LL_colname) # keep the lake level column, could be useful

df_save.to_csv( os.path.join(out_path,dir_name, 'dLdt_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt, run_save)), index=True)



#%% Extract annual blockmax 
#====================
# Extract the maximum dL/dt for each year and for each dt window 

# Get years for which you can compute block max
list_years = np.unique(LL_df.index.year)

if run_type == 0: 
    list_years_loop = list_years[1:len(list_years)-1] # get rid of first and last
else:
    list_years_loop = list_years[1:len(list_years)] # get rid only of first 

# Get column names for daymax 
list_colnames = []
for dt in dt_list:
    name = 'daymax_{}'.format(dt)
    list_colnames.append(name)


# initialise df to put the max dL/dt for each year for each dt
d = np.zeros( (len(list_years_loop), len(LL_df.columns)) )
df_blockmax = pd.DataFrame(d, columns = ['year'] + list(LL_df.columns)[1:len(LL_df.columns)])
df_blockmax['year'] = list_years_loop

# initialise df to put the day on which that max dL/dt happens
df_daymax = pd.DataFrame(d, columns = ['year'] + list_colnames) #list(LL_df.columns)[1:len(LL_df.columns)]
df_daymax['year'] = list_years_loop

# Loop over full years
for i in list_years_loop:
    
    filter_condition = (LL_df.index.year == i)
    df_filter = LL_df.loc[filter_condition] # filter each year into a df
    if i == 2020 and (run_type == 'observed' or run_type == 'run_observational'):
        df_filter_2020 = df_filter.copy() 
        
    for dt in dt_list: 
        
        # max dL/dt for each year, put in df_blockmax
        max_dL = max(df_filter['dLdt_{}'.format(dt)]) # store what day this is on
        df_blockmax.loc[df_blockmax['year'] == i, 'dLdt_{}'.format(dt) ] = max_dL
        
        # lwhen it happened, ist of max days
        max_day = list(df_filter[df_filter['dLdt_{}'.format(dt)] == max_dL].index)

        # just one max day to put in the df_daymax
        max_day = df_filter[df_filter['dLdt_{}'.format(dt)] == max_dL].index[0]
        df_daymax.loc[df_daymax['year'] == i, 'daymax_{}'.format(dt) ] = max_day #dLdt_{}

# join dataframes for blockmax and save
df_join_blockmax = df_blockmax.join(df_daymax.set_index('year'), on='year', how='left', rsuffix='_rx')

#%% 

# SAVE block max csv
df_join_blockmax.to_csv( os.path.join(out_path, dir_name, 'blockmax_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt, run_save)), index=False)

#%%

# Save txt file for KX
df_join_blockmax.drop(columns='daymax_{}'.format(dt)).to_csv( os.path.join(out_path, dir_name, 'blockmax_{}_{}_dt{}_{}.txt'.format(startYEAR, endYEAR, dt, run_save)), index=False, sep='\t')


