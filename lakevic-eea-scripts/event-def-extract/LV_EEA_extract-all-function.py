# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:20:38 2023

@author: rpietroi

Extract climate variable (dL/dt for dt=180 days) from timeseries of lake levels
Extract blockmax
Save, for analysis in climate explorer or R 

"""

#%%
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

inDIR_Obs = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\data\data-modified\lakelevels'
inDIR_PyObs = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\WBM-git\output\2022\output\obs' 
inDIR_PyGCM = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\WBM-git\output\2022\output'

out_path = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\output\event_def_dt180' 

#%% define functions

def extract_dLdt(data, LL_colname, dt_list):
    # extract data and clean
    if isinstance(data, str) == True:
        WBM_df = pd.read_csv(data,
                         index_col = 'date',
                         parse_dates=True)
    elif isinstance(data, pd.DataFrame):
        WBM_df = data
    else:
        print('error check datatype')
    LL_df = pd.DataFrame(WBM_df.iloc[:,0]) # Select only first column, this should have lake levels
    LL_df.columns = [LL_colname] # give the LL_colname as a string and here i turn it into list?

    # if dt is not a list make it a list 
    if isinstance(dt_list, list) == True:
        dt_list = dt_list
    else:
        dt_list = [dt_list]
        
    # extract dLdt 
    for dt in dt_list:
        LL_df['dLdt_{}'.format(dt)] = np.zeros(len(LL_df))
        # loop over each day from dt to the end, calculate dL/dt, insert right (date indicates last day)
        for i in list(LL_df.index[dt:len(LL_df)]):  
            deltaL =  LL_df[LL_colname].loc[i] - LL_df[LL_colname].loc[i - pd.Timedelta(days=dt)]
            LL_df['dLdt_{}'.format(dt)].loc[i] = deltaL.round(6) #deltaL[0]
    
    return LL_df

def save_dLdt(data, outpath, dt, run_save):
    startYEAR = data.index[0].year 
    endYEAR = data.index[len(LL_df)-1].year
    data.to_csv( os.path.join(outpath, 'dLdt_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt, run_save)), index=True)
    
def extract_dLdt_blockmax(data, dt):
    LL_df = data
    
    # remove first year, and if last year is incomplete remove that one too
    list_years = np.unique(LL_df.index.year)
    if len(LL_df.index.year[LL_df.index.year == list_years[-1]]) < 365:
        list_years_loop = list_years[1:len(list_years)-1] # get rid of first and last year
    else:
        list_years_loop = list_years[1:len(list_years)] # get rid only of first year
    
    # if dt is not a list make it a list 
    if isinstance(dt, list) == True:
        dt_list = dt
    else:
        dt_list = [dt]
        
    # get column names for daymax
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
    df_join_blockmax = df_blockmax.join(df_daymax.set_index('year'), on='year', how='left', rsuffix='_rx').set_index('year')
    
    return df_join_blockmax
        
def save_dLdt_blockmax(data, outpath, dt, run_save):
    startYEAR = data.index[0]
    endYEAR = data.index[len(data)-1]
    data.to_csv( os.path.join(outpath, 'blockmax_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt, run_save)), index=True)
    

def downscale_and_intrp(df):
    df_out = df[df.index.day == 1].resample('D').asfreq().interpolate(method='linear').round(3) 

    return df_out 

#%%

def remove_overlap_blocks(df_blockmax, dt):
    # remove all blocks that end between october and december

    list_drop = []
    for i in range(len(df_blockmax)):
        if df_blockmax['daymax_{}'.format(dt)].iloc[i].month >= 10:
            list_drop.append(df_blockmax.index[i])

    df_drop = df_blockmax.drop(index=list_drop)
    print('dropped {} rows'.format(len(list_drop)))

    return df_drop



#%% run on  extended levels 1800s-2020

"""
OBSERVATIONAL
Extract dLdt and annual blockmax on extended Lake Levels 
"""

# set params
filepath = os.path.join(inDIR_Obs, 'lakelevel_ext_intr_HCDH22.csv')
sim_type = run_save = dir_name = 'obs'
dt = 180

# extract dLdt
LL_df = extract_dLdt(filepath, 'LL', dt)

# save dLdt 
outpath = os.path.join(out_path,dir_name)
#save_dLdt(LL_df, outpath, dt, run_save)

# extract and save blockmax
df_blockmax = extract_dLdt_blockmax(LL_df, dt)
#save_dLdt_blockmax(df_blockmax, outpath, dt, run_save)

#%%


"""
OBSERVATIONAL SENSITIVITY TEST
Downscale/Upscale extended timeseries late 1800s-2020 to monthly resolution to see if this introduces a bias in the PR you calculate 
"""

filepath = os.path.join(inDIR_Obs, 'lakelevel_ext_intr_HCDH22.csv')
dt = 180
run_save = 'obs_downscaled'
dir_name = 'obs'

df = pd.read_csv(filepath,
                     index_col = 'date',
                     parse_dates=True)
# replace this with function above 
df_firstday_intr = df[df.index.day == 1].resample('D').asfreq().interpolate(method='linear').round(3) 

LL_df = extract_dLdt(df_firstday_intr, 'water_level', dt)

# save dLdt 
outpath = os.path.join(out_path,dir_name)
save_dLdt(LL_df, outpath, dt, run_save)

# extract and save blockmax
df_blockmax = extract_dLdt_blockmax(LL_df, dt)
#save_dLdt_blockmax(df_blockmax, outpath, dt, run_save)

#%% 

"""
OBSERVATIONAL SENSITIVITY TEST
Remove blockmax of years that overlap (when windows end between October and December) to rerun attribution code and
check sensitivity to this 
"""

# set params
filepath = os.path.join(inDIR_Obs, 'lakelevel_ext_intr_HCDH22.csv')
run_save = 'obs_rm_overlap'
dir_name = 'obs'
dt = 180

# extract dLdt
LL_df = extract_dLdt(filepath, 'LL', dt)

# save dLdt 
outpath = os.path.join(out_path,dir_name)

# extract and save blockmax
df_blockmax = extract_dLdt_blockmax(LL_df, dt)

# remove overlapping blocks
df_blockmax_drop = remove_overlap_blocks(df_blockmax, dt)
#save_dLdt_blockmax(df_blockmax_drop, outpath, dt, run_save)



# %%
