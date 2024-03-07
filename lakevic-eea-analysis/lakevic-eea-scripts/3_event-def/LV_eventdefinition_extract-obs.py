#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

THESIS : LAKE VICTORIA EVENT ATTRIBUTION: EVENT DEFINITION

    EEA script 1. 
    
Extract event as dL/dt (for many different dts) for WBM-OBSERVATIONAL RUN 
As an evaluation of the WBM

Created on Wed May 18 17:58:30 2022

@author: rosapietroiusti

Update March 2023 (INCOMPLETE UPDATE)

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

inDIR_Obs = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/output/lakelevels/'
inDIR_Py = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output'

fig_path = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/figures/figures_10_may22_eventdef'
out_path = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/output/event_definition'


#LL_observed, interpolated 
filepath = os.path.join(inDIR_Obs, 'lakelevel_all_intr_DH22.csv') #  Lake levels linearly interpolated to daily resolution, merged HM and DH22 (downloaded April 2022)
lakelevels_obs_intr_raw = pd.read_csv(filepath)

#LL WBM_observational run PERSIANN 1983-2020
# old version, agreed curve from 2011: filepath = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output/v03_obs/WBM_run_v03r02_mm.csv'
filepath = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion_HPC/local/output/PERSIANN/v01r01_out/WBM_run_PERSIANN_obs_v01r01_1983_2020_mm.csv'
lakelevels_WBM_observational_raw = pd.read_csv(filepath)

#LL_WBM GCMs - hist
# 1) GFDL-ESM4
filepath = '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output/v03_hist/WBM_run_ISIMIP3b_hist_GFDL-ESM4_v03r01_1950_2020_mm.csv'
lakelevels_WBM_hist_raw = pd.read_csv(filepath)

#LL_WBM GCMs - hist-nat
#1) GFDL-ESM4

# Choose data 
#=============================================================================
#df = lakelevels_WBM_hist_raw #lakelevels_WBM_observational_raw 
#df = lakelevels_obs_intr_raw 
df = lakelevels_WBM_observational_raw 

# Clean data sets to datetime
#=============================================================================
df['date'] = pd.to_datetime(df['date'])
LL_raw = df.set_index(['date'])


#=============================================================================
# Settings
#=============================================================================
run_type = 'run_observational' # ''observed'#  #  'run_GCM_hist' # 'run_observational' #'run_observational'#  'run_ISIMIP_histnat_GCM'...
dt_type = 'list_long_final' # 'list_short' # 'list_long_final' #'list_longer'#'list_long2' 
model_name = 'GFDL-ESM4'

# Extract column of interest
if run_type == 'observed':
    LL_colname = 'water_level'
else:
    LL_colname = 'L_obs' # \\\\ change to L_wb !! \\\\
    
# Take only column of interest (lake levels)
LL_df = pd.DataFrame(LL_raw[LL_colname].copy())

# Get years 
startYEAR = LL_df.index[0].year #1949 # cut to full years?
endYEAR = LL_df.index[len(LL_df)-1].year # 2021

# Set list of dts
if dt_type == 'list_long':
    dt_list = list([30,60,90,120,135,180,270, 365,455, 545,635, 730]) # OLD LIST: list([90,135,180,270, 365,455, 545,635, 730]) # OLDER LIST: list([90,120,180,270, 365,455, 545,635, 730])
if dt_type == 'list_long2':
    dt_list = list([30,60,90,120,150,180,270, 365,455, 545,635, 730]) # OLD LIST: list([90,135,180,270, 365,455, 545,635, 730]) # OLDER LIST: list([90,120,180,270, 365,455, 545,635, 730])
if dt_type == 'list_longer': # only for observed and WBM observational run
    dt_list = list([30,60,90,120,150,180,210,240,270, 365,455, 545,635, 730]) # OLD LIST: list([90,135,180,270, 365,455, 545,635, 730]) # OLDER LIST: list([90,120,180,270, 365,455, 545,635, 730])
if dt_type == 'list_long_final': # for plotting
    dt_list = list([30,60,90,120,180,270, 365,455, 545,635, 730]) # the same as in the rank plot and other script, use this ! 
if dt_type == 'list_short':
    dt_list = list([60,90,120,150,180,365]) # FOR GCM-DRIVEN WBM
if dt_type == 'list_shortest':
    dt_list = list([60,180,365])
    
    # WHY DO I NEED 150??? Get rid of 150 !! 
    
    
# Settings for plots and saving
if run_type == 'run_observational':
    run_type_title = 'WBM observational run'
if run_type == 'observed':
    run_type_title = 'observed lake levels'
if run_type == 'run_GCM_hist' :
    run_type_title = 'WBM GCM hist run ({})'.format(model_name)
    run_type_save = '{}_{}'.format(run_type,model_name)
else: 
    run_type_save = run_type

# -Plot-
LL_df[LL_colname].plot()
plt.title('{}'.format(run_type_title))
plt.show()


#%% Plot lake levels if observed or observational
#====================
# PLOT LAKE LEVELS 
#====================

if (run_type == 'observed' or run_type == 'run_observational'):
    
    # Get summary stats
    meanL_all = LL_df[LL_colname].mean()
    maxL_all = LL_df[LL_colname].max()
    minL_all = LL_df[LL_colname].min()
    rangeL_all = maxL_all-minL_all
    p10_L = scipy.stats.scoreatpercentile(LL_df[[LL_colname]], 10) # 10th percentile 
    p90_L = scipy.stats.scoreatpercentile(LL_df[[LL_colname]], 90) # 90th percentile 

    
    
    # Get extreme high levs
    high_levs = LL_df[LL_df[LL_colname] > p90_L]
    list_years = np.unique(high_levs.index.year)
    
    high_dates = []
    high_levs_list =[]
    
    for i in list_years:
        #print(i)
        filter_condition = (high_levs.index.year == i)
        df_filter = high_levs.loc[filter_condition]
        maxL = max(df_filter[LL_colname])
        maxday = df_filter[df_filter[LL_colname] == maxL].index[0]
        high_levs_list.append(maxL)
        high_dates.append(maxday)
    
    # Select only extreme levels I want
    data = {'date': high_dates,
            'level': high_levs_list
            }
    df = pd.DataFrame(data)
    df = df.sort_values(by='level', ascending=False).reset_index(drop=True)
    if run_type == 'observed':
        select = [0,2,3,9] # THIS IS HARD CODED !! 
    else: 
        select = [0,1,8] # THIS IS HARD CODED !! 
    df_select = df.iloc[select]
    
    # Plot
    left = LL_df.index[0] # '{}-01-01'
    right = LL_df.index[len(LL_df)-1] # '{}-03-03'
    
    # Plot data
    fig, ax = plt.subplots() #figsize=(7,5)
    LL_df[LL_colname].plot(ax=ax, )
    
    # Plot settings
    ax.grid(True)
    plt.xlabel("date")
    plt.ylabel("lake level [m.a.s.l.]")
    if run_type == 'observed':
        plt.title("observed lake levels (WMO 1948-1992, DAHITI 1992-2022)")
    else:
        plt.title("WBM {} ({}-{})".format(run_type, left.year, right.year)) # add more info, better title
    ax.set_ylim(minL_all-0.2,maxL_all+0.7)
    plt.axhline(p90_L, ls = '--')
    plt.text(left+pd.Timedelta(days=(60)),p90_L+.05,'p90={:.1f} m'.format(p90_L), color='C0', weight='bold')
    #ax.set_xlim(left, right)
    
    # Add labels on extremes 
    for tick,label,i in zip(df_select['date'] , df_select['level'], range(4)):
        print(i)
        print(tick)
        print(label)
        if i == 0:
            if run_type == 'observed':
                shift = 5
            else: 
                shift = 3
            tick_plot = tick - pd.Timedelta(days=(365*shift))
        else: 
            tick_plot = tick - pd.Timedelta(days=(365*3.5))
        ax.text(tick_plot, #x
                label+0.1 ,  #y
                '{} \n {:.2f} m'.format(tick.strftime('%Y-%m-%d'), label), #text
                horizontalalignment='center', weight='bold', color='C0') # 
    
    
    fig.tight_layout()
    
    #plt.savefig(os.path.join(fig_path,'LL_plot_{}_{}_{}_v04r01.png'.format(left.year, right.year, run_type)),dpi=300)
    



#%% calc dLdt for dtlist
#====================
# Make loop for delta L between day of interest and dt days before that day 
# Plot these on top of each other or as suplots! and save!! 


for dt in dt_list:
    LL_df['dLdt_{}'.format(dt)] = np.zeros(len(LL_df))
    
    # loop over each day from dt to the end, calculate dL/dt, insert right
    for i in list(LL_df.index[dt:len(LL_df)]):  
        #print(i)
        #print(LL_df['dLdt_{}'.format(dt)].loc[i])
        deltaL =  LL_df[LL_colname].loc[i] - LL_df[LL_colname].loc[i - pd.Timedelta(days=dt)]
        LL_df['dLdt_{}'.format(dt)].loc[i] = deltaL


#%%
# Plot : subplots of dL/dt for different dts 
#========================================

dt_list_plot = [60,180,365] # short
#dt_list_plot = dt_list # all 
label_list = ['a)', 'b)', 'c)']

# Amount of plots
rows = math.ceil(len(dt_list_plot) /4)
if len(dt_list_plot) > 3:
    cols = 4 
    titlesize=20
    figsize=(15,10 )
else:
    cols = 3
    titlesize=12
    figsize=(6,6)

if len(dt_list_plot) < 4:
    fig, axes2d = plt.subplots(nrows=cols, ncols=rows, figsize=figsize, sharex=True ) #  sharex=True, , figsize=(20,20 )
else: 
    fig, axes2d = plt.subplots(nrows=rows, ncols=cols, figsize=figsize ) #
if (rows+cols) > 1:
    axes2d = axes2d.flatten()

for dt,ax, lab in zip(dt_list_plot, enumerate(axes2d), label_list):
    ax = ax[1]
    plot = LL_df['dLdt_{}'.format(dt)].plot(label='dLdt_{}'.format(dt), ax=ax)
    ax.set_title('$\Delta$L/$\Delta$t, $\Delta$t={} days'.format(dt))
    ax.set_ylabel("[m]")
    ax.text(-.15,1.1, lab, transform=ax.transAxes, fontsize=15)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        
plt.xlabel("")
#if run_type == 'observed':
    #fig.suptitle('event definition: \n {} ({}-{})'.format(run_type_title, startYEAR, endYEAR), fontsize=titlesize )
#else: 
    #fig.suptitle('{} ({}-{})'.format(run_type_title, startYEAR, endYEAR), fontsize=titlesize )
fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'dLdt_eventdef_{}_{}_{}_dt{}_small_v2.png'.format(run_type_save,startYEAR,endYEAR,dt_type)),dpi=300)



#%% Save dLdt
#====================
df_save = LL_df 
df_save = df_save.drop(columns=[LL_colname])

#df_save.to_csv( os.path.join(out_path, 'dLdt_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt_type, run_type_save)))



#%% Extract annual blockmax 
#====================
# Extract the maximum dL/dt for each year and for each dt window 

# Get years for which you can compute block max
list_years = np.unique(LL_df.index.year)
if run_type == 'observed': 
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
    #print(i)
    filter_condition = (LL_df.index.year == i)
    df_filter = LL_df.loc[filter_condition] # filter each year into a df
    if i == 2020 and (run_type == 'observed' or run_type == 'run_observational'):
        df_filter_2020 = df_filter.copy() 
        
    for dt in dt_list: 
        #print(dt)
        
        # max dL/dt for each year, put in df_blockmax
        max_dL = max(df_filter['dLdt_{}'.format(dt)]) # store what day this is on
        df_blockmax.loc[df_blockmax['year'] == i, 'dLdt_{}'.format(dt) ] = max_dL
        
        # lwhen it happened, ist of max days
        max_day = list(df_filter[df_filter['dLdt_{}'.format(dt)] == max_dL].index)
        # if dt == 180:
        #     print(dt)
        #     print(max_day)
        
        # just one max day to put in the df_daymax
        max_day = df_filter[df_filter['dLdt_{}'.format(dt)] == max_dL].index[0]
        df_daymax.loc[df_daymax['year'] == i, 'daymax_{}'.format(dt) ] = max_day #dLdt_{}

#%%
# join dataframes for blockmax and save
df_join_blockmax = df_blockmax.join(df_daymax.set_index('year'), on='year', how='left', rsuffix='_rx')
#df_join_blockmax.to_csv( os.path.join(out_path, 'blockmax_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt_type, run_type_save)), index=False)


#%% Sensitivity test of 2020 event : RANK
# for observed or observational run
#========================================

dt_list_plot = dt_list # [30,60,90,120,180,270,365,455,545,635,730]

if (run_type == 'observed' or run_type == 'run_observational'):
    
    # Get rank of each block max
    df_ranks = df_blockmax.copy()
    for dt in dt_list :       
        df_ranks = df_ranks.sort_values(by='dLdt_{}'.format(dt), ascending=False)
        df_ranks['dLdt_{}_rank'.format(dt)] = range(1, len(df_ranks)+1)

    # Extract for 2020
    d = np.zeros( (len(dt_list), 2) )
    df_2020 = pd.DataFrame(d, columns = ['dt', 'rank'])
    df_2020['dt'] = dt_list
    df_2020['rank'] = list(df_ranks[df_ranks['year'] ==2020].iloc[0, len(dt_list)+1:len(df_ranks.columns) ]) 
    
    df_2020 = df_2020[df_2020['dt'].isin(dt_list_plot)]
    
    # Plot dL/dt of 2020 as a function of dt 
    fig, ax = plt.subplots() # figsize=(5.5,4.5)figsize=(7,5) figsize=(5.5,5) figsize=(6,4.5)
    ax.grid(True,  linestyle='dashed')
    ax.set_axisbelow(True)
    ax.scatter(df_2020['dt'], df_2020['rank'], s=200, c=df_2020['rank'], cmap='coolwarm_r') # , alpha=0.5
    ax.set_ylim(-1,20)
    ax.set_xticks(df_2020['dt'])
    list_xticklabels = df_2020['dt']
    ax.set_xticklabels( list_xticklabels)
    plt.gca().invert_yaxis()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel("$\Delta$t (days)")
    plt.ylabel("rank of 2020 event")
    #plt.title('{}: \n rank of 2020 event ($\Delta$L/$\Delta$t) in {}-year ({}-{}) \n annual block maxima timeseries as a function of $\Delta$t '
    #          .format(run_type_title, (list_years_loop[len(list_years_loop)-1]-list_years_loop[0]+1), list_years_loop[0], list_years_loop[len(list_years_loop)-1]),  fontsize=11) # 
    plt.xticks(rotation = 45, ha="right",rotation_mode="anchor" )
    
    fig.tight_layout()
    #plt.savefig(os.path.join(fig_path,'dLdt_rank_dt{}_{}_v3.png'.format(dt_type,run_type)),dpi=300)


