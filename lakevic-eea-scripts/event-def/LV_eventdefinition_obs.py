#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LAKE VICTORIA EVENT ATTRIBUTION: 
    
Event definition script from observational lake levels 
Extract from obs the dL/dt and make plots for event definition 
Figures 8-9 Pietroiusti et al 2023 ESD
Created on Wed May 18 17:58:30 2022
Update March 2023
@author: rosapietroiusti

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
import csv
import math
import pandas as pd
from matplotlib import transforms
from matplotlib.ticker import MaxNLocator
import scipy as scipy
from scipy import stats
from matplotlib.pyplot import cm
import matplotlib.patches as mpatches

#%%

#===========#
#== DATA  ==#
#===========#

inDIR_Obs = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\data\data-modified\lakelevels'
fig_path = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\figures\figures_4_mar23_eventdef' 
out_path = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\output\event-def'

# Observational lakelevels 
#=============================================================================
filepath = os.path.join(inDIR_Obs, 'lakelevel_all_intr_DH22.csv') #  Lake levels linearly interpolated to daily resolution, merged HM and DH22 (downloaded April 2022)
lakelevels_obs_intr_raw = pd.read_csv(filepath)

filepath = os.path.join(inDIR_Obs, 'lakelevel_all_raw_DH22.csv') # raw 1949-202
lakelevels_obs_raw = pd.read_csv(filepath)

filepath = os.path.join(inDIR_Obs, 'lakelevel_ext_intr_HCDH22.csv') #  Lake levels linearly interpolated to daily resolution 1897-2020
lakelevels_obs_intr_ext_raw = pd.read_csv(filepath)


# Clean data and set datetime
#=============================================================================
# LL_observed, interpolated 
df = lakelevels_obs_intr_raw
df['date'] = pd.to_datetime(df['date'])
LL_obs_raw = df.set_index(['date'])
# -Plot-
LL_obs = LL_obs_raw.copy()
LL_obs['water_level'].plot()
plt.show()

# LL_observed, gaps 
df = lakelevels_obs_raw
df['date'] = pd.to_datetime(df['date'])
LL_obs_gaps_raw = df.set_index(['date'])
# -Plot- 
LL_obs_gaps = LL_obs_gaps_raw.copy()
LL_obs_gaps['water_level'].plot()
plt.show()

#LL_observed, ext 
df = lakelevels_obs_intr_ext_raw
df['date'] = pd.to_datetime(df['date'])
LL_obs_raw = df.set_index(['date'])
# -Plot-
LL_obs = LL_obs_raw.copy()
LL_obs['water_level'].plot()
plt.show()



# Settings
#=============================================================================
startYEAR = 1897
endYEAR = 2021
run_type = 'observed'
dt_type = 'list'


#%% calc dLdt
# Make loop for delta L between day of interest and dt days before that day i.e. it is centered RIGHT
#=============================================================================

# List of dts
dt_list = list([60,180,365]) 
#dt_list = list([30,60,90,120,180,270, 365,455, 545,635, 730]) # shorter : list([90,135,180,270, 365,455, 545,635, 730]) # list([90,120,180,270, 365,455, 545,635, 730])
#dt_list = list([30,60,90,120,180,270, 365,455, 545,635, 730]) # for Figure 8a (rank of event as a function of dt)

for dt in dt_list:
    LL_obs['dLdt_{}'.format(dt)] = np.zeros(len(LL_obs))
    
    for i in list(LL_obs.index[dt:len(LL_obs)]):  
        deltaL =  LL_obs['water_level'].loc[i] - LL_obs['water_level'].loc[i - pd.Timedelta(days=dt)]
        LL_obs['dLdt_{}'.format(dt)].loc[i] = deltaL

    # plot to check     
    LL_obs['dLdt_{}'.format(dt)].plot(label='dLdt_{}'.format(dt))
    plt.title('dL/dt, dt={} days'.format(dt))
    plt.show()



#%%
# Plot : subplots of dL/dt for different dts 
#========================================

from datetime import date
from datetime import timedelta

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
    #plot = LL_obs['dLdt_{}'.format(dt)].plot(label='dLdt_{}'.format(dt), ax=ax)
    ax.plot(LL_obs['dLdt_{}'.format(dt)])
    ax.set_title('$\Delta$L/$\Delta$t, $\Delta$t={} days'.format(dt))
    ax.set_ylabel("(m)")
    ax.text(-.15,1.1, lab, transform=ax.transAxes, fontsize=15)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    
    if lab == 'b)':
        # Add red box
        rect=mpatches.Rectangle((date(2018,3,1) , .9),timedelta(365*3.5),.5, 
                                fill=False,
                                color="red",
                               linewidth=2)
                               #facecolor="red")
        ax.add_patch(rect)
        
plt.xlabel("")
plt.xlim(LL_obs.index[0],LL_obs.index[-1])
fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'dLdt_eventdef_{}_{}_dt{}_small.png'.format(startYEAR,endYEAR,dt_type)),dpi=300)
#plt.savefig(os.path.join(fig_path,'dLdt_eventdef_{}_{}_dt{}_small.pdf'.format(startYEAR,endYEAR,dt_type)),dpi=300)


#%% Save dLdt
#=============================================================================
df_save = LL_obs 
df_save = df_save.drop(columns=['water_level'])

#df_save.to_csv( os.path.join(out_path, 'dLdt_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt_type, run_type)))



#%% Extract blockmax (yearly)
#=============================================================================

# Extract the maximum dL/dt for each year and for each dt window 
list_years = np.unique(LL_obs.index.year)

# df to put the max dL/dt for each year for each dt
d = np.zeros( (len(list_years)-2, len(LL_obs.columns)) )
df_blockmax = pd.DataFrame(d, columns = ['year'] + list(LL_obs.columns)[1:len(LL_obs.columns)])
df_blockmax['year'] = list_years[1:len(list_years)-1]

# List column names 
list_colnames = []
for dt in dt_list:
    name = 'daymax_{}'.format(dt)
    list_colnames.append(name)

    
# df to put the day on which that max dL/dt happens
df_daymax = pd.DataFrame(d, columns = ['year'] + list_colnames) #list(LL_obs.columns)[1:len(LL_obs.columns)]
df_daymax['year'] = list_years[1:len(list_years)-1]


# Loop from 1949 to 2021 (full years), doesn't include 2022 
for i in list_years[1:len(list_years)-1]:
    #print(i)
    filter_condition = (LL_obs.index.year == i)
    df_filter = LL_obs.loc[filter_condition] # filter each year into a df
    
    for dt in dt_list: 
        
        # max dL/dt for each year, put in df_blockmax
        max_dL = max(df_filter['dLdt_{}'.format(dt)]) # store what day this is on
        df_blockmax.loc[df_blockmax['year'] == i, 'dLdt_{}'.format(dt) ] = max_dL
        
        # when it happened, list of max days
        max_day = list(df_filter[df_filter['dLdt_{}'.format(dt)] == max_dL].index)
        if dt == 135:
            if len(max_day) > 1: # where you have record on more than one day
                print(dt)
                print(max_day)
        
        # just one max day to put in the df_daymax
        max_day = df_filter[df_filter['dLdt_{}'.format(dt)] == max_dL].index[0]
        df_daymax.loc[df_daymax['year'] == i, 'daymax_{}'.format(dt) ] = max_day #dLdt_{}


# join dataframes for blockmax and save
#=============================================================================
df_join = df_blockmax.join(df_daymax.set_index('year'), on='year', how='left', rsuffix='_rx')
#df_join.to_csv( os.path.join(out_path, 'blockmax_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt_type, run_type)), index=False)

# Clean up datetime
df_join['year'] = pd.to_datetime(df_join['year'], format="%Y")
df_join = df_join.set_index(['year'])



#%% Get rank of each max event 
#=============================================================================
df_ranks = df_blockmax.copy() # for next plot
#df_ranks = df_join.copy() # for saving

for dt in dt_list : #list([90,180,365,450, 550, 730])
    
    df_ranks = df_ranks.sort_values(by='dLdt_{}'.format(dt), ascending=False)
    df_ranks['dLdt_{}_rank'.format(dt)] = range(1, len(df_ranks)+1) #[1:len(df)+1]

#df_ranks.to_csv( os.path.join(out_path, 'blockmax_ranks_{}_{}_dt{}_{}.csv'.format(startYEAR, endYEAR, dt_type, run_type)), index=False)


#%% Figure 8 paper: 
# (a) Sensitivity test of 2020 event definition by rank and 
# (b) annual block maxima timeseries for dt=180
#=============================================================================


# subplot 1: Extract rank of 2020 event
d = np.zeros( (len(dt_list), 2) )
df_2020 = pd.DataFrame(d, columns = ['dt', 'rank'])
df_2020['dt'] = dt_list
df_2020['rank'] = list(df_ranks[df_ranks['year'] ==2020].iloc[0, len(dt_list)+1:len(df_ranks.columns) ]) 

df_2020 = df_2020[df_2020['dt'].isin(dt_list)]

# subplot 2: plot timeseries of annual block max
# select dt=180
dt_select = 180
df_select = df_join[['dLdt_{}'.format(dt_select)]].copy()
# clean up
df_select['year'] = pd.to_datetime(df_select.index, format="%Y")
df_select = df_select.set_index(['year'])
#calculate rolling mean for 10-years
df_select['rolling'] = df_select['dLdt_{}'.format(dt_select)].rolling(10, center=True).mean()


# n rs cs
rs=2
cs=1

# initiate figure
fig = plt.figure(figsize=(6,7)) 

# 1) Plot dL/dt rank of 2020 as a function of dt 
#=============================================================================

ax = plt.subplot(rs, cs, 1)

ax.grid(True,  linestyle='dashed')
ax.set_axisbelow(True)
ax.scatter(df_2020['dt'], df_2020['rank'], s=200, c=df_2020['rank'], cmap='coolwarm_r') # , alpha=0.5
ax.set_ylim(0,27) #20 #-1,24
ax.set_xticks(df_2020['dt'])
list_xticklabels = df_2020['dt']
ax.set_xticklabels( list_xticklabels)
plt.gca().invert_yaxis()
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel("$\Delta$t (days)")
#plt.ylabel("rank of 2020 event")
#plt.ylabel("$rank_{2020}$")
plt.ylabel("rank $\Delta$L/$\Delta$t$_{max,2020}$")
plt.title('(a)', loc='left')
#plt.title('{}: \n rank of 2020 event ($\Delta$L/$\Delta$t) in {}-year ({}-{}) \n annual block maxima timeseries as a function of $\Delta$t '
#          .format(run_type_title, (list_years_loop[len(list_years_loop)-1]-list_years_loop[0]+1), list_years_loop[0], list_years_loop[len(list_years_loop)-1]),  fontsize=11) # 
plt.xticks(rotation = 45, ha="right",rotation_mode="anchor" )

# Add red box
rect=mpatches.Rectangle((155,1),50,25, 
                        fill=False,
                        color="red",
                       linewidth=2)
                       #facecolor="red")
plt.gca().add_patch(rect)

# 2) Plot annual blockmax timeseries for dt = 180
#=============================================================================

ax = plt.subplot(rs, cs, 2)

ax.step(df_select.index,df_select['dLdt_{}'.format(dt_select)], color='r', label='$\Delta$t = 180 days')
ax.plot(df_select['rolling'], color='lime', label='10-year rolling mean') 
ax.axhline(df_select.loc["2020", 'dLdt_{}'.format(dt_select)].values, color='gray', ls='--', label='$\Delta$L/$\Delta$t$_{max,2020}$') # label='2020 magnitude')
plt.legend()
plt.ylabel("$\Delta$L/$\Delta$t$_{max}$ (m)")
plt.title('(b)', loc='left')

fig.tight_layout()
#plt.savefig(os.path.join(fig_path,'dLdt_rank_dtlist_blockmax_tseries_{}_{}.png'.format(startYEAR, endYEAR)),dpi=300)
#plt.savefig(os.path.join(fig_path,'dLdt_rank_dtlist_blockmax_tseries_{}_{}.pdf'.format(startYEAR, endYEAR)),dpi=300)



#%%

# 2) Plot annual blockmax timeseries for dt = dtlist
#=============================================================================

#dt_list_plot = [60,180,365] # short
dt_list_plot = dt_list # all 

# Amount of plots
rows = math.ceil(len(dt_list_plot) /4)
if len(dt_list_plot) > 3:
    cols = 4 
    titlesize=20
    figsize=(15,10 )
else:
    cols = 3
    titlesize=12
    figsize=(4,6)

if len(dt_list_plot) < 4:
    fig, axes2d = plt.subplots(nrows=cols, ncols=rows, figsize=figsize, sharex=True ) #  sharex=True, , figsize=(20,20 )
else: 
    fig, axes2d = plt.subplots(nrows=rows, ncols=cols, figsize=figsize ) #
if (rows+cols) > 1:
    axes2d = axes2d.flatten()

for dt,ax in zip(dt_list_plot, enumerate(axes2d)):
    ax = ax[1]
    # Plot
    # plot = df_join['dLdt_{}'.format(dt)].plot(label='dLdt_{}'.format(dt), ax=ax)
    ax.step(df_join.index, df_join['dLdt_{}'.format(dt)], color='r')
    ax.plot(df_join['dLdt_{}'.format(dt)].rolling(10, center=True).mean(), color='lime')
    ax.set_title('$\Delta$L/$\Delta$t (m), $\Delta$t={} days'.format(dt))
    #ax.axhline(y=0, color='gray', ls='--')
    val=float(df_join.loc["2020", 'dLdt_{}'.format(dt)].values)
    ax.axhline(val, color='gray', ls='--', label='$\Delta$L/$\Delta$t$_{max,2020}$='+'{:.3f} m'.format(val)) 
    ax.legend()
    # Make years fit
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        
fig.suptitle('annual block max: {}$-${}'.format(startYEAR, endYEAR), fontsize=titlesize )


fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'annualblockmax_obs_{}_{}_dt{}.pdf'.format(startYEAR,endYEAR,dt_type)),dpi=300)
#plt.savefig(os.path.join(fig_path,'annualblockmax_obs_{}_{}_dt{}.png'.format(startYEAR,endYEAR,dt_type)),dpi=300)




#%% Plot 180 days of each year for sentisitivity analysis (rectangles): Figure in paper
#=============================================================================

# select dt=180
dt_select = 180
df_select = df_join[['dLdt_{}'.format(dt_select), 'daymax_{}'.format(dt_select)]].copy()
df_select['daymax_{}'.format(dt_select)] = pd.to_datetime(df_select['daymax_{}'.format(dt_select)])

df_select['day_end'] = df_select['daymax_{}'.format(dt_select)].dt.dayofyear
df_select['day_start'] = df_select['day_end'] - dt_select

df_select = df_select.reset_index()

# make x axis labels with name of month
startdate = np.datetime64('2000-01-01')
enddate = np.datetime64('2000-12-31')
alldays = pd.date_range(start=startdate,end=enddate)

list_ticks = np.zeros(12)
for i in range(1,13):
    sel = alldays[alldays.month == i].dayofyear[0] -1
    list_ticks[i-1] = sel

list_ticks_pre = list_ticks - 365
list_monthnames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

list_ticks_label = np.append(list_ticks_pre[7:12],list_ticks)
list_monthnames_labels = np.append(list_monthnames[7:12],list_monthnames)

# Plot 
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# make colormap
#length_cmap = max(df_select.iloc[:,3]) - min(df_select.iloc[:,3])
length_cmap = 365
colors = cm.inferno(np.linspace(0, 1, length_cmap))
colors = cm.inferno(np.linspace(0.3, 0.9, length_cmap))

j = 10
length_cmap = 365//j
colors = cm.inferno(np.linspace(0.3, 0.9, length_cmap))

#linewidth
lw=2 # 3

# plot
fig, ax = plt.subplots(figsize=(5.5,6))


j=0
k=0
for i in range(len(df_select)):

    #c = colors[df_select.iloc[i,3] // j -1]
    
    c='C0'
    if k<=0:
        ax.plot([df_select.iloc[i,3], df_select.iloc[i,4]] , #end, start
             [df_select.iloc[i,0], df_select.iloc[i,0]], # year
             c=c,
             lw=lw, label='180-day time window')
    else:
        ax.plot([df_select.iloc[i,3], df_select.iloc[i,4]] , #end, start
             [df_select.iloc[i,0], df_select.iloc[i,0]], # year
             c=c,
             lw=lw)
    k=+1
    
    if i < len(df_select)-1:
        if df_select.iloc[i,4] + 365 <= df_select.iloc[i-1,3] : 
        #if start day of subsequent year is before end day of previous year
        # or if end day of previous year is after start day of subsequent year=
            c='purple'
            if j<=0:
                ax.plot([df_select.iloc[i-1,3]  - 365, df_select.iloc[i,4]] , #end, start
                     [df_select.iloc[i,0], df_select.iloc[i,0]], # year
                     c=c,
                     lw=lw,
                     label='overlap with previous \nor subsequent year')
            else:
                ax.plot([df_select.iloc[i-1,3]  - 365, df_select.iloc[i,4]] , #end, start
                     [df_select.iloc[i,0], df_select.iloc[i,0]], # year
                     c=c,
                     lw=lw)
            j=+1
            print(df_select.iloc[i,2])
                  
        elif (df_select.iloc[i,3] >= df_select.iloc[i+1,4] + 365): 
            c='purple'
            ax.plot([df_select.iloc[i,3], df_select.iloc[i+1,4] + 365] , #end, start
                     [df_select.iloc[i,0], df_select.iloc[i,0]], # year
                     c=c,
                     lw=lw,
                     ls='-')
            
        else:
            c='C0'
    else:
        c='C0'
    
    
    #plt.axvline(x=0, ls='--', c='lightgray', lw=1)
    plt.xlim(-180,365)
    #plt.xlabel("\n annual block max $\Delta$L/$\Delta$t 180-day interval")
    plt.xticks(list_ticks_label, list_monthnames_labels, rotation = 45)
# plt.text(0.08, -0.18, 'previous year', transform=ax.transAxes, c='gray')
# plt.text(0.55, -0.18, 'block max year', transform=ax.transAxes, c='gray')

# if specifying portrait view
plt.text(0.08, -0.1, 'year n-1', transform=ax.transAxes, c='gray')
plt.text(0.55, -0.1, 'year n', transform=ax.transAxes, c='gray')
ax.axvspan(-180, 0, alpha=0.2, color='gray') #, hatch='//'

#plt.ylabel('block max year n', c='gray')
plt.ylabel('$\longleftarrow$', c='gray', size=18)
plt.text(-.175, .6, 'year n', transform=ax.transAxes, c='gray', rotation=90,horizontalalignment='center',
     verticalalignment='center', size=12)
#plt.arrow(.1,.5,0,.1, transform=ax.transAxes)

ax.yaxis.set_minor_locator(MultipleLocator(365))
plt.ylim(df_select.iloc[0,0]-np.timedelta64(1,'Y'),df_select.iloc[len(df_select)-1,0]+np.timedelta64(1,'Y'))
plt.legend(loc=(.73,.82), fontsize=9, framealpha=1) # .73,.88
plt.gca().invert_yaxis()

fig.tight_layout()
#plt.savefig(os.path.join(fig_path,'blockmax_intervals_{}_{}_dt{}_v3.png'.format(startYEAR,endYEAR,dt_select)),dpi=300)
#plt.savefig(os.path.join(fig_path,'blockmax_intervals_{}_{}_dt{}_v3.pdf'.format(startYEAR,endYEAR,dt_select)),dpi=300)


