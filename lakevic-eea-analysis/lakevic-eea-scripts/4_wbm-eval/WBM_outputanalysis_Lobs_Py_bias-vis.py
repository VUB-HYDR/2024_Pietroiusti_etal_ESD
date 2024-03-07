#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:35:48 2022

@author: rosapietroiusti

Pietroiusti et al 2023 ESD

WBM evaluation 
Plotscript: Python WBM output - look at BIAS cfr to OBS


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

workDIR = os.path.join(r"C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-clean" ) #"C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea"
inDIR_Py = os.path.join(workDIR,r'lakevic-eea-wbm\output\2022\output\obs\PERSIANN\v01r01_out')  #WBM-git\output\2022\output\obs\PERSIANN\v01r01_out
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

#%%

# Settings  
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


# settings, flags

col_obs = 'C0'
col_mod = 'C1'

flag_save=True


#%%
"""
PART 1: Paper figures 8 and C10: observed vs. modelled lake levels
""" 


#%% 1) Lake level plot observed vs. modelled : Fig 8 paper 
#=============================================================================


left = WBM_obs.index[0]
right = WBM_obs.index[len(WBM_obs)-1]
right= pd.to_datetime('2021-01-01') # to get nicer labels 

fig, ax = plt.subplots(figsize=(7,4)) 

# plot observed levels 
ax.plot(LL_obs['water_level'],label="observed", color=col_obs) #(1983-2020)
# plot modelled 
ax.plot(WBM_obs['L_obs'], label="modelled ", c=col_mod ) #(observational run) 
ax.set_ylabel("lake level (m a.s.l.)")
ax.set_ylim(1133,1137)
ax.grid(True)
ax.legend(frameon=False, loc='upper left')
ax.set_axisbelow(True)

#xaxis labels - 3 options ! 
xlist =  pd.to_datetime(list(range(1983,2022,2)), format='%Y')

# Set xaxis labels - options ! 
ax.xaxis.set_ticks(xlist)
ax.xaxis.set_ticklabels(xlist.strftime("%Y"), rotation = 45, ha="right", rotation_mode="anchor")

#zoom in to 2018-2020
#from matplotlib.dates import DateFormatter
#left= pd.to_datetime('2018-01-01')
#date_form = DateFormatter("%b-%Y")
#ax.xaxis.set_major_formatter(date_form)
#plt.xticks(rotation=30, ha='right', rotation_mode="anchor")

# Add box around legend manually
import matplotlib.patches as mpatches

rect=mpatches.Rectangle(( WBM_obs.index[0] + np.timedelta64(5,'M') ,1136.1),np.timedelta64(10,'Y'),.8, 
                        fill=True,
                       alpha=0.8,
                       facecolor="w", 
                       edgecolor='lightgray')
plt.gca().add_patch(rect)

# Plot residual
c= '0.6' 
c1= 'gray'
ax2=ax.twinx()
residual = WBM_obs.loc[str(left):str(right)]['L_obs'] - LL_obs.loc[str(left):str(right)]['water_level']
df_res = pd.DataFrame(residual)
df_res.columns = ['residual']
ax2.axhline(color='0.6', ls='--')
ax2.plot(df_res, color=c, label='bias')
ax2.set_ylim(-1,7)
ax2.set_ylabel(" lake level bias (m)", color=c1)
ax2.tick_params(axis='y', color=c1, labelcolor=c1)
ax2.yaxis.set_label_coords(1.07,0.13)
ax2.legend(loc=(0.015,0.78), frameon=False)

ax2.legend(loc=(0.015,0.77), frameon=False)

#plt.legend()
plt.xlim(left, right)

fig.tight_layout()

if flag_save == True:
    plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_newcolors.png'),dpi=300)
    plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_newcolors.pdf'),dpi=300)







#%% version for EGU23 - no grid and no bias line 
#=============================================================================

# Lake level plots and difference model - observed
left = WBM_obs.index[0]
right = WBM_obs.index[len(WBM_obs)-1]
right= pd.to_datetime('2021-01-01') # to get nicer labels 

fig, ax = plt.subplots(figsize=(7,4)) #figsize=(7,5)

ax.plot(LL_obs['water_level'],label="observed", color=col_obs) #(1983-2020) #C2
ax.plot(WBM_obs['L_obs'], label="modelled ", color=col_mod ) #(observational run) #C0

ax.set_ylabel("lake level (m a.s.l.)")
ax.grid(False)
ax.legend(loc='upper left')
ax.set_axisbelow(True)

#xaxis labels
xlist =  pd.to_datetime(list(range(1983,2022,2)), format='%Y')

# Set xaxis labels  
ax.xaxis.set_ticks(xlist)
ax.xaxis.set_ticklabels(xlist.strftime("%Y"), rotation = 45, ha="right", rotation_mode="anchor")

plt.xlim(left, right)

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_EGU_v2.png'),dpi=300)
#plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_EGU_v2.pdf'),dpi=300)







#%% analyse first diff of lake level bias to see when the two are becoming v different
# Fig C10 paper and response letter
#=============================================================================

from numpy import diff
from matplotlib.dates import DateFormatter

flag_version = 'paper' #paper or response

# set timebounds
if flag_version == 'response':
    left= pd.to_datetime('2003-01-01')
    right= pd.to_datetime('2010-01-01') # to get nicer labels 
    freq_label='6MS'
elif flag_version == 'paper':
    left= pd.to_datetime('2018-01-01')
    right= pd.to_datetime('2021-01-01') # to get nicer labels   
    freq_label='2MS' 

# Data (a)
df_obs = LL_obs.loc[str(left):str(right)]['water_level']
df_mod = WBM_obs.loc[str(left):str(right)]['L_obs']

# Data (b)
residual = df_mod - df_obs 
df_res = pd.DataFrame(residual)
df_res.columns = ['residual']
df_res_sm = df_res.rolling(80, center=True).mean() #smooth
badrows = np.any(np.isnan(df_res_sm), axis=1) #if along the row any of the cells has a nan
df_res_sm = df_res_sm[~badrows] #remove Nans

# Data (c)
df_res_diff = df_res_sm['residual'].diff()

# plot
fig = plt.figure(figsize=(7,9))

#(a)
ax1 = plt.subplot(3, 1, 1)
ax1.plot(df_obs,label="observed", color=col_obs) #(1983-2020)
ax1.plot(df_mod, label="modelled ", c=col_mod ) #(observational run)
ax1.set_ylabel("lake level (m a.s.l.)")
ax1.grid(True)
plt.legend()
plt.title('(a)', loc='left')

#(b)
c='C0'
ax2 = plt.subplot(3, 1, 2, sharex=ax1)
df_res_sm['residual'].plot(ax=ax2, legend=False, label="bias (smoothed)", c=c)
plt.grid(True)
plt.legend()
ax2.set_ylabel('lake level bias (m)')
plt.title('(b)', loc='left')

#(c)
#c='C0'
ax3 = plt.subplot(3, 1, 3, sharex=ax1)
df_res_diff.plot(ax=ax3, legend=False, label='rate of change', c=c)
plt.grid(True)
ax3.set_ylabel('lake level bias rate of change (m day$^{-1}$)')
plt.legend()
plt.title('(c)', loc='left')

ax3.set_xlabel('')
xlist =  pd.date_range(str(left), end=str(right), freq=freq_label) # 2 MS for original SI figure
ax3.xaxis.set_ticks(xlist)
ax3.xaxis.set_ticklabels(xlist.strftime("%b-%Y"), rotation = 30, ha="right", rotation_mode="anchor")

fig.tight_layout()

if flag_save == True:
    plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_bias_newcolors.png'),dpi=300)
    plt.savefig(os.path.join(fig_path,'LL_WBMPy_evaluation_lakelevs_bias_newcolors.pdf'),dpi=300)









#%% 
"""
PART 2: Model evaluation subpanels (Fig 9)
""" 

# obs green, mod blue
#col_obs = 'C2'
#col_mod = 'C0'


# try also: obs blue, mod green or orange


#%%
# NOT IN PAPER 
#=============================================================================

#data 
left= pd.to_datetime('1983-01-01')
right= pd.to_datetime('2021-01-01') # to get nicer labels 
df_obs = LL_obs.loc[str(left):str(right)]['water_level']
df_mod = WBM_obs.loc[str(left):str(right)]['L_obs']

residual = df_mod - df_obs 
df_res = pd.DataFrame(residual)
df_res.columns = ['residual']


# make data bias levels Fig 36 a)
df_roll = df_res.rolling(180, center=True).mean()

# bias dL/dt daily Fig 36 b)
WBM_diff = WBM_obs['L_obs'].diff()
LL_diff = LL_obs.diff()
diff_res = WBM_diff - LL_diff['water_level'].loc['1983-01-01':'2020-12-31 ']
diff_roll = diff_res.rolling(365, center=True).mean()

#diff_month = diff_res.resample('M', label='right').mean()
diff_climato = diff_res.groupby(diff_res.index.dayofyear).mean()
diff_climato_roll = diff_climato.rolling(30, center=True).mean()


# Plot Figure 36
fig, axes = plt.subplots(1,2, figsize=(10,4)) 

# Panel a) lake levels 1983-2020
ax=axes[0]
ax.fill_between(df_roll.index, df_roll['residual'], where=(df_roll['residual'] <= 0), color='red', alpha=0.5, 
                interpolate=True)
ax.fill_between(df_roll.index, df_roll['residual'], where=(df_roll['residual'] > 0), color='blue', alpha=0.5, 
                interpolate=True)
#ax.set_ylabel('$\Delta$level [m]')
ax.set_ylabel('lake level bias (m)')
ax.set_xlim(left, right)
#ax.set_title('(a) Lake level bias, 180-day rolling window ', loc='left') # 
#ax.set_title('(a) Lake level bias 1983$-$2021 ', loc='left') # 
ax.set_title('(a) ', loc='left') # 

# Panel b) daily diff climatology
ax=axes[1]
ax.fill_between(diff_climato_roll.index, diff_climato_roll, where=(diff_climato_roll <= 0), color='red', alpha=0.5, 
                interpolate=True)
ax.fill_between(diff_climato_roll.index, diff_climato_roll, where=(diff_climato_roll > 0), color='blue', alpha=0.5, 
                interpolate=True)
plt.xlim(0, 365)
dates = pd.date_range(start='2014-01-01', end='2014-12-31')
ax.xaxis.set_major_locator(mdates.MonthLocator())
fmt = mdates.DateFormatter('%b')
ax.xaxis.set_major_formatter(fmt)
plt.setp(ax.get_xticklabels(), rotation=30);
#plt.title('b) Climatology bias, 30-day rolling window', loc='left')
#plt.title('(b) Climatology bias', loc='left')
plt.title('(b)', loc='left')
#plt.ylabel('$\Delta$ ($\Delta$L/$\Delta$t) [m/d]')
plt.ylabel('$\Delta$L/$\Delta$t bias (m/day)')
fig.tight_layout()


#plt.savefig(os.path.join(fig_path,'WBM_bias_LL_dLdt_fill.pdf'),dpi=300)

#%% plot residual not smoothed
# NOT IN PAPER 
#=============================================================================


plt.fill_between(diff_roll.index, diff_roll, where=(diff_roll <= 0), color='red', alpha=0.5, 
                interpolate=True)
plt.fill_between(diff_roll.index, diff_roll, where=(diff_roll > 0), color='blue', alpha=0.5, 
                interpolate=True)

plt.ylabel('$\Delta$($\Delta$L / $\Delta$t) [m/d]')
plt.xlim(left, right)
plt.title('daily dL/dt bias, rolling window')

#%%
# Paper Fig. 9 subpanels
#=============================================================================

# Calc bias of dL/dt for 180 days 

# Panel a) b) 
fig, axes = plt.subplots(2,1, figsize=(5,6), sharex=(True)) #

ax=axes[0]
LL_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 '].plot(ax=ax, label='observations', c=col_obs)
WBM_obs_dLdt['dLdt_180'].plot(ax=ax, label='modelled', c=col_mod)
ax.legend()
ax.set_ylabel(' $X$ [m]') #$\Delta$L/$\Delta$t for $\Delta$t=180 days
ax.set_ylabel('$\Delta$L/$\Delta$t (m)')
ax.set_title('(a)', x=-.1, fontsize=16)
ax.xaxis.set_tick_params(which='both', labelbottom=True)

#xlist = pd.to_datetime([1983] + list(range(1990,2020,5)) + [2021], format='%Y')
#ax.xaxis.set_ticks(xlist)
#ax.xaxis.set_ticklabels(xlist.strftime("%Y"))

ax=axes[1]
diff_dLdt = (WBM_obs_dLdt - LL_obs_dLdt.loc['1983-01-01':'2020-12-31 '])['dLdt_180']
data = diff_dLdt.rolling(3).mean()
ax.fill_between(data.index, data, where=(data <= 0), color='red', alpha=0.5, 
                interpolate=True)
ax.fill_between(data.index, data, where=(data > 0), color='blue', alpha=0.5, 
                interpolate=True)
#ax.set_ylabel(' bias $\Delta X$ [m]') 
ax.set_ylabel('bias $\Delta$L/$\Delta$t (m)')
ax.set_title('(b)', x=-.1, fontsize=16)


fig.tight_layout()

if flag_save == True:
    plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_tseries_newcolors.pdf'),dpi=300)
    plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_tseries_newcolors.png'),dpi=300)


#%%

# Figure 37 panel d) 

# Histogram of dLdt in WBM and obs 

plt.figure(figsize=(5.5,4))
LL_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 '].hist(density=True, alpha=0.9, bins=30, color=col_obs,label='observations')
WBM_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 '].hist(density=True, alpha=0.7, bins=30, color=col_mod,label='modelled')
#plt.xlabel('$\Delta$L/$\Delta$t for $\Delta$t=180 days [m]')
#plt.xlabel('$X$ [m]') #($\Delta$L/$\Delta$t for $\Delta$t=180 days)
plt.xlabel('$\Delta$L/$\Delta$t (m)')
plt.ylabel('density')
plt.legend()
plt.title('(d)', x=-.1, fontsize=16)

fig.tight_layout()

if flag_save == True:
    plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_hist_newcolors.pdf'),dpi=300)
    plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_hist_newcolors.png'),dpi=300)

#%% for EGU

# Histogram of dLdt in WBM and obs 

plt.figure(figsize=(5.5,4))
LL_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 '].hist(density=True, alpha=0.9, bins=30, color=col_obs,label='observations')
WBM_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 '].hist(density=True, alpha=0.7, bins=30, color=col_mod,label='modelled')
#plt.xlabel('$\Delta$L/$\Delta$t for $\Delta$t=180 days [m]')
#plt.xlabel('$X$ [m]') #($\Delta$L/$\Delta$t for $\Delta$t=180 days)
plt.xlabel('$\Delta$L/$\Delta$t (m)')
plt.ylabel('density')
plt.legend()
plt.title('(d)', x=-.1, fontsize=16)
plt.grid(False)

fig.tight_layout()
#plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_hist_EGU_v3.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_hist_EGU_v3.png'),dpi=300)




#%%

# c) Climatology of the variable X

fig, ax = plt.subplots(figsize=(5,3)) #

data1 = LL_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 ']
data1 = data1.groupby(data1.index.dayofyear).mean()
data1.plot(ax=ax, label='observations', c=col_obs)
data2 = WBM_obs_dLdt['dLdt_180'].groupby(WBM_obs_dLdt.index.dayofyear).mean()
data2.plot(ax=ax, label='modelled', c=col_mod)
# X axis 
plt.xlim(0, 365)
dates = pd.date_range(start='2014-01-01', end='2014-12-31')
ax.xaxis.set_major_locator(mdates.MonthLocator())
fmt = mdates.DateFormatter('%b')
ax.xaxis.set_major_formatter(fmt)
plt.setp(ax.get_xticklabels(), rotation=30);
# axis labels
ax.set_xlabel("")
ax.set_ylabel("$X$ [m]")
ax.set_ylabel('$\Delta$L/$\Delta$t (m)')
# Title+legend
plt.title('(c)', x=-0.1, fontsize=16) #Climatology
plt.legend()

fig.tight_layout()
if flag_save == True:
    plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_climatology_newcolors.pdf'),dpi=300)
    plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_climatology_newcolors.png'),dpi=300)

#%%

# b) Climatology of the daily variation in lake levels 

n = 30

fig, ax = plt.subplots() #

data1 = LL_diff.loc['1983-01-01':'2020-12-31 '].rolling(n, center=True).mean() # if n=180 and center=False this is the same climatologuy as X -means maximum rise is centered in March going from Nov to May (in obs)
data1 = data1.groupby(data1.index.dayofyear).mean()
data1.plot(ax=ax, label='observations', c=col_obs)
data2= WBM_diff.rolling(n, center=True).mean()
data2=data2.groupby(data2.index.dayofyear).mean()
data2.plot(ax=ax, label = 'WBM')
data3= data2 - data1['water_level']
data3.plot(ax=ax)
plt.xlim(0, 365)
dates = pd.date_range(start='2014-01-01', end='2014-12-31')
ax.xaxis.set_major_locator(mdates.MonthLocator())
fmt = mdates.DateFormatter('%b')
ax.xaxis.set_major_formatter(fmt)
plt.setp(ax.get_xticklabels(), rotation=30);
plt.title('b) Climatology daily $dL/dt$', loc='left')
plt.legend()


#%%

# c) Climatology of the daily variation in lake levels - turned to timescale of X

n = 180

fig, ax = plt.subplots() #

data1 = LL_diff.loc['1983-01-01':'2020-12-31 '].rolling(n, center=True).mean() # if n=180 and center=False this is the same climatologuy as X -means maximum rise is centered in March going from Nov to May (in obs)
data1.groupby(data1.index.dayofyear).mean().plot(ax=ax, label='observations')
data2= WBM_diff.rolling(n, center=True).mean()
data2.groupby(data2.index.dayofyear).mean().plot(ax=ax, label = 'modelled')
plt.xlim(0, 365)
dates = pd.date_range(start='2014-01-01', end='2014-12-31')
ax.xaxis.set_major_locator(mdates.MonthLocator())
fmt = mdates.DateFormatter('%b')
ax.xaxis.set_major_formatter(fmt)
plt.setp(ax.get_xticklabels(), rotation=30);
plt.title('b) Climatology daily $dL/dt$', loc='left')
plt.legend()



#%%

# d) climatology of lake levels 

n = 3

fig, ax = plt.subplots() #

data1 = LL_obs['water_level'].loc['1983-01-01':'2020-12-31 '].rolling(n, center=True).mean()
data1.groupby(data1.index.dayofyear).mean().plot(ax=ax, label='observations')

data2= WBM_obs['L_obs'].rolling(n, center=True).mean()
data2.groupby(data2.index.dayofyear).mean().plot(ax=ax, label = 'WBM')

plt.xlim(0, 365)
dates = pd.date_range(start='2014-01-01', end='2014-12-31')
ax.xaxis.set_major_locator(mdates.MonthLocator())
fmt = mdates.DateFormatter('%b')
ax.xaxis.set_major_formatter(fmt)
plt.setp(ax.get_xticklabels(), rotation=30);
plt.title('c) Climatology lake levels')
plt.legend()

#%% 

# scatterplot of X in observations v WBM

data1 = LL_obs_dLdt['dLdt_180'].loc['1983-01-01':'2020-12-31 ']
data2 = WBM_obs_dLdt['dLdt_180']
data = pd.DataFrame({'observations':data1,
                     'WBM':data2})

#%%

sns.kdeplot(x=data.observations, y=data.WBM, cmap="Reds", shade=True)
plt.show()

#%%

#lin reg model 

X = data.observations
Y = data.WBM
X = sm.add_constant(X)
model = sm.OLS(Y,X)
results = model.fit()
print( results.params)
#print(results.tvalues)
#print(results.t_test([2, 0]))

#%%

fig, ax = plt.subplots()

g = sns.jointplot(x=data.observations, y=data.WBM, cmap="Reds", kind='kde',
                  marginal_kws=dict(fill=True,color='orangered'),
                  joint_kws=dict(fill=True),
                  height=5)

# change labels
g.ax_joint.set_xlabel('observations, $\Delta$L/$\Delta$t (m)')
g.ax_joint.set_ylabel('WBM,  $\Delta$L/$\Delta$t (m)')

# plot 1-1 line
x0, x1 = g.ax_joint.get_xlim()
y0, y1 = g.ax_joint.get_ylim()
lims = [max(x0, y0), min(x1, y1)]
g.ax_joint.plot(lims, lims, '-r', label='x=y')

# plot regression line
xvals = np.linspace(lims[0],lims[1])
yvals = xvals * results.params[1] + results.params[0]
g.ax_joint.plot(xvals, yvals, '--r', label='OLS regression')

#add legend for lines
legendMain=g.ax_joint.legend(loc='upper left')

# title 
plt.suptitle('(e)',x=0,fontsize=16)

#fig.tight_layout()
#plt.savefig(os.path.join(fig_path,'WBM_bias_dLdt_jointplot.pdf'), bbox_inches='tight',pad_inches=0.1,dpi=300)











