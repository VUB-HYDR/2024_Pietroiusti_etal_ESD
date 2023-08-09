#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:27:18 2022

@author: rosapietroiusti

Script 5: Analysing Python WBM - look at drivers of 2020 event 

Update March 2023 

"""


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
from matplotlib.colors import ListedColormap
from matplotlib import cm
import matplotlib.ticker as ticker

#%%

#===========#
#== DATA  ==#
#===========#

workDIR = os.path.join(r"C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea" )

inDIR_Py = os.path.join(workDIR,r'WBM-git\output\2022\output\obs\PERSIANN\v01r01_out\\') # '/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion/output'

fig_path = os.path.join(workDIR,r'lakevic-eea-scripts\figures\figures_6_mar23_wbmdrivers') #'/Users/rosapietroiusti/Documents/Thesis/WBM_pythonscripts/figures/figures_9_may22_WBMoutput_py'



# Python version of WBM-obs
#=============================================================================
filepath = os.path.join(inDIR_Py, 'WBM_run_PERSIANN_obs_v01r01_1983_2020_mm.csv') # 'WBM_run_PERSIANN_obs_v04r01_1983_2020_mm.csv') # --> 1983-2020
WBM_obs = pd.read_csv(filepath,
                          index_col = 'date',
                          parse_dates=True)

# Settings  
#=============================================================================
ver_n = 1 #4
run_n = 1

# Settings  
#=============================================================================
name_obs = 'WBM obs (1983-2020)'
name_obs_plot = 'WBM observational run \n(1983-2020)'


startYEAR = WBM_obs.index[0].year
endYEAR = WBM_obs.index[len(WBM_obs)-1].year


#%%===========================================================================
"""
PART 1: 
    Figure 10. a) Climatology of WB terms and b-c) what they were in 2019-2020
=============================================================================
"""

# Data
df_mm = WBM_obs[["P_lake", "E_lake", "Q_in", "Q_out"]]
df_mm_resample = df_mm.resample('M', label='right').sum()
df_mm_climatology = df_mm_resample.groupby(df_mm_resample.index.month).mean() # climatology
df_mm_stdev = df_mm_resample.groupby(df_mm_resample.index.month).std()

# Calculate max and min
df_mm_max = df_mm_resample.groupby(df_mm_resample.index.month).max()
df_mm_min = df_mm_resample.groupby(df_mm_resample.index.month).min()

# plotting function
def plot_climato_year(data_climato, data_stdevs, data_monthly, year, subplot, ax, data_maxmin=0, plot_maxmin=False, plot_Plake=True,
                      plot_Qin=True,plot_Qout=True, plot_E=False, labels=False, output=True):
    if plot_maxmin == True:
        labelslist = ['mean climatology'] *3
    elif labels==True:
        labelslist = ['precipitation', 'inflow', 'outflow']
    else:
        labelslist = ['', '', '']
            
    if plot_Plake == True:
        P, = ax.plot(data_climato["P_lake"], c='#1f77b4', label=labelslist[0] )
    if plot_Qin == True:
        Qin, = ax.plot(data_climato["Q_in"], c='#2ca02c', label =labelslist[1] ) # 'inflow'
    if plot_Qout == True:
        Qout, = ax.plot(data_climato["Q_out"], c='#d62728', label =labelslist[2]) 
    if plot_E == True:
        E, = ax.plot(data_climato["E_lake"], c='#ff7f0e')    
    
    # Plot max and min
    if plot_maxmin == True:  # add options of different ones
        ax.plot(data_maxmin["P_lake"], c='#1f77b4', ls=':', label='records' )
        ax.plot(data_maxmin["Q_in"], c='#2ca02c', ls=':', label='records' )
        ax.plot(data_maxmin["Q_out"], c='#d62728', ls=':')
    if plot_maxmin == True:
        ax.plot(data_maxmin["P_lake"], c='#1f77b4', ls=':')
        ax.plot(data_maxmin["Q_in"], c='#2ca02c', ls=':')
    
        # Plot 2019 or 2020 
    datay = data_monthly.loc['{}-01-01'.format(year):'{}-12-31'.format(year)].reset_index(drop=True)
    datay.index = datay.index + 1
    if plot_Plake == True:
        Py, = ax.plot(datay["P_lake"], c='#1f77b4', ls='--', label='2019 / 2020'  )
    if plot_E == True:
        Ey, = ax.plot(datay["E_lake"], c='#ff7f0e', ls='--') 
    if plot_Qin == True:
        Qiny, = ax.plot(datay["Q_in"], c='#2ca02c', ls='--') 
    if plot_Qout == True:
        Qouty, = ax.plot(datay["Q_out"], c='#d62728', ls='--') 

    # Shade standard devs
    elslist = []
    colorlist = []
    plotstd = [None] *4
    if plot_Plake == True:
        elslist.append("P_lake")
        colorlist.append('#1f77b4')
    if plot_Qin == True:
        elslist.append("Q_in")
        colorlist.append('#2ca02c')
    if plot_Qout == True:
        elslist.append("Q_out")
        colorlist.append('#d62728')
    if plot_E == True:    
        elslist.append("E_lake")
        colorlist.append('#ff7f0e')
    for i,k, color in zip(elslist, range(len(elslist)), colorlist): # "P_lake", "E_lake", , "Q_out"
        print(i)
        y = data_climato[i]
        error = data_stdevs[i]
        plotstd[k] = ax.fill_between(data_climato.index, y-error, y+error, alpha = 0.3, color=color)
        #return i
        
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(which='minor', width=0.75, length=2.5)
    ax.set_title(subplot, loc='left')
    #ax.set_ylabel('l.l.e. (mm/month)')
    ax.set_xlabel('months')
    
    if output==True:
        if plot_E == True:
            return P,Qin,Qout,E,Py,Qiny,Qouty,Ey,plotstd 
        else:
            return P,Qin,Qout,Py,Qiny,Qouty,plotstd 


#%%
#===========================================================================
# Plot
#=============================================================================

fig = plt.figure(figsize=(9,6),dpi=300)  # initiate the figure figsize = (15, 12), figsize=(10,6),

# a ) 
ax = plt.subplot2grid((2, 4), (0, 1), colspan=2, rowspan=1)

#stuff = ax.plot(df_mm_climatology["P_lake", "E_lake", "Q_in", "Q_out"],label=['precipitation', 'evaporation', 'inflow', 'outflow'])

stuff = df_mm_climatology.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"], ax=ax, label=['precipitation', 'evaporation', 'inflow', 'outflow'], legend=None)


for i in ["P_lake", "E_lake", "Q_in", "Q_out"]:
    y = df_mm_climatology[i]
    error = df_mm_stdev[i]
    ax.fill_between(df_mm_climatology.index, y-error, y+error, alpha = 0.3)

#ax.legend(loc="upper right", bbox_to_anchor=(1.7, 1.03))

# Label x axis 
ax.set_title('(a) ', loc='left')
#ax.set_ylabel('l.l.e. (mm/month)')
ax.set_ylabel('l.l.e. (mm month$^{-1}$)')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
#ax.tick_params(which='minor', width=0.75, length=2.5)
ax.set_xlabel("")
ax.set_xticks(range(1,13,1))
ax.set_xticklabels(pd.date_range('2019-01-01','2019-12-31', 
              freq='MS').strftime("%b").tolist(), fontsize=9, rotation=45,ha='right', rotation_mode='anchor')

# b ) 
ax1 = plt.subplot2grid((2,4), (1, 0), colspan=2, rowspan=1)
P,Qin,Qout,Py,Qiny,Qouty,plotstd = plot_climato_year(df_mm_climatology, df_mm_stdev, df_mm_resample,  2019, '(b) 2019', ax1)
#P,Qin,Qout,E,Py,Qiny,Qouty,Ey,plotstd = plot_climato_year(df_mm_climatology, df_mm_stdev, df_mm_resample,  2019, '(b) 2019', ax1, plot_E=True)
#ax1.set_ylabel('l.l.e. (mm/month)')
ax1.set_ylabel('l.l.e. (mm month$^{-1}$)')
ax1.set_xlabel("")
ax1.set_xticks(range(1,13,1))
ax1.set_xticklabels(pd.date_range('2019-01-01','2019-12-31', 
              freq='MS').strftime("%b").tolist(), fontsize=9, rotation=45,ha='right', rotation_mode='anchor')

# c ) 
ax2 = plt.subplot2grid((2, 4), (1, 2), colspan=2, rowspan=1, sharey = ax1)
plot_climato_year(df_mm_climatology, df_mm_stdev, df_mm_resample,  2020, '(c) 2020', ax2)
#ax2.set_ylabel('l.l.e. (mm/month)')
ax2.set_ylabel('l.l.e. (mm month$^{-1}$)')
ax2.set_xlabel("")
ax2.set_xticks(range(1,13,1))
ax2.set_xticklabels(pd.date_range('2020-01-01','2020-12-31', 
              freq='MS').strftime("%b").tolist(), fontsize=9, rotation=45,ha='right', rotation_mode='anchor')

# legend (turn on if I plotE=True to make nice legend)
#lgnd = fig.legend([(P,plotstd[0]),(Qin,plotstd[1]),(Qout,plotstd[2]), (E,plotstd[3])], 
#                  ["precipitation","inflow","outflow", "evaporation"], title = "Climatology",
#                  bbox_to_anchor=(.95,.95),frameon=False)
#lgnd._legend_box.align = "left"

lgnd2 = fig.legend([Py,Qiny,Qouty], 
                  ["precipitation","inflow","outflow"], title = "2019$-$2020 values", 
                  bbox_to_anchor=(.95,.75),frameon=False) #(1,.95)
lgnd2._legend_box.align = "left"

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_climatology_std_years_noE_v3.pdf'.format(startYEAR, endYEAR)),dpi=300, transparent=False)
#plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_climatology_std_years.png'.format(startYEAR, endYEAR)),dpi=300, transparent=False)



#%% 
"""
==========================================================================
 PART 2: Supplementary figure accumulated anomalies 
 update 2023
=============================================================================
"""

# Data
# ==============================================================

# (a) anomaly monthly 2019-2020
df_mm_climatology_2y = df_mm_climatology.append(df_mm_climatology).reset_index(drop=True)
df_mm_climatology_2y.index = df_mm_climatology_2y.index +1

df_anomaly = df_mm_resample.loc['{}-01-01'.format(2019):'{}-12-31'.format(2020)].reset_index(drop=True)
df_anomaly.index = df_anomaly.index +1
df_anomaly = df_anomaly - df_mm_climatology_2y

# (b)  accumulated anomaly over n months (data = WBM_obs)

def plot_rollingterms(data, n): 
    # data is monthly accumulations! n is number of months over which you want the accumulation calc'd
    # First i resample to monthly (sum) and then I take a 12-month moving window cumulation (sum), 
    # assign the value to the final month, so its a cumulative anomaly in the previous n months

    df = data[["P_lake", "E_lake", "Q_in", "Q_out"]].resample('M', label='right').sum().rolling(n, center=False).sum()
    df_residual = data["P_lake"] - data["E_lake"] + data["Q_in"] - data["Q_out"]
    df_res_rolling = df_residual.resample('M', label='right').sum().rolling(n, center=False).sum()
    df['res'] = df_res_rolling
    
    y0 = df["P_lake"] / 1000
    y1 = (df["P_lake"] + df["Q_in"]) / 1000
    y3 = df["E_lake"] / 1000
    y4 = (df["E_lake"] + df["Q_out"]) / 1000
    res = df["res"] / 1000

    ax.fill_between(df.index, 0, y0, alpha=0.9, edgecolor='C0', hatch='///', label='precipitation')
    ax.fill_between(df.index, y0, y1, alpha=0.9, color='C2', hatch='///', label='inflow' ) # inflow
    ax.fill_between(df.index, 0, -y3, alpha=0.8, color='C1', hatch='///', label='evaporation') # evaporation
    ax.fill_between(df.index, -y3, -y4, alpha=0.8, color='C3', hatch='///', label='outflow') # outflow

    ax.fill_between(df.index, res, where=(res > 0), color='blue', alpha=0.9, 
                    interpolate=True, label='pos residual') 
    ax.fill_between(df.index, res, where=(res <= 0),  alpha=0.9, 
                    interpolate=True, color='red', label='neg residual') 


# Plot 
#=============================================================================


fig = plt.figure(figsize=(9,9))

# panel a) anomaly monthly 2019-2020

ax = plt.subplot(3,1,1)
ax.plot(df_anomaly["P_lake"], marker='o', label='precipitation')
ax.plot(-df_anomaly["E_lake"], marker='o', label='evaporation')
ax.plot(df_anomaly["Q_in"], marker='o', label='inflow')
ax.plot(-df_anomaly["Q_out"], marker='o',label='outflow')
ax.set_title('(a) ',loc='left', fontsize=14)
ax.set_title('(a) 2019-2020 monthly',loc='left', fontsize=14)
ax.set_xlim(0.5,24.5)
ax.set_ylim(-70,200)
ax.set_xticks(range(1,25,1))
ax.set_ylabel('anomaly (mm/month l.l.e)')
ax.set_xlabel("month")

ax.legend()
ax.grid()

# panel b) accum over 6 months

ax = plt.subplot(3,1,2)
plot_rollingterms(WBM_obs, 6)

ax.set_title('(b)',loc='left', fontsize=14)
ax.set_title('(b) 1983-2020 6-monthly',loc='left', fontsize=14)

ax.margins(x=0)
ax.legend(loc="upper center", ncol =6, bbox_to_anchor=(.5,-.1))
ax.set_ylabel('input/output (m/6 months l.l.e)')
years=mdates.YearLocator(1)
ax.xaxis.set_minor_locator(years)

# panel c) accum over 12 months 

ax = plt.subplot(3,1,3)
plot_rollingterms(WBM_obs, 12)

ax.set_title('(c)',loc='left', fontsize=14)
ax.set_title('(c) 1983-2020 12-monthly',loc='left', fontsize=14)

ax.margins(x=0)
ax.legend(loc="upper center", ncol =6, bbox_to_anchor=(.5,-.1))
ax.set_ylabel('input/output (m/year l.l.e)')
years=mdates.YearLocator(1)
ax.xaxis.set_minor_locator(years)

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_anomalies_v1.pdf'.format(startYEAR, endYEAR, ver_n, run_n)),dpi=300, transparent=False)

#%% 

"""
==========================================================================
 PART 3: Supplementary figure accumulated anomalies only 2019-2020
 update 2023
=============================================================================
"""

df_residual_1920 = df_anomaly["P_lake"] - df_anomaly["E_lake"] + df_anomaly["Q_in"] - df_anomaly["Q_out"]

fig , axes = plt.subplots(2, 1,figsize=(6,6))

ax = axes[0]
ax.plot(np.cumsum(df_anomaly["P_lake"]), marker='o', label='precipitation')
ax.plot(np.cumsum(-df_anomaly["E_lake"]), marker='o', label='evaporation')
ax.plot(np.cumsum(df_anomaly["Q_in"]), marker='o', label='inflow')
ax.plot(np.cumsum(-df_anomaly["Q_out"]), marker='o',label='outflow')
ax.legend()
ax.grid()
#ax.set_xticks(list(range(1,25,2)) + [24])
ax.set_xticks(range(1,25,1))
ax.set_xticklabels([])
ax.set_xlim(0.5,25)
ax.set_title('(a)',loc='left')
ax.set_ylabel('cumulative anomaly (mm l.l.e)')

ax = axes[1]
ax.plot(np.cumsum(df_residual_1920), label='residual anomaly',marker='o',c='gray')
ax.legend()
#ax.set_xticks(range(1,25,2))
#ax.set_xticks(list(range(1,25,2)) + [24])
ax.set_xticks(range(1,25,1))
ax.set_xlim(0.5,25)
ax.set_xticklabels(pd.date_range('2019-01-01','2020-12-31', 
              freq='MS').strftime("%b-%Y").tolist(), fontsize=8, rotation=45,ha='right', rotation_mode='anchor')

ax.grid()
ax.set_title('(b)',loc='left')
ax.set_ylabel('cumulative residual (mm l.l.e)')

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'WBM_terms_accumulated_anomalies_2019-2020.pdf'),dpi=300, transparent=False)

#%% numbers
print(np.cumsum(df_anomaly["P_lake"]))
print(np.cumsum(df_residual_1920))


#%%

"""

Extra plots and calculations
==========================================================================
"""

#===========================================================================
# Extra: Climatology plot with data
#=============================================================================


# Plot
fig, ax = plt.subplots(figsize=(6,4)) #, figsize=(14,4 )

# Plot
df_mm_climatology.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"], ax=ax)

summarystats = pd.DataFrame(round(df_mm_climatology.mean(),1))
summarystats.columns = ['']

ax.text(13, 200, 'monthly mean \nclimatology \n(mm/month):')
ax.text(13, 200, str(summarystats),
          horizontalalignment='left',
          verticalalignment='top')
ax.set_ylabel('mm/month')
#mean_bias = df_mm_climatology.mean()[0] - df_mm_climatology.mean()[1] + df_mm_climatology.mean()[2] - df_mm_climatology.mean()[3]
#ax.text(13, 5,'{}'.format(newname_plot))
ax.grid(True)
plt.title('water balance terms: climatology')

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'WBM_terms_Pythonv0{}r0{}_{}_{}_climatology_EGU.png'.format(ver_n, run_n, startYEAR, endYEAR)),dpi=300, transparent=False)



#%% ==========================================================================
# Extra: Plot yearly acc for EGU22
#=============================================================================

# Data
df_mm_resampley = WBM_obs[["P_lake", "E_lake", "Q_in", "Q_out"]].resample('Y', label='right').sum()

# Choose data - careful you are over-writing the monthly resampled dataframe 
#df_mm_resample = df_mm_resampley

# Plot
fig, ax = plt.subplots(figsize=(6,4)) #, 
df_mm_resampley.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"], ax=ax)
df_residual = df_mm_resampley["P_lake"] - df_mm_resampley["E_lake"] + df_mm_resampley["Q_in"] - df_mm_resampley["Q_out"]
df_residual.plot(color='lightgray', ax=ax, label = 'residual' )
pos_x=51
ax.text(pos_x, 1700, 'yearly mean \n(mm/year):')
summarystats = pd.DataFrame(round(df_mm_resample.mean(),1))
summarystats.columns = ['']
ax.text(pos_x, 1680, str(summarystats),
          horizontalalignment='left',
          verticalalignment='top')
ax.set_ylabel('mm/year')
#mean_bias = df_mm_resample.mean()[0] - df_mm_resample.mean()[1] + df_mm_resample.mean()[2] - df_mm_resample.mean()[3]
#ax.text(pos_x, 600,'{} v0{}r0{} \nmean bias {:.2f} mm/yr'.format(newname,ver_n, run_n, mean_bias))
ax.grid(True)
plt.legend(loc='lower left')
plt.axhline(y=0, color='gray', ls='--')
#plt.axhline(y=1505.9, color='gray', ls='--')
plt.title('water balance terms: yearly accumulated')

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_Pythonv0{}r0{}_Y_res_EGU_v2.png'.format(startYEAR, endYEAR, ver_n, run_n)),dpi=300, transparent=False)

#%% Same but cleaner

# Plot
fig, ax = plt.subplots(figsize=(7,4)) #, 
df_mm_resampley.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"], ax=ax, label=['precipitation', 'evaporation', 'inflow', 'outflow'])
df_residual = df_mm_resampley["P_lake"] - df_mm_resampley["E_lake"] + df_mm_resampley["Q_in"] - df_mm_resampley["Q_out"]
df_residual.plot(color='lightgray', ax=ax, label = 'residual' )

ax.set_ylabel('mm/year (l.l.e.)')
ax.grid(True)
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.axhline(y=0, color='gray', ls='--')

fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_Pythonv0{}r0{}_Y.png'.format(startYEAR, endYEAR, ver_n, run_n)),dpi=300, transparent=False)

#%%


# Plot yearly acc anomalies - maybe as little subplots?

df_mm_anomalyy = df_mm_resampley - df_mm_resampley.mean()
df_mm_anomalyy.plot()

#%%
# Plot 6mo anomaly / yearly anomaly with rolling window

n=6

WBM_rolling6 = WBM_obs[["P_lake", "E_lake", "Q_in", "Q_out"]].resample('M', label='right').sum().rolling(n, center=False).sum()
#WBM_rolling6.plot()
df_residual = WBM_obs["P_lake"] - WBM_obs["E_lake"] + WBM_obs["Q_in"] - WBM_obs["Q_out"]
df_res_rolling6 = df_residual.resample('M', label='right').sum().rolling(n, center=False).sum()
#df_res_rolling6.plot()
WBM_rolling6['res'] = df_res_rolling6

df_anomaly6_roll = ( WBM_rolling6 - WBM_rolling6.mean()) / WBM_rolling6.mean()
df_anomaly6_roll.plot()

df_anomaly6_abs = WBM_rolling6 - WBM_rolling6.mean()
df_anomaly6_abs.plot()

#%% analysis: 
# =========================

# seasonal mean and anomaly for the accumulated rolling amount thing

WBM_roll_climato = WBM_rolling6.groupby(WBM_rolling6.index.month).mean()
WBM_roll_climato.plot()

data1 = WBM_rolling6.loc['2020-01-01':'2020-12-31'].reset_index(drop=True)
data2 = WBM_roll_climato.reset_index(drop=True)

df_anomaly_climato = data1 - data2

df_anomaly_climato_rel = df_anomaly_climato / WBM_roll_climato.reset_index(drop=True)

