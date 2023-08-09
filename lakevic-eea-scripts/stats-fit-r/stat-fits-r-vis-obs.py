# -*- coding: utf-8 -*-
"""
Created on Fri May 19 13:10:24 2023

@author: rpietroi
"""

import sys, os 
import numpy as np
import statsmodels
import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy
from matplotlib.legend_handler import HandlerTuple

from scipy import stats
from scipy.stats import genextreme as gev
import scikits.bootstrap as boot 


# set working directory
os.chdir(os.path.join (os.path.dirname(__file__), '..', '..' )) #os.getcwd() to check

#set figures-path
fig_path = os.path.join( "figures/figures_15_may23_gevfits-r") 
if os.path.exists(fig_path):
    print("The directory", fig_path, "exists!")
else:
    os.makedirs(fig_path)
    print("The directory", fig_path, "was made!")
    


"""
Plotting functions
"""


"""
1. Observations
"""

fit_path = os.path.join( "output/stat-fits-output/obs/lakevic-fit-obs-gmst-1897-2020.csv") 
df_fit = pd.read_csv(fit_path, sep=',', index_col=0)

dLdt_path =  os.path.join( "output/event_def_dt180/obs/blockmax_1897_2021_dt180_obs.csv") 
df_dLdt = pd.read_csv(dLdt_path, sep=',', index_col=0).drop(columns='daymax_180').rename(columns={"dLdt_180": "dLdt"})

gmst_path = os.path.join("data/input-data/gmst-obs/igiss_al_gl_a_4yrlo.csv")
df_gmst = pd.read_csv(gmst_path, sep=',', index_col=0)

data = pd.merge(df_dLdt, df_gmst, left_index=True, right_index=True)
data= data.loc[data.index[0]:2020]



#%% plot w subplots 


def add_interval(ax, xdata, ydata, c=None, caps="  ", width=0.01):
    line = ax.add_line(mpl.lines.Line2D(xdata, ydata, color=c))
    left = xdata[0] - width / 2
    right = xdata[0] + width / 2
    a0 = ax.add_line(mpl.lines.Line2D([left,right], [ydata[0], ydata[0]], color=c))
    a1 = ax.add_line(mpl.lines.Line2D([left,right], [ydata[1], ydata[1]], color=c))
    return (line,(a0,a1))

# make figure
fig = plt.figure(figsize=(10,4.5))

# a) plot linear model of loc parameter 
#fig, ax = plt.subplots()
ax = plt.subplot(1, 2, 1)
ax.set_title("(a)", loc='left')

ax.scatter(data["gmst"], data["dLdt"])

intercept, slope = df_fit["params.mu0"], df_fit["params.mu1"]

gmst_1900 = df_fit["gmst_1900"][0]
gmst_2020 = df_fit["gmst_2020"][0]
xvals = [data["gmst"].min(), data["gmst"].max()]
yvals = [intercept[0] + xvals[0] * slope[0]  , intercept[0] + xvals[1] * slope[0] ] 

# plot line and 95 CI of loc in 1900 and 2020 
ax.plot(xvals, yvals, color = 'red', label='$\mu$ = $\mu_0$ + $\mu_1$ GMST \n$\mu_0$= {} ({},{})\n$\mu_1$= {} ({},{})'.format(round(intercept[0],2), round(intercept[1],2), round(intercept[2],2), round(slope[0],2), round(slope[1],3), round(slope[2],2)))
add_interval(ax, [gmst_1900, gmst_1900], [df_fit["loc_1900"][1], df_fit["loc_1900"][2]], caps="--", c='red', width=0.02)
add_interval(ax, [gmst_2020, gmst_2020], [df_fit["loc_2020"][1], df_fit["loc_2020"][2]], caps="--", c='red', width=0.02)

ax.set_xlabel("GMST anomaly wrt. 1951-1980 ($\degree$C)")
ax.set_ylabel("$\Delta$L/$\Delta$t (m)")

# modify the legend
legend = plt.legend(handleheight=11, loc='upper left', bbox_to_anchor=(0, 1.06), bbox_transform=ax.transAxes, frameon=False)

rect=mpatches.Rectangle((0.015,0.985), 0.5, -0.2, #(0.015,0.985), 0.44, -0.25, 
                        fill=True,
                       alpha=0.8,
                       facecolor="w", 
                       edgecolor='lightgray', transform=ax.transAxes, capstyle='round')
plt.gca().add_patch(rect)



# b) plot return time plot 

#fig, ax = plt.subplots()
ax = plt.subplot(1, 2, 2)
ax.set_title("(b)", loc='left')

timesteps = np.r_[2:10000:1] # this is just how many years to plot along the x axis to plot as return time plot
shape, scale, loc_1900, loc_2020 = -df_fit["params.shape"], df_fit["params.scale"], df_fit["loc_1900"], df_fit["loc_2020"] # change sign of shape for use in python 

#plot observed value
p1 = ax.hlines(df_fit['int_obs'],  0, 1e4, color='fuchsia',lw=0.8,alpha=1,ls='--', label='2020 observed value')

# calculate return periods (inverse survival function) - check this is correct combination of parameters !! 
GEV_surv1900 = gev.isf(1./timesteps, c=shape[0],loc=loc_1900[0],scale=scale[0])
GEV_surv1900_CI1 = gev.isf(1./timesteps, c=shape[2],loc=loc_1900[1],scale=scale[2]) # check which shape param is in which CI in this package # check this is right way to plot CIs ?? I plotted one that looked best to me...
GEV_surv1900_CI2 = gev.isf(1./timesteps, c=shape[1],loc=loc_1900[2],scale=scale[1])

GEV_surv2020 = gev.isf(1./timesteps, c=shape[0],loc=loc_2020[0],scale=scale[0])
GEV_surv2020_CI1 = gev.isf(1./timesteps, c=shape[2],loc=loc_2020[1],scale=scale[2])
GEV_surv2020_CI2 = gev.isf(1./timesteps, c=shape[1],loc=loc_2020[2],scale=scale[1]) # why does the CI shrink around the observed val in 2020??

# plot best estimates
p2, = ax.plot(timesteps, GEV_surv1900, c='blue', label='GEV shift fit 1900')
p3, = ax.plot(timesteps, GEV_surv2020, c='red', label='GEV shift fit 2020')

# plot CIs
ax.fill_between(timesteps, GEV_surv1900_CI1, GEV_surv1900_CI2, color='blue', alpha=0.3)
ax.fill_between(timesteps, GEV_surv2020_CI1, GEV_surv2020_CI2, color='red', alpha=0.3)
ax.set_xlim(1, 1e4)

# calc empirical return periods
data["rank"] = scipy.stats.rankdata(-data["dLdt"])
data["rp_emp"] = 1 / (data["rank"] / len(data["rank"] +1))
alph=0.8
# plot points with 1900 shift
shift = loc_1900[0] - df_fit["params.mu0"][0]
p4 = ax.scatter(data["rp_emp"], data["dLdt"] + shift, c="blue", marker="x", alpha=alph)
# plot points with 2020 shift
shift = loc_2020[0] - df_fit["params.mu0"][0]
p5 = ax.scatter(data["rp_emp"], data["dLdt"] + shift, c="red", marker="x", alpha=alph)

# label axes and legend
plt.xscale('log')
plt.xlabel('return period (years)')
plt.ylabel('$\Delta$L/$\Delta$t (m)')
plt.legend([p1, p2, p3, (p4, p5)], ['2020 observed magnitude', 'GEV shift fit 1900', 'GEV shift fit 2020', 'empirical return periods'], 
           handler_map={tuple: HandlerTuple(ndivide=None)}, loc='upper left')


fig.tight_layout()

# =============================================================================
# plt.savefig(os.path.join(fig_path,'gev-shiftfit-obs-2.pdf'),dpi=300)
# plt.savefig(os.path.join(fig_path,'gev-shiftfit-obs-2.png'),dpi=300)
# =============================================================================



#%% Extra script: test scatter points

#rps_2020 = 1 / gev.sf(data["dLdt"],  c=-df_fit["params.shape"][0], loc=loc_2020[0], scale=scale[0])

data["rank"] = scipy.stats.rankdata(-data["dLdt"])
data["rp_emp"] = 1 / (data["rank"] / len(data["rank"] +1))

#%% test figure out sign convention of shape param - is opposite to R package, this way gives same res as R
rp_1900 = 1 / gev.sf(df_fit['int_obs'],  c=-df_fit["params.shape"][0], loc=loc_1900[0], scale=scale[0])
print(rp_1900)

# test figure out plotting CIs, not 100% sure about which ones i should combine...
rp_1900_1 = 1 / gev.sf(df_fit['int_obs'],  c=-df_fit["params.shape"][2], loc=loc_1900[1], scale=scale[2])
print(rp_1900_1)


# gives same res as in R
rp_2020 = 1 / gev.sf(df_fit['int_obs'],  c=-df_fit["params.shape"][0], loc=loc_2020[0], scale=scale[0])
print(rp_2020)