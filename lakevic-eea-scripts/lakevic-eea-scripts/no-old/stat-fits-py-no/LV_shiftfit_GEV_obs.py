# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 13:47:31 2023

@author: rpietroi

what this script does:
    1. fit gev to dLdt from observations (1949-2021) and see how it compares to KX -> stationary - ok
    2. do non-stationary shift fit and see how it compares to KX results -> similar for 1949-2020 incl - ok
    3. calculate PR and delta-I based on these results and their CI, compare to KX -> the CI are not the same, understand how to calc these!
        
extension:
    1. do this to longer timeseries, with data from late 1800s, including and excluding 2020/2021
    2. do this to WBM dLdt (obs and GCM-deriveds) and longer GCM results to 2100 
    3. KS statistics for observations and model results
    4. IOD analysis, regress out the IOD (results are unstable otherwise)
    5. autocorrelation timescale on dLdt

credit: Lucy Temple for her attribution and shift fit code

questions
    1. why detrend the timeseries compared to GMST? I think it has to do with the shift fitting? Try without this
    2. for Clair - how to calculate PR CI and delta I CI
    3. dealing with IOD, look at how they dealt with ENSO and do the same (covariate? regress it out? causal effect as regression coefficents)

notes:
    1. it seems I accidentally included 2020 in my GEV fit in the KX ! Slightly neg trend if I exlclude it. See if including extended timeseries OR accounting for IOD improves this ! 


"""


import sys, os 
import numpy as np
import statsmodels
import pandas as pd 
import matplotlib.pyplot as plt
import xarray as xr

from scipy import stats
from scipy.stats import genextreme as gev
import scikits.bootstrap as boot 

#%%

# set working directory
os.chdir(os.path.dirname(__file__)) #os.getcwd() to check

#set figures-path
fig_path = os.path.join('..', '..', "figures/figures_11_apr23_gevfits") 
if os.path.exists(fig_path):
    print("The directory", fig_path, "exists!")
else:
    os.makedirs(fig_path)
    print("The directory", fig_path, "was made!")

#gmst data 
gmst_path = os.path.join('..', '..', "data/input-data/gmst-obs/igiss_al_gl_a_4yrlo.nc") 

#observational dLdt data 
data_path = os.path.join('..', '..', "output/event_def_dt180/obs/for_KX/blockmax_1948_2022_dt180_obs.txt") 

print('\nthe gmst and data paths are \n{} \n{} \n'.format(gmst_path,data_path))

    
#%%

#open gmst as dataarray 1879-2021 (does this also have 2022?)
gmst_raw = xr.open_dataarray(gmst_path, decode_times=False)
gmst_raw = gmst_raw.assign_coords(time=np.arange(1879, 1879+len(gmst_raw)))
gmst_raw = gmst_raw.sel(time=slice(1879, 2021))
gmst_raw = gmst_raw.rename({'time': 'year'})
gmst_raw.plot()

# open dLdt from obs 1949-2021
data_raw = pd.read_csv(data_path, sep='\t', index_col=0)
plt.scatter(data_raw.index, data_raw, label='annual block max dLdt 1949-2021', c='C1') # need to add a second axis!

plt.legend()

#%%

data_raw.hist(bins=30)

#%% 1) stationary GEV of data including/ excluding 2020 and 2021

startyear = 1949
endyear = 2020

data = data_raw.loc[startyear:endyear]

# fit stationary GEV to data 
parameters = gev.fit(data) #-shape, loc, scale - similar to results from CX # -shape = -0.025866 loc = 0.27414 scale = 0.21870
#print('data is from {}-{} \n-shape, loc, scale stationary GEV non detrended tseries\n'.format(startyear,endyear), *parameters)

print('data is from {}-{} \n-shape, loc, scale params stationary GEV non detrended tseries'.format(startyear,endyear), 
      tuple([round(x,4) if isinstance(x, float) else x for x in parameters]))

# from 1st percentile of theoretical GEV to 99th percentile of theoretical GEV, 100 steps
x = np.linspace(gev.ppf(0.01, *parameters),
                gev.ppf(0.99, *parameters), 100)

fig, ax = plt.subplots(1, 1)
rv = gev(*parameters)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='fitted gev pdf')

ax.hist(data, density=True, histtype='stepfilled', alpha=0.2, bins=20, label='original data')
ax.legend(loc='best', frameon=False)
#plt.text(0,0, ['-shape, loc, scale\n', *parameters])
plt.show()


# calculate 2020 return period with stationary fit
rp2020_stat = 1 / gev.sf(data_raw.loc[2020], *parameters) # survival function is 1 - cdf
print('\nreturn period of 2020 event based on observations and stationary non detrended fit {}-{}'.format(startyear,endyear), rp2020_stat)


#%% 2) try to do non-stationarity  of loc parameter

'''
1. Linear fit between GMST = X and dLdt = Y
2. get loc_2020 from GMST_2020 and loc_1900 from GMST_1900
2. Do GEV fit with loc parameter from both of the above
'''

gmst = gmst_raw.sel(year=slice(1949, endyear)).to_dataframe() # cut to same period as dLdt
data = data_raw.loc[1949:endyear]

plt.scatter(gmst,data)
plt.xlabel('gmst anomaly')
plt.ylabel('dLdt blockmax')
plt.show()

#%%

lm = stats.linregress(gmst.iloc[:,0], data.iloc[:,0]) #.values
slope, intercept = lm[0], lm[1]
mu = slope * gmst + intercept
mu = mu.rename(columns={'Ta': 'mu'}) #to_dataframe().

mu_2020 =  slope * gmst_raw.loc[2020] + intercept

mu_2 = slope * gmst
mu_2 = mu_2.rename(columns={'Ta': 'mu'}) #.to_dataframe()

#%%

plt.scatter(gmst,data)
#plt.plot(gmst, mu, c='red')
#plt.axline([0, 0], slope=1, color='black')
x_vals = [ gmst['Ta'].min() , gmst['Ta'].max() ]
y_vals = [mu[gmst['Ta'] == gmst['Ta'].min()]['mu'] , mu[gmst['Ta'] == gmst['Ta'].max()]['mu']]
plt.plot(x_vals, y_vals, color = 'red', label='linear regression gmst and dLdt \nyears 1949 - {} \nslope = {} m/'.format(endyear, round(slope,4)) )
plt.xlabel('gmst anomaly')
plt.ylabel('dLdt blockmax')
plt.legend()
plt.show()


#mu_factual = slope * 1.2 + intercept
#mu_counter_factual = slope * 0 + intercept
# not correct ! this isnt exactly the mu, what I will use is the slope of the linreg ! 

#%%

'''
GEV shift fit with Lucy's method
============================

including data 1949 - 2020 (i.e. including the 2020 event, which is a mistake)
I absolutely need to condition on IOD for this to make any sense at all !!
The timeseries of dLdt to GMST shows no linear trend at all 

'''

gmst_df = gmst_raw.to_dataframe()
gmst = gmst_df.loc[1949:endyear] # # gmst anomaly 
data_slice = data_raw.loc[1949:endyear] # data I am fitting
data = pd.merge(gmst,data_slice, left_index=True, right_index=True).rename(columns={'Ta': 'gmst', 'dLdt_180':'dLdt'})

#%% linear regression for detrending data  - WHY ARE YOU DETRENDING? is it necessary to do the shift fit? 

linfit = stats.linregress(data['gmst'], data['dLdt'])
print('\nalpha parameter of linear shift fit', round(linfit[0],4)) # the 0th value stored in our linear regression is the slope/gradient - this is the same as what we did above (lm)

# detrended values = observed values - (gradient * GMST)
data['dLdt_detrend'] = data['dLdt'] - (data['gmst'] * linfit[0]) # almost no difference 

data.plot() # almost no difference  in detrended timeseries. linear shift assumption is absolutely not met. 

#%%

# GEV fit
shift_fit1900 = gmst_df['Ta'].loc[gmst_df.index == 1900]
shift_fit2020 = gmst_df['Ta'].loc[gmst_df.index == 2020]
print('\ngmst anomaly in 1900:', float(shift_fit1900))
print('gmst anomaly in 2020:', float(shift_fit2020))

# Fit GEV on detrended data and get parameters

detrend_dLdt = np.array(data['dLdt_detrend'] )

dLdtGEV = (gev.fit(detrend_dLdt )) # parameters -shape, loc, scale  -> loc is the same as stationary fit, shape and scale are different (means KX doesn't detrend your data when it does the fit)

# Test without detrending - why does this look horrible?? DO THIS AGAIN, UNDERSTAND

#dLdt =  np.array(data['dLdt'] )
#dLdtGEV = (gev.fit(dLdt )) 

print('\nshape, mu_0, scale params GEV detrended tseries best estimates {}-{} \n'.format(startyear,endyear), 
      tuple([round(x,4) if isinstance(x, float) else x for x in dLdtGEV]))

#KS test
fitKS = stats.kstest(detrend_dLdt, 'genextreme',args=(dLdtGEV[0],dLdtGEV[1],dLdtGEV[2]))
print('\n', fitKS)

#%%
#
timesteps = np.r_[2:10000:1] # this is just how many years to plot along the x axis to plot as return time plot

dLdtGEV_surv1900=gev.isf(1./timesteps,dLdtGEV[0],loc=(dLdtGEV[1]+(shift_fit1900*linfit[0])),scale=dLdtGEV[2]) # 'inverse survival function' = 1 / sf = 1 / (1 - cdf) = return period

dLdtGEV_surv2020=gev.isf(1./timesteps,dLdtGEV[0],loc=(dLdtGEV[1]+(shift_fit2020*linfit[0])),scale=dLdtGEV[2]) 

#%% confidence interval 

seedn = 1
print('\nstart bootstrapping with seed number {} \n'.format(seedn))
ci = boot.ci(detrend_dLdt, gev.fit, n_samples=100, seed=seedn) 


# turn to 1000 later (-shape, loc, scale), 95% CI, and specify a seed! 

#%%
print('\n-shape, loc, scale 95 CI \n', ci)
print('\n-shape 95 CI', ci[0,0], ci[1,0])
print('scale 95 CI',ci[0,2],ci[1,2])
print('mu_0 95 CI',ci[1,0],ci[1,1])

dLdtGEV_surv1900_CI1=gev.isf(1./timesteps,c=ci[0,0],loc=(ci[0,1]+(shift_fit1900*linfit[0])),scale=ci[0,2]) # 'inverse survival function' = 1 / sf = 1 / (1 - cdf) = return period
dLdtGEV_surv1900_CI2=gev.isf(1./timesteps,c=ci[1,0],loc=(ci[1,1]+(shift_fit1900*linfit[0])),scale=ci[1,2])


dLdtGEV_surv2020_CI1=gev.isf(1./timesteps,c=ci[0,0],loc=(ci[0,1]+(shift_fit2020*linfit[0])),scale=ci[0,2]) # 'inverse survival function' = 1 / sf = 1 / (1 - cdf) = return period
dLdtGEV_surv2020_CI2=gev.isf(1./timesteps,c=ci[1,0],loc=(ci[1,1]+(shift_fit2020*linfit[0])),scale=ci[1,2]) # 'inverse survival function' = 1 / sf = 1 / (1 - cdf) = return period


#%%

fig, ax = plt.subplots()

#plot observed value
ax.hlines(data_raw.loc[2020,'dLdt_180'],  0, 1e4, color='aqua',lw=0.8,alpha=1,ls='--', label='2020 observed value')

# plot best estimates
ax.plot(timesteps, dLdtGEV_surv1900, c='blue', label='GEV shift fit 1900')
ax.plot(timesteps, dLdtGEV_surv2020,c='red', label='GEV shift fit 2020')

# plot CIs
ax.fill_between(timesteps, dLdtGEV_surv1900_CI1,dLdtGEV_surv1900_CI2, color='blue', alpha=0.3)
ax.fill_between(timesteps, dLdtGEV_surv2020_CI1,dLdtGEV_surv2020_CI2, color='red', alpha=0.3)

# label axes
plt.xscale('log')
plt.xlabel('return period (years)')
plt.ylabel('dLdt annual max (m)')
plt.legend()

# understand why the 1-year return period can't be estimated ? 

# estimate return period of 2020 magnitude and CI (and PR?)

#%%

# calculate 2020 return period with non-stationary fit

parameters2020_BE = dLdtGEV[0], (dLdtGEV[1]+( float(shift_fit2020)*linfit[0])) , dLdtGEV[2]
parameters2020_CI1 = ci[0,0], (ci[0,1]+(float(shift_fit2020)*linfit[0])), ci[0,2]
parameters2020_CI2 = ci[1,0], (ci[1,1]+(float(shift_fit2020)*linfit[0])) , ci[1,2]

prob2020_CI = np.array( [float(gev.sf(data_raw.loc[2020], *parameters2020_BE)), float(gev.sf(data_raw.loc[2020], *parameters2020_CI1)), float(gev.sf(data_raw.loc[2020], *parameters2020_CI2)) ]) 

rp2020_CI = 1 / prob2020_CI   # survival function is 1 - cdf

print('\nreturn period of 2020 event based on observations 1949-{} in 2020 climate and CI \n'.format(endyear), rp2020_CI)

print('probabilities in 2020 climate \n', prob2020_CI)

print('loc param in 2020 climate and CI \n', dLdtGEV[1]+( float(shift_fit2020)*linfit[0]), ci[0,1]+(float(shift_fit2020)*linfit[0]), ci[1,1]+(float(shift_fit2020)*linfit[0]) )

#%%

# calculate 2019 return period with non-stationary fit

parameters1900_BE = dLdtGEV[0], (dLdtGEV[1]+( float(shift_fit1900)*linfit[0])) , dLdtGEV[2]
parameters1900_CI1 = ci[0,0], (ci[0,1]+(float(shift_fit1900)*linfit[0])), ci[0,2]
parameters1900_CI2 = ci[1,0], (ci[1,1]+(float(shift_fit1900)*linfit[0])) , ci[1,2]

prob1900_CI = np.array( [float(gev.sf(data_raw.loc[2020], *parameters1900_BE)), float(gev.sf(data_raw.loc[2020], *parameters1900_CI1)), float(gev.sf(data_raw.loc[2020], *parameters1900_CI2)) ]) 

rp1900_CI = 1 / prob1900_CI   # survival function is 1 - cdf

print('\nreturn period of 2020 event based on observations 1949-{} in 1900 climate and CI \n'.format(endyear), rp1900_CI)

print('probabilities in 1900 climate \n', prob1900_CI)

print('loc param in 1900 climate and CI \n', dLdtGEV[1]+( float(shift_fit1900)*linfit[0]), (ci[0,1]+(float(shift_fit1900)*linfit[0])), (ci[1,1]+(float(shift_fit1900)*linfit[0])) )

#%% test calculate probability ratio 

PR_test = prob2020_CI / prob1900_CI
print('\ntest calculating a PR - BE is ok, CI no!', PR_test)

#%% calculate intensity change

int2020_check = gev.isf(prob2020_CI[0], *parameters2020_BE)
#print(int2020_check) # check magnitude of 2020 event in 2020 climate (event with rp=obs_rp)
int2020_CI1 = gev.isf(prob2020_CI[0], *parameters2020_CI1)
#print(int2020_CI1)
int2020_CI2 = gev.isf(prob2020_CI[0], *parameters2020_CI2)
#print(int2020_CI2)

int2020_CI = [int2020_check, int2020_CI2, int2020_CI1]
print('\nintensity 2020',
      tuple([round(x,5) if isinstance(x, float) else x for x in int2020_CI]))



int1900 = gev.isf(prob2020_CI[0], *parameters1900_BE)
#print(int1900)
int1900_CI1 = gev.isf(prob2020_CI[0], *parameters1900_CI1)
#print(int1900_CI1)
int1900_CI2 = gev.isf(prob2020_CI[0], *parameters1900_CI2)
#print(int1900_CI2)

int1900_CI = [int1900, int1900_CI2, int1900_CI1]
print('intensity 1900',
      tuple([round(x,5) if isinstance(x, float) else x for x in int1900_CI]))

# test intensity change

deltaI_BE = int2020_CI[0] - int1900_CI[0]
#print(deltaI_BE)

deltaI_CI = [int2020_CI[0] - int1900_CI[0], int2020_CI[1] - int1900_CI[2], int2020_CI[2] - int1900_CI[1] ]
print('test deltaI, check that CI is ok', deltaI_CI )

#%%

# to do: print results, CI etc. in a tidy way and do this clean to be able to use it on all different type of data 

# to do: check calulation and BE of PROBABILITY RATIO and INTENSITY CHANGE 