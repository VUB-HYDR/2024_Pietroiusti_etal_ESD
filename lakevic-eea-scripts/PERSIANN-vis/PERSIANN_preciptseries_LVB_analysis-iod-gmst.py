# -*- coding: utf-8 -*-
"""
Created on Fri May  5 19:38:10 2023

@author: rpietroi

to do:
    1. gmst and precip relationship stratifying on iod 
    2. detrend IOD? Like they did for ENSO in WWA EA drought? Use detrended as covariate in GEV fit 
    3. Detrend and standardize IOD and precip and see correlation / causal effect ? see reading stats course


"""

#%%
import sys, os 
import numpy as np
import statsmodels
import pandas as pd 
import matplotlib.pyplot as plt
import statsmodels
import statsmodels.api as sm

# set working directory
os.chdir(os.path.dirname(__file__)) #os. getcwd() to check

#set data-path
data_path = os.path.join('..', '..', "data") 

#set figures-path
fig_path = os.path.join('..', '..', "figures/figures_10_mar23_precipiod") 
if os.path.exists(fig_path):
    print("The directory", fig_path, "exists!")
else:
    os.makedirs(fig_path)
    print("The directory", fig_path, "was made!")
    

iod_path = os.path.join(data_path, 'data-modified/iod-obs/dmi_seasonal_1870_2021.csv')
#precip_path = os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_seasonal_ppn_intrp_1983_2020.csv')
precip_path = os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_seasonal_ppn_fillclimato_1983_2020.csv')
gmst_path = os.path.join(data_path, 'input-data/gmst-obs/igiss_al_gl_a_4yrlo.csv')

iod_df = pd.read_csv(iod_path, index_col=0)
precip_df = pd.read_csv(precip_path,  index_col=0)
gmst_df = pd.read_csv(gmst_path,  index_col=0)

# cut to same years as precip
iod_df = iod_df.loc[precip_df.index]
gmst_df = gmst_df.loc[precip_df.index]

#%%

# extract variables
iod_ond = iod_df.iod_ond
precip_ond = precip_df.p_ond
precip_jf = precip_df.p_jf
precip_mam = precip_df.p_mam
precip_yr = precip_df.p_year
gmst = gmst_df.gmst
years = precip_df.index

# correlation matrices

corr_gmst_precip = precip_df.apply(lambda s: gmst_df.corrwith(s))
corr_gmst_iod = gmst_df.apply(lambda s: iod_df.corrwith(s))
corr_iod_precip =  precip_df.apply(lambda s: iod_df.corrwith(s))
corr_gmst_years = np.corrcoef(years,gmst)[0,1]

print('corr_gmst_precip', '\n', corr_gmst_precip, '\n\n', 
      'corr_iod_precip', '\n',corr_iod_precip, '\n\n', 
      'corr_gmst_iod', '\n', corr_gmst_iod, '\n\n', 
      'corr_gmst_years',corr_gmst_years)

#%%
#scatterplots

def plot_scatterplot(x,y,ax,txtpx,txtpy):
    ax.scatter(x, y)
    ax.set_xlabel(x.name)
    ax.set_ylabel(y.name)
    ax.text(txtpx,txtpy,'r = {:.3f}'.format(np.corrcoef(x, y)[0,1]), transform=ax.transAxes, c='C0', fontsize=12)


fig, axes = plt.subplots(3,3, sharex=False, figsize=(10,6))
txtpx = 0.4
txtpy = 0.85

ax=axes[0,0]
x = gmst
y = precip_yr 
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[0,1]
x = gmst
y = precip_ond
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[0,2]
x = gmst
y = precip_mam
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[1,0]
x = iod_ond
y = precip_ond
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[1,1]
x = iod_ond
y = precip_yr
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[1,2]
x = iod_df.iod_year
y = precip_yr
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[2,0]
x = gmst
y = iod_ond
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[2,1]
x = gmst
y = iod_df.iod_year
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[2,2]
x = years
y = gmst
plot_scatterplot(x,y,ax,txtpx,txtpy)



fig.tight_layout()

#plt.savefig(os.path.join(fig_path,'DMI_precip_gmst_corr-fillclimato.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'DMI_precip_gmst_corr-fillclimato.png'),dpi=300)

#%%

# plot in time 

fig, axes = plt.subplots(3,1, sharex=False, figsize=(6,8))
txtpx = 0.4
txtpy = 0.85

ax=axes[0]
x = precip_yr.index
y = precip_yr 
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[1]
x = precip_ond.index
y = precip_ond
plot_scatterplot(x,y,ax,txtpx,txtpy)

ax=axes[2]
x = precip_mam.index
y = precip_mam
plot_scatterplot(x,y,ax,txtpx,txtpy)

fig.tight_layout()

# %%

# multiple linear regression precip, iod, gmst/time

# a) trend of precip (start with OND and yearly) and time 

def fit_trend(data):
    # if its a series
    x=range(len((data.index)))
    y=data
    fit=sm.OLS(y, sm.add_constant(x)).fit()    
    return fit

def plot_trend(data,fit,color, ax, str):
    x=range(len((data.index)))
    y=data
    ax.scatter(x, y)
    ax.axline([0, fit.params[0]], slope=fit.params[1], ls='--', label="$P_{}$ = {:.2f} + {:.2f} YR (R$^2$={:.3f})".format({str},fit.params[0], fit.params[1], fit.rsquared), c=color)

def detrend(year, data):
         fit = sm.OLS(data, sm.add_constant(year)).fit()
         return fit.resid.rename(data.name)

def scale(data):
         return (data - data.mean()) / data.std()
    
#%%
fig = plt.figure(figsize=(9,7)) #
ax = plt.subplot()

fit_ond = fit_trend(precip_ond)
plot_trend(precip_ond, fit_ond, 'C0', ax=ax, str='OND')

fit_mam = fit_trend(precip_mam)
plot_trend(precip_mam, fit_mam, 'C1', ax=ax, str='MAM')

plt.legend()
ax.set_ylabel('Precipitation (mm)')
ax.set_xticks(range(0,38,5), range(1983, 1983+len(range(38)),5))

#print(fit_ond.params, '\n', fit_ond.rsquared, '\n')
#print(fit_mam.params, '\n', fit_mam.rsquared, '\n')

#%%


# trend in iod in these years

fig = plt.figure(figsize=(9,7)) #
ax = plt.subplot()

fit= fit_trend(iod_ond)
print(fit.params, '\n', fit.rsquared, '\n')
print(fit.summary())

ax.scatter(iod_ond.index, iod_ond)
ax.axline([iod_ond.index[0], fit.params[0]], slope=fit.params[1], ls='--', label="$IOD_{}$ = {:.2f} + {:.2f} YR (R$^2$={:.3f})".format({'OND'},fit.params[0], fit.params[1], fit.rsquared))

plt.legend()
ax.set_ylabel('IOD')

#plt.savefig(os.path.join(fig_path,'DMI_trend_1983-2020.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'DMI_trend_1983-2020.png'),dpi=300)


# %%

# do this again with standardized data, and detrended IOD 


# b) fit precip ond and iod ond

fit = sm.OLS(precip_ond, sm.add_constant(iod_ond)).fit()
print(fit.params, '\n', fit.rsquared, '\n')
print(fit.summary())

plt.scatter(iod_ond, precip_ond)
plt.xlabel("IOD OND index")
plt.ylabel("Precip OND")
plt.axline([0, fit.params[0]], slope=fit.params[1], color='black', ls='--', 
           label="P_OND = {:.2f} + {:.2f} IOD (R$^2$={:.2f})".format(fit.params[0], fit.params[1], fit.rsquared))
plt.legend()

#plt.savefig(os.path.join(fig_path,'precip-iod-ond-regression.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'precip-iod-ond-regression.png'),dpi=300)


#%% b2) fit precip on detrended iod 

iod_ond_detr = detrend(iod_ond.index, iod_ond)
#iod_ond_detr.plot()
#iod_ond.plot()

fit = sm.OLS(precip_ond, sm.add_constant(iod_ond_detr)).fit()
print(fit.params, '\n', fit.rsquared, '\n')
print(fit.summary())

plt.scatter(iod_ond_detr, precip_ond)
plt.xlabel("IOD OND index (detrended)")
plt.ylabel("Precip OND")
plt.axline([0, fit.params[0]], slope=fit.params[1], color='black', ls='--', 
           label="P_OND = {:.2f} + {:.2f} IOD (R$^2$={:.2f})".format(fit.params[0], fit.params[1], fit.rsquared))
plt.legend()


#plt.savefig(os.path.join(fig_path,'precip-iod-detr-ond-regression.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'precip-iod-detr-ond-regression.png'),dpi=300)


# %%

# c) MLR precip OND on year and IOD OND


data_year = precip_ond.index
fit_mlr = sm.OLS(precip_ond, sm.add_constant(np.stack([data_year, iod_ond], axis=1))).fit()
print(fit_mlr.params, '\n', fit_mlr.rsquared, '\n')
print(fit_mlr.summary())

print(np.corrcoef(data_year, precip_ond)[0,1])

#%%

# c2) MLR precip OND on year and IOD OND detrended
data_year = precip_ond.index
fit_mlr = sm.OLS(precip_ond, sm.add_constant(np.stack([data_year, iod_ond_detr], axis=1))).fit()
print(fit_mlr.params, '\n', fit_mlr.rsquared, '\n')
print(fit_mlr.summary())

print(np.corrcoef(data_year, precip_ond)[0,1])

# about the same, trend with year is not robust... similar regression coefficent to when I regressed precip just on YEAR



# %%

# c) MLR precip OND on gmst and IOD OND


fit_mlr = sm.OLS(precip_ond, sm.add_constant(np.stack([gmst, iod_ond], axis=1))).fit()
print(fit_mlr.params, '\n', fit_mlr.rsquared, '\n')
print(fit_mlr.summary())

print(np.corrcoef(gmst, precip_ond)[0,1])

