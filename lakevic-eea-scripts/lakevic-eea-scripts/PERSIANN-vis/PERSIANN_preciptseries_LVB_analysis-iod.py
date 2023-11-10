# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:11:08 2023

@author: rpietroi

this script: 
    Analysis Precip and IOD to see which IOD months to use
    clean the data of seasonal precip and seasonal iod, incl. fill nan in daily precip 

result: use iod OND, its the most positively correlated with precip_iod, precip_yearly and even precip_mam, and IOD_yearly 

to do elsewhere:
    1. detrend IOD? Like they did for ENSO in WWA EA drought? Use detrended as covariate? 
    2. Detrend and standardize IOD and precip and see correlation / causal effect ? 
    3. Corr with temperature stratified/conditioned on IOD 


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
    

#%%
# iod data

iod_path = os.path.join(data_path, "input-data/IOD/dmi.had.long.data")
iod_raw = np.loadtxt(iod_path,  skiprows=1, max_rows=2021-1870+1) # from jan 1870 to dec 2021
years = iod_raw[:,0]

plt.plot(iod_raw[:,1:13].flatten())
plt.scatter(range(len(iod_raw[:,1:13].flatten())), iod_raw[:,1:13].flatten() )

#%%

# extract seasonal iod 

iod_ond = np.mean(iod_raw[:,10:13], axis=1 )
iod_jf = np.mean(iod_raw[:,1:3], axis=1 )
iod_mam = np.mean(iod_raw[:,3:6], axis=1 )
iod_year = np.mean(iod_raw[:,1:13], axis=1 )

# linear regressions with time = trend 

x=years 
y=iod_ond
fit_ond=sm.OLS(y, sm.add_constant(x)).fit()

y=iod_jf
fit_jf=sm.OLS(y, sm.add_constant(x)).fit()

y=iod_mam
fit_mam=sm.OLS(y, sm.add_constant(x)).fit()

y=iod_year
fit_year=sm.OLS(y, sm.add_constant(x)).fit()


# plot 

alpha = 0.4
fig, ax = plt.subplots(figsize=(9, 5))

plt.scatter(years, iod_ond, label = 'IOD OND', alpha=alpha)
xvals = years[[0,-1]]
yvals = fit_ond.params[0] + xvals * fit_ond.params[1]
plt.plot(xvals, yvals, ls='--', label="IOD = {:.2f} + {:.4f} YR (R$^2$={:.3f})".format(fit_ond.params[0], fit_ond.params[1], fit_ond.rsquared))

plt.scatter(years, iod_jf, label = 'IOD JF', alpha=alpha)
yvals = fit_jf.params[0] + xvals * fit_jf.params[1]
plt.plot(xvals, yvals, ls='--', label="IOD = {:.2f} + {:.4f} YR (R$^2$={:.3f})".format(fit_jf.params[0], fit_jf.params[1], fit_jf.rsquared))

plt.scatter(years, iod_mam, label = 'IOD MAM', alpha=alpha)
fit = fit_mam
yvals = fit.params[0] + xvals * fit.params[1]
plt.plot(xvals, yvals, ls='--', label="IOD = {:.2f} + {:.4f} YR (R$^2$={:.3f})".format(fit.params[0], fit.params[1], fit.rsquared))

plt.scatter(years, iod_year, label = 'IOD year', alpha=alpha)
fit = fit_year
yvals = fit.params[0] + xvals * fit.params[1]
plt.plot(xvals, yvals, ls='--', label="IOD = {:.2f} + {:.4f} YR (R$^2$={:.3f})".format(fit.params[0], fit.params[1], fit.rsquared))


plt.ylabel('Dipole Mode Index')
plt.legend()

# save
#plt.savefig(os.path.join(fig_path,'DMI_seasonal_trend.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'DMI_seasonal_trend.png'),dpi=300)


#%% Save the raw and detrended IODs, and calculate correlation between seasonal and yearly

d = {'years':years.astype(int), 'iod_ond': iod_ond, 'iod_jf':iod_jf, 'iod_mam': iod_mam, 'iod_year':iod_year}
iod_df = pd.DataFrame(data=d).set_index('years')

corr_matrix_iod = iod_df.corr()

# save
#iod_df.to_csv(os.path.join(data_path, 'data-modified/iod-obs/dmi_seasonal_1870_2021.csv'), index=True, float_format='%.5f')

# to do: detrend and save 

#%%

# read in ppn data # from 1983 to 2020
precip_path = os.path.join(data_path, "data-modified/PERSIANN")
data_monthly = pd.read_csv(os.path.join(precip_path, 'PERSIANN_basinlake_monthly_ppn_1983_2020.csv'))
data_daily = pd.read_csv(os.path.join(precip_path, 'PERSIANN_basinlake_daily_ppn_1983_2020.csv'))

#monthly data
data_monthly['time'] = pd.to_datetime(data_monthly['time'])
data_monthly = data_monthly.set_index(['time'])

#daily data
data_daily['time'] = pd.to_datetime(data_daily['time'])
data_daily = data_daily.set_index(['time'])
data_daily['doy'] = data_daily.index.dayofyear

# get daily climatology
daily_climato = data_daily.groupby(data_daily.index.dayofyear).mean('time')

# deal with nan values in data_daily to make own resampled data_monthly and check difference

# find all nans
badrows = np.isnan(data_daily['precipitation'])
print(badrows.sum())

# replace nans with climatology 
data = data_daily.copy()
data.loc[badrows,'precipitation'] = daily_climato.loc[data.loc[badrows,'doy'].values, 'precipitation'].values
data_daily_fill = data.drop(columns=['doy'])
data_daily_fill['precipitation'].plot()

# check how many nans are left 
badrows = np.isnan(data_daily['precipitation'])
print(badrows.sum())

# resample to monthly
data_monthly_rsp = data_daily_fill.resample('M', label='right').sum()
data_yearly_rsp = data_daily_fill.resample('Y', label='right').sum()


# alternative - linearly interpolate - NO ! Use the climatology one
# this makes quite a big difference in the early years (1980s)
data = data_daily.copy()
data.loc[badrows,'precipitation'] = np.interp(np.flatnonzero(badrows), np.flatnonzero(~badrows), data.loc[~badrows, 'precipitation'])
data_daily_intr = data.drop(columns=['doy'])
data_daily_intr['precipitation'].plot()

data_monthly_rsp_intr = data_daily_intr.resample('M', label='right').sum()
data_yearly_rsp_intr = data_daily_intr.resample('Y', label='right').sum()

#%%

fig, ax = plt.subplots()
data_monthly_rsp['precipitation'].plot(ax=ax, label='filled climatology')
data_monthly_rsp_intr['precipitation'].plot(ax=ax, label='filled interpolate')
data_monthly['precipitation'].plot(ax=ax, label='raw')
plt.legend()

fig, ax = plt.subplots()
data_monthly_rsp['precipitation'].resample('Y', label='right').sum().plot(ax=ax, label='filled climatology')
data_yearly_rsp_intr['precipitation'].plot(ax=ax, label='filled interpolate')
data_monthly['precipitation'].resample('Y', label='right').sum().plot(ax=ax, label='raw')
plt.legend()

#%% 

# Make plot nicer 

fig, ax = plt.subplots()
data_monthly_rsp['precipitation'].resample('Y', label='right').sum().plot(ax=ax, label='PERSIANN-CDR after filling missing days with climatology')
data_monthly['precipitation'].resample('Y', label='right').sum().plot(ax=ax, label='PERSIANN-CDR before filling missing days')
plt.legend()
plt.xlabel('')
plt.ylabel("Precipitation (mm yr$^{-1}$)")

# save
#plt.savefig(os.path.join(fig_path,'PERSIANN_basinlake_fillnan_climato_SI_v2.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'PERSIANN_basinlake_fillnan_climato_SI_v2.png'),dpi=300)


#%%

# save the data 
#data_daily_intr.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_daily_ppn_intrp_1983_2020.csv'), index=True, float_format='%.5f')
#data_monthly_rsp_intr.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_monthly_ppn_intrp_1983_2020.csv'), index=True, float_format='%.5f')
#data_yearly_rsp_intr.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_yearly_ppn_intrp_1983_2020.csv'), index=True, float_format='%.5f')


# save the data - error was here !!! i saved interpolated yearly data instead of climatology filled
#data_daily_fill.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_daily_ppn_climato_fill_1983_2020.csv'), index=True, float_format='%.5f')
#data_monthly_rsp.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_monthly_ppn_climato_fill_1983_2020.csv'), index=True, float_format='%.5f')
#WRONG ::: data_yearly_rsp_intr.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_yearly_ppn_climato_fill_1983_2020.csv'), index=True, float_format='%.5f')
#data_yearly_rsp.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_yearly_ppn_climato_fill_1983_2020_correct.csv'), index=True, float_format='%.5f')






#%% calc seasonal precips 

#precip_monthly_df = data_monthly_rsp_intr # interpolated - old
#precip_monthly_df = data_monthly_rsp # climatology - new
precip_monthly_df = data_monthly["1985-01-31": "2020-12-31"] # not corrected, now cut it to 1985

# OND
precip_monthly_OND = precip_monthly_df[precip_monthly_df.index.month > 9]
precip_yearly_OND = precip_monthly_OND.resample('Y', label='right').sum().rename(columns={'precipitation':'p_ond'})
# JF
precip_monthly_JF = precip_monthly_df[precip_monthly_df.index.month < 3]
precip_yearly_JF = precip_monthly_JF.resample('Y', label='right').sum().rename(columns={'precipitation':'p_jf'})
#MAM
precip_monthly_MAM = precip_monthly_df[(precip_monthly_df.index.month >2) & (precip_monthly_df.index.month < 6)]
precip_yearly_MAM = precip_monthly_MAM.resample('Y', label='right').sum().rename(columns={'precipitation':'p_mam'})
#JJAS
precip_monthly_JJAS = precip_monthly_df[(precip_monthly_df.index.month >5) & (precip_monthly_df.index.month < 10)]
precip_yearly_JJAS = precip_monthly_JJAS.resample('Y', label='right').sum().rename(columns={'precipitation':'p_jjas'})

#%%
precip_df = precip_yearly_OND.merge(precip_yearly_JF,left_index=True, right_index=True).merge(precip_yearly_MAM,left_index=True, right_index=True).merge(precip_yearly_JJAS,left_index=True, right_index=True)
# precip_df = precip_df.merge(data_yearly_rsp_intr,left_index=True, right_index=True).rename(columns={'precipitation':'p_year'})
precip_df = precip_df.merge(data_yearly_rsp,left_index=True, right_index=True).rename(columns={'precipitation':'p_year'})

precip_df.index = pd.to_datetime(precip_df.index).strftime('%Y').rename('years')

# save the data 
#precip_df.to_csv(os.path.join(data_path, 'data-modified/PERSIANN/PERSIANN_basinlake_seasonal_ppn_fillclimato_1983_2020.csv'), index=True, float_format='%.5f')

#%%

"""

Analysis Precip and IOD to see which IOD months to use 

result: use OND, its the most positively correlated with precip_iod, precip_yearly and even precip_mam, and IOD_yearly 

to do: 
    1. detrend IOD? Like they did for ENSO? Use detrended as covariate? See correlation with precip
    2. Detrend and standardize IOD and precip and see correlation / causal effect ? 
    3. Corr with temperature stratified/conditioned on IOD 

"""
iod = iod_df.loc[1985:2020, :]
iod.index = pd.to_datetime(iod.index, format='%Y').strftime('%Y') #.rename('years')
precip = precip_df

#%%

# calculate correlations of seasonal precip and seasonal iod

corr_matrix_iod = iod.corr()
corr_matrix = iod.apply(lambda s: precip.corrwith(s))

#%%

filepath = r"C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\output\precip-iod-corr"

#corr_matrix_iod.to_csv(os.path.join(filepath, 'iod-corr.txt'), sep='\t',float_format='%.5f')
#corr_matrix.to_csv(os.path.join(filepath, 'precip-iod-corr-fillnanclimato.txt'), sep='\t',float_format='%.5f')

# go for IOD OND 

#%%

fig, axes = plt.subplots(3,1, figsize=(6,8))

ax=axes[0]
ax.set_title('(a)', loc='left', size=13)
ax.scatter(iod['iod_ond'], precip['p_ond'], label='r = {:.3f}'.format(corr_matrix.iloc[0,0]))
ax.set_xlabel('Dipole Mode Index (OND)')
ax.set_ylabel('Precipitation$_{OND}$ (mm)')
ax.text(.82,.6,'r = {:.3f}'.format(corr_matrix.iloc[0,0]), transform=ax.transAxes, c='C0', fontsize=12)
#ax.legend()

ax=axes[1]
ax.set_title('(b)', loc='left', size=13)
ax.scatter(iod['iod_ond'], precip['p_year'], label='r = {:.3f}'.format(corr_matrix.iloc[0,3]))
ax.set_xlabel('Dipole Mode Index (OND)')
ax.set_ylabel('Precipitation$_{year}$ (mm)')
ax.text(.82,.6,'r = {:.3f}'.format(corr_matrix.iloc[0,3]), transform=ax.transAxes, c='C0', fontsize=12)
#ax.legend()

ax=axes[2]
ax.set_title('(c)', loc='left', size=13)
ax.scatter(iod['iod_year'], precip['p_year'], label='r = {:.3f}'.format(corr_matrix.iloc[4,3]))
ax.set_xlabel('Dipole Mode Index (yearly)')
ax.set_ylabel('Precipitation$_{year}$ (mm)')
ax.text(.82,.6,'r = {:.3f}'.format(corr_matrix.iloc[4,3]), transform=ax.transAxes, c='C0', fontsize=12)
#ax.legend()

fig.tight_layout()

# save
#plt.savefig(os.path.join(fig_path,'DMI_precip_corr_nan_cut1985.pdf'),dpi=300)
#plt.savefig(os.path.join(fig_path,'DMI_precip_corr_nan_cut1985.png'),dpi=300)


# %%
