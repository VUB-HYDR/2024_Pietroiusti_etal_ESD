# Functions from PERSIANN-nan-figures-1985 ipynb
# July 2023 rosa.pietroiusti@vub.be
# ---------------------------------------------------------------

import numpy as np
import pandas as pd
import os, glob 
import math
import xarray as xr
import geopandas as gpd
import regionmask as regionmask
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4
import statsmodels
import statsmodels.api as sm
import matplotlib.pyplot as plt


# general 
# ---------------------------------------------------------------

def open_ncfile(ncfile, var):

    with xr.open_dataset(ncfile, engine='netcdf4') as ds:
        print(str(list(ds.keys())))
        print(str(list(ds.dims)))
        da = ds[str(var)] 
        
        return da

    
def summary(da):
    print('min', da.min(), '\nmax', da.max(), '\nmean', da.mean())


def mask_da(da, shp_path):
    # mask data with a shapefile 
    
    shp = gpd.read_file(shp_path, crs="epsg:4326")
    mask = regionmask.mask_geopandas(shp, np.array(da.lon), np.array(da.lat)) + 1
    mask_da = da.where(mask == 1 )
    return mask_da

def cut_full_years(da):
    
    if da.time[-1].dt.dayofyear < 365:
        i_end = -2
    else:
        i_end = -1
    if da.time[0].dt.dayofyear > 1:
        i_start = 1
    else:
        i_start = 0 
        
        da_out = da.sel(time=slice("{}-01-01".format(np.unique(da.time.dt.year)[i_start]), "{}-12-31".format(np.unique(da.time.dt.year)[i_end])))
    
    print(i_start, i_end)
    return da_out

def find_nan(da):
    # find missing days in a dataframe and output the indices 
    
    indices = []
    data = da
    for t in data.time:
        if np.isnan(data.sel(time=t).values).all() == True:
            indices.append(t.values)
        
    return pd.DatetimeIndex(indices) 

def count_daysperyear(dlist):
    # in a list of pd.datetimeindex find how many days are in each year (e.g. list of nan days)
    
    years = np.unique(dlist.year)
    df_out = pd.DataFrame({'year': years, 'ndays': np.zeros(len(years))}).set_index('year')

    for year in years:
        slice_days = dlist[dlist.year == year]
        ndays = len(slice_days)
        df_out.loc[year] = ndays
    
    return df_out
    
def count_dayspermonth(dlist):
    # in a list of pd.datetimeindex find how many days are in each month 
    
    months = range(1,13)
    df_out = pd.DataFrame({'month': months, 'ndays': np.zeros(len(months))}).set_index('month')
    
    for i in months:
        slice_days = dlist[dlist.month == i]
        ndays = len(slice_days)
        df_out.loc[i] = ndays
        
    return df_out

# fill missing data
# ---------------------------------------------------------------

def fill_climato_spatial(da):
    # fill spatial dataarray with climatology (original script: PERSIANN_fillnan-spatial.py)
    
    # get missing days
    indices = []
    data = da
    for t in data.time:
        if np.isnan(data.sel(time=t).values).all() == True:
            indices.append(t.values)
    print('nan days=', len(indices))
    
    # calc climatology
    climatology_daily = data.groupby("time.dayofyear").mean('time')
    
    # fill nan days
    da_fill = data.copy()
    for i in indices:
        da_fill.loc[dict(time=i)] = climatology_daily[da_fill.sel(time=i).time.dt.dayofyear.values - 1].values
    
    # check if any missing days are left
    check = []
    data = da_fill
    for t in data.time:
        if np.isnan(data.sel(time=t).values).all() == True:
            check.append(t.values)
    print('nan days left =', len(check))
    
    return da_fill



def fill_climato_timeseries(df_daily, varname):
    # fill timeseries dataframe (original script: PERSIANN_preciptseries_LVB_analysis-iod.oy)
    
    #calc climatology
    daily_climato = df_daily.groupby(df_daily.index.dayofyear).mean('time')
    
    # find all nans
    badrows = np.isnan(df_daily[varname])
    print('nan days = ', badrows.sum())

    # replace nans with climatology 
    data = df_daily.copy()
    data.loc[badrows,varname] = daily_climato.loc[data.loc[badrows,'doy'].values, varname].values
    data_daily_fill = data.drop(columns=['doy'])
    #data_daily_fill[varname].plot()

    # check how many nans are left 
    badrows = np.isnan(data_daily[varname])
    print('nan days left =', badrows.sum())
    
    return data_daily_fill
    

    
    
# statistics 
# ---------------------------------------------------------------
def fit_trend(data, var):
    x=range(len((data.index)))
    y=data[var]
    fit=sm.OLS(y, sm.add_constant(x)).fit()    
    return fit, x, y

def f_test(fit):
    print('F test - F value:', fit.fvalue, 'p-value', fit.f_pvalue) # why not mann kendal or walden test?
    


# Extract seasons
# ---------------------------------------------------------------
def extract_seasonal_precip_basin(precip_data, basin_shp):
    #Monthly
    precip_monthly_df = mask_da(precip_data, basin_shp).mean(dim=['lat','lon']).resample(time='M', label='right').sum().to_dataframe()
    # Yearly 
    precip_yearly_df = mask_da(precip_data, basin_shp).mean(dim=['lon','lat']).resample(time='Y', label='right').sum().to_dataframe()
    # OND
    precip_yearly_OND = precip_monthly_df[precip_monthly_df.index.month > 9].resample('Y', label='right').sum()
    # JF
    precip_yearly_JF = precip_monthly_df[precip_monthly_df.index.month < 3].resample('Y', label='right').sum()
    #MAM
    precip_yearly_MAM = precip_monthly_df[(precip_monthly_df.index.month >2) & (precip_monthly_df.index.month < 6)].resample('Y', label='right').sum()
    #JJAS
    precip_yearly_JJAS = precip_monthly_df[(precip_monthly_df.index.month >5) & (precip_monthly_df.index.month < 10)].resample('Y', label='right').sum()
    
    precip_seasons = pd.DataFrame({'year': precip_yearly_df[precip_yearly_df.columns[0]] ,
                                  'JF': precip_yearly_JF[precip_yearly_JF.columns[0]],
                                  'MAM': precip_yearly_MAM[precip_yearly_MAM.columns[0]],
                                  'JJAS': precip_yearly_JJAS[precip_yearly_JJAS.columns[0]],
                                  'OND': precip_yearly_OND[precip_yearly_OND.columns[0]]})
    
    return precip_seasons #[precip_yearly_df, precip_yearly_OND, precip_yearly_JF, precip_yearly_MAM, precip_yearly_JJAS]

# Plot trends in seasonal precipitation
# ----------------------------------------

def plot_trends_seasonal_precip(precip_data_seasons, varname, start_year, end_year):
    #set up figure
    rs = 2
    cs = 1
    fig = plt.figure(figsize=(9,7)) 
    ax = plt.subplot(rs, cs, 1)
    
    # yearly precip
    plt.scatter(range(len((precip_data_seasons.index))), precip_data_seasons.year)

    # fit trend 
    data_fit = precip_data_seasons.loc["{}-01-01".format(start_year):"{}-12-31".format(end_year)]
    print('trends fitted between {} - {}'.format(data_fit.index.year[0], data_fit.index.year[-1])   )

    shift = len(precip_data_seasons.loc[:"{}-01-01".format(start_year)])
    fit, x, y = fit_trend(data_fit, 'year')

    #plot trendlint
    plt.plot([x[0] + shift, x[-1] + shift], 
             [x[0]*fit.params[1]+fit.params[0] , x[-1]*fit.params[1]+fit.params[0]],
            ls='--', label="$P_{}$ = {:.2f} + {:.2f} YR (R$^2$={:.3f})".format('{YR}',fit.params[0], fit.params[1], fit.rsquared))
    
    # setup situation
    ax.set_title('(a)', loc='left')
    ax.legend(loc=(1.01,.88),frameon=False)
    ax.set_ylabel('Precipitation (mm)')
    ax.set_xticks(range(0,len((precip_data_seasons.index)),5), range(precip_data_seasons.index.year[0], precip_data_seasons.index.year[0]+len((precip_data_seasons.index)),5))
    
    # seasonal precip
    ax = plt.subplot(rs, cs, 2)
    

    for period in ['OND', 'JF', 'MAM', 'JJAS']:
        fit, x, y = fit_trend(data_fit, period)
        ax.scatter(range(len((precip_data_seasons.index))), precip_data_seasons[period])
            #plot trendlint
        ax.plot([x[0] + shift, x[-1] + shift], 
                 [x[0]*fit.params[1]+fit.params[0] , x[-1]*fit.params[1]+fit.params[0]],
                ls='--', label="$P_{}$ = {:.2f} + {:.2f} YR (R$^2$={:.3f})".format(period,fit.params[0], fit.params[1], fit.rsquared))
    ax.set_title('(b)', loc='left')
    ax.legend(loc=(1.01,.58),frameon=False)
    ax.set_ylabel('Precipitation (mm)')
    ax.set_xticks(range(0,len((precip_data_seasons.index)),5), range(precip_data_seasons.index.year[0], precip_data_seasons.index.year[0]+len((precip_data_seasons.index)),5))





    
# Plotting functions (Figure 2 paper) 
# ---------------------------------------------------------------
def plot_map(ax, data, title, cmap, levels, label, norm=False):
    gl = ax.gridlines(draw_labels=True,x_inline=False,y_inline=False, linewidth=0.8, color='gray', alpha=0.6, linestyle='dotted')
    gl.ylabels_right = False
    gl.xlabels_top = False
    if norm==True:
        plot = ax.contourf(lons, lats, data,  cmap=cmap, norm=mcolors.CenteredNorm(), levels=levels,  transform=ccrs.PlateCarree())
    if norm==False:
        plot = ax.contourf(lons, lats, data,  cmap=cmap, levels=levels,  transform=ccrs.PlateCarree())
    ax.set_title(title, loc = 'left', size=12 )
    ax.coastlines(color='dimgray', linewidth=0.5)
    ax.add_feature(cfeature.LAKES, facecolor="none", edgecolor='black', linewidth=0.8 )
    return plot, cmap

def plot_basin(ax, col):
    basin_shp.boundary.plot(ax=ax, edgecolor=col, linewidth=1.5, label='basin')
    

def make_month_xaxis_labels_ticks():
    startdate = np.datetime64('2000-01-01')
    enddate = np.datetime64('2000-12-31')
    alldays = pd.date_range(start=startdate,end=enddate)

    list_ticks = np.zeros(12)
    for i in range(1,13):
        sel = alldays[alldays.month == i].dayofyear[0] -1
        list_ticks[i-1] = sel

    list_monthnames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
    return list_ticks, list_monthnames

