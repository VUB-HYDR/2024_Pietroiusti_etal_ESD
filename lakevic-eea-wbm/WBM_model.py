#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------
Lake Victoria Water Balance Model 
-------------------------------------
    - Read in precipitation data and divide into Plake and Pbasin --> 1 per day for each year
    - Convert to 1D P_lake
    - Read in evaporation data, reproject to own grid --> 1 per day for 1 year, copied for each year
    - Convert to 1D E_lake
    - Read in outflow data --> 1 per day for each year 
    - Convert to 1D Qout
    - Calculate runoff per grid cell using Curve Number method (see paper for references)  
    - Convert to 1D Qin 
    - Calculate WB

Created on Mon Apr 25 15:28:27 2022
Update 20 June 2023

@author: rosa.pietroiusti@vub.be

Translated and adapted from Vanderkelen et al. 2018 MATLAB model

"""

#=============#
#== IMPORTS ==#
#=============#

# basics
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys, os , glob, gc
import pandas as pd
import rioxarray
from osgeo import gdal

# For vector data
import regionmask

# keep this - used it to reproject evap
# python-cdo  
#from cdo import *
#cdo=Cdo()

#=========================================================#
#== INITIALISE FILEPATHS, PHYS CONSTANTS AND GEOMETRY ====#
#=========================================================#

# append path to call other scripts from here
WBM_path = os.path.dirname(__file__) # /$VSC_DATA_VO_USER/vsc10419/Thesis/WBM-git/lakevic-eea-wbm
sys.path.append(WBM_path)
from WBM_settings import * # imports all variables - put all user defined settings and paths and flags here (?) like a namelist 
import WBM_inicon as inicon # initialise physical constants - e.g. print(inicon.Lvap)
import WBM_inigeom as inigeom # lake victoria boundaries, size of a grid cell and resolution, shapefiles - e.g. print(inigeom.res_m, inigeom.A_cell)
from WBM_myfunctions import compare_lake_areas, check_zeroarray, check_valuesarray, check_dimorders, check_punits, fig_showsave, update_fign # functions I defined 


# # For raster data - check this, delete if not necessary
# if run_where == 0 :
#     from osgeo import gdal # local computer
# if run_where == 1 :
#     import gdal # for HPC 


# Start message from settings file
print(start_message_1)
print(start_message_2)

#%% read in geo data / read in WB terms data 
print('...reading in data')

#=========================#
#== READ WBM TERMS DATA ==#
#=========================#

# Precipitation 
#=========================#
if run_type == 0:
    var = 'precipitation'
if run_type == 1 or run_type == 2:
    var = 'pr' #SOFT CODE THIS
with xr.open_dataset(filepath_precip, decode_coords="all") as ds:
    print('precipitation data')
    print(str(list(ds.keys())))
    print(str(list(ds.dims)))
    ds_dims_precip = list(ds.dims) # new: time, lat, lon (reordered, WBM is made to work for this order)   # old: time,lon, lat
    ds_varnames = list(ds.keys()) # get names of variables precipitation
    time_precip = np.array(ds.get('time'))
    lats_precip = np.array(ds.get('lat'))
    lons_precip = np.array(ds.get('lon'))
    precip_raw = ds[str(var)].sel(time=slice(startDATE, endDATE))  # slice to modelling timeperiod
    precip_units = precip_raw.units
    precip_raw.rio.write_crs("epsg:4326", inplace=True) # \\\\ see if this should be done here or not and that CRS is correct !!\\\\
    check_dimorders(ds)
    precip_raw = check_punits(precip_raw) # this function modifies ppn, check it is ok


# Evaporation (Latent Heat Flux)
#=========================#
# Cut and remapped to own grid -  using cdo-python, different remapping algorithms are possible: see https://stackoverflow.com/questions/61806343/regrid-netcdf-file-in-python
    #remapbic : bicubic interpolation
    #remapbil : bilinear interpolation
    #remapnn : nearest neighbour interpolation
    #remapcon : first order conservative remapping
    #remapcon2 : 2nd order conservative remapping
#=========================#

# cdo.remapbil(filepath_grid, input=filepath_evap, output=filepath_evap_remap) # the remapping done with bilinear interpolation 

var = 'ALHFL_S' 
with xr.open_dataset(filepath_evap_remap, decode_coords="all") as ds: # open the remapped file
    print('evaporation COSMO-CLM data')
    print(str(list(ds.keys())))
    print(str(list(ds.dims)))
    ds_dims_evap = list(ds.dims) # time,lon, lat ---> check this is OK 
    ds_varnames = list(ds.keys()) # get names of variables
    time_evap = np.array(ds.get('time'))
    lats_evap = np.array(ds.get('lat'))
    lons_evap = np.array(ds.get('lon'))
    evap_raw = ds[str(var)] # order: time, y, x
    evap_raw.rio.write_crs("epsg:4326", inplace=True) # \\\\ see if this should be done here or not and that CRS is correct \\\\
    check_dimorders(ds)
    
# Outflow multi-source (m^3 /day)
#=========================#
outflow_raw = pd.read_csv(filepath_outflow)
outflow = outflow_raw.copy()
outflow.columns = ['date', 'outflow']
outflow['date'] = pd.to_datetime(outflow['date'])
outflow = outflow.set_index(['date'])

# Observed lake levels 
#=========================#
df = pd.read_csv(filepath_lakelevels_hist)
df['date'] = pd.to_datetime(df['date'])
lakelevels = df.set_index(['date'])

#===================#
#== READ GEO DATA ==# 
#===================#

lake_shp = inigeom.lake_shp
basin_shp = inigeom.basin_shp 

# Make lake and basin masks - I moved this here because it was giving problems in the inigeom file, see if better to truncate lat/lon values in Ppn xarray/get higher tolerance level for masking so that I can move this to inigeom script
#=========================#
rmask_lake = regionmask.mask_geopandas(lake_shp, lons_precip, lats_precip) + 1       # only lake set = 1
rmask_basinlake = regionmask.mask_geopandas(basin_shp, lons_precip, lats_precip) + 1 # lake+ basin set = 1
rmask_basin = rmask_basinlake * np.where(np.isnan(rmask_lake),1,np.nan) # only basin set = 1


#%% Precipitation
print('...manipulating precipitation')

#=====================#
#== 1) MANIP PRECIP ==#
#=====================#
# Convert units and separate precip over lake and precip over basin over modelling time-period 

if flag_plotprecip == 1:
    # --- Test --- 
    precip_raw.mean('time').plot()
    plt.title('P_mean mm/day')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# Select to modelling period and convert from mm/day (kg/m2/day) to m/day (1000 kg/m2/day)
precip = precip_raw * 1e-3

# Clear some memory
del precip_raw 
gc.collect()

if flag_plotprecip == 1:
    # --- Test --- 
    precip.mean('time').plot()
    plt.title('P_mean m/day')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# Check the dates are all there and there are no negative numbers 
if list(set(np.array(DATEs)) - set(np.array(precip['time']))): 
    print("model dates and precip (sliced) dates do not match")
if np.any(precip < 0): 
    print("there are some negative values in your precipitation input file \n ... setting them to 0")
    precip = precip.where(precip>0, 0) #0 or neg are turned to 0
if np.any(precip != precip): 
    print("there are some Nan values in your precipitation input file \n ... go check them")

# Mask : Precipitation over the lake (using regionmask)     
precip_lake = precip.where(rmask_lake == 1 )

# Mask: Precipitation over the basin (regionmask)
precip_basinlake = precip.where(rmask_basinlake == 1 )
precip_basin = precip.where(rmask_basin == 1 )

if flag_plotprecip == 1:
    # --- Test --- 
    rmask_lake.plot()
    plt.title('lake mask (regionmask)')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    
    # --- Test --- 
    rmask_basinlake.plot()
    plt.title('basinlake mask (regionmask)')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    rmask_basin.plot()
    plt.title('basin mask (regionmask)')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    
    # --- Test --- 
    precip_lake[0].plot() # .transpose() not necessary if dimensions are reordered to time,lat,lon
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    precip_lake.mean('time').plot() 
    plt.title('mean precip lake')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    precip_basin[6].plot()
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    precip_basin.mean('time').plot()
    plt.title('mean precip basin')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    precip_basinlake[0].plot()
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    precip_basinlake.mean('time').plot()
    plt.title('mean precip basinlake')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# Check masking is ok (make this more foolproof)
print('number of lake pixels {}'.format((~np.isnan(rmask_lake)).sum().values))
print('number of basin+lake pixels {}'.format( (~np.isnan(rmask_basinlake)).sum().values)  )
print('number of basin pixels {}'.format( (~np.isnan(rmask_basin)).sum().values)  )
if (~np.isnan(precip_lake[0])).sum().values != 1298: 
    print("Your masked P_lake might not have the right number of lake pixels on all days")
if (~np.isnan(precip_basinlake[0])).sum().values != 4934: 
    print("Your masked basinlake might not have the right number of pixels on all days")
if (~np.isnan(precip_basin[0])).sum().values != 3637: 
    print("Your masked P_basin might not have the right number of pixels on all days")

  

# /////////// Can insert modifications to P_lake or P_basin here //////////// 

#=================#
#== 1D: P_lake ===# # xarray dataarray
#=================#

# Get 1D timeseries : P_lake

#P_wb = precip_lake.mean('lat').mean('lon') # xarray dataarray ---> Incorrect way of taking mean, overestimates P_lake 
P_wb = precip_lake.mean(dim=('lat','lon'), skipna=True) # ---> I think this is more correct way of taking spatial mean, mean P_lake becomes about 126 mm/mo (for period 1993-2014) which is approx what Inne had, skipna is True by default so not really necessary to specify
# --- Test --- 
P_wb.plot()
plt.title('P_wb [m/d]')
fig_showsave(plt, flag_savefig, fig_path, fig_n)
fig_n = update_fign(flag_savefig, fig_n)

# to get time: np.array(P_wb['time'])
# to get variables: np.array(P_wb) or P_wb.values

# Clear memory
del precip_basinlake , precip_basin , precip_lake
gc.collect()  

#%% Evaporation
print('...manipulating evaporation')
#=====================#
#== 2) MANIP EVAP ====# 
#=====================#
# Convert from W m-2 to mm/day and from negative to positive 
# Repeat climatology to fill modelling time-period

if flag_plotevap == 1:
    # --- Test --- 
    #evap_raw[0].plot() # 366 days (2018 was leap year, repeat the daily value for each E term in WBM based on calendar day!)
    evap_raw.mean('time').plot() 
    plt.title('evap raw mean W m-2')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# Convert from W m-2 to mm/day (kg m-2 day-1) and from neg to pos
evap_conv = - evap_raw / inicon.Lvap * inicon.sec_per_day

if flag_plotevap == 1:
    # --- Test --- 
    evap_conv.mean('time').plot()
    plt.title('evap raw mean mm/day (pos)')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# Convert from mm/day to m/day
evap_conv = evap_conv * 1e-3

if flag_plotevap == 1:
    # --- Test --- 
    evap_conv.mean('time').plot()
    plt.title('evap raw mean m/day')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    evap_conv[0].plot() # 366 days (2018 was leap year, repeat the daily value for each E term in WBM)
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# Mask: evaporation over lake 
evap_lake = evap_conv * (rmask_lake.values) # \\ CHECK THIS CLIPPING !! \\\ IT WORKS BUT WHY??

if flag_plotevap == 1:
    # --- Test --- 
    evap_lake[0].plot()
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

print('number of lake evap pixels {}'.format( (~np.isnan(evap_lake[0])).sum().values)  )
if (~np.isnan(evap_lake[0])).sum().values != 1298: 
    print("Your masked E_lake might not have the right number of lake pixels on all days")
  
#============#
#== 1D: E ===#   # numpy array 
#============#

# Get 1D timeseries : turn 1D for one year 

#E_mean_year = evap_lake.mean('lat').mean('lon') # --> to check agreement with MATLAB model, overestimates a bit
E_mean_year = evap_lake.mean(dim=('lat','lon'), skipna=True) # --> more correct method
# --- Test --- 
E_mean_year.plot()
plt.title('E_wb')
fig_showsave(plt, flag_savefig, fig_path, fig_n)
fig_n = update_fign(flag_savefig, fig_n)

# Numpy method of getting E_lake for all years
# Slicing and hard-coding index of Feb 29th (couldn't find a way of soft-coding that day's index)
# See alternative method using for loop and xarray dates in extrascripts
j = 0 
E_lake = []
for i in range(len(np.unique(DATEs.year))):
    year = np.unique(DATEs.year)[i]
    year_length = sum(DATEs.year == year)
    start_slice = j
    j = j + year_length
    end_slice = j-1
    if year_length == 366: 
        #df[start_slice:end_slice, 'E_lake'] = E_mean_year.values
        E_slice = E_mean_year.values
    elif year_length == 365:
        #df[start_slice:end_slice, 'E_lake'] = E_mean_year.drop_sel(time="28-02-2008", axis=time) # doesn't work ! 
        #df[start_slice:end_slice, 'E_lake'] = E_mean_year.values[np.r_[0:59,60:366]] # get rid of Feb 29th, index 59, doenst work   
        E_slice = E_mean_year.values[np.r_[0:59,60:366]]
    else : 
        print('**ERROR** in length of years ')
    E_lake.extend(E_slice.flatten())
   
# \\\\\ Add a check that length of E is the same as modelling timeperiod \\\\

# Get 1D timeseries as numpy array 

E_wb = np.array(E_lake)
# --- Test --- 
plt.plot(E_wb)
plt.title('E_wb')
fig_showsave(plt, flag_savefig, fig_path, fig_n)
fig_n = update_fign(flag_savefig, fig_n)


# Clear memory 
del evap_raw , evap_lake , evap_conv
gc.collect()


#%% Inflow CN

#=====================#
#== 3) CALC INFLOW ===# 
#=====================#
# Determine CN pixel values of grid based on land cover and hydrologic soil class 
# USDA curve number method: Chapter 10, https://www.nrcs.usda.gov/wps/portal/nrcs/detailfull/national/water/manage/hydrology/?cid=stelprdb1043063 
# Natural Resources Conservation Service (NRCS) method of estimating direct runoff from storm rainfall


# Read map of soil types
dataset = gdal.Open(filepath_soilclass_map, gdal.GA_ReadOnly)
soilclass_map = band = dataset.GetRasterBand(1).ReadAsArray().astype('int64') # numpy array

# Read hydrologic soil classes 
soilclass_hydro = np.loadtxt(filepath_soilclass_txt)

# Read map of land cover
dataset = gdal.Open(filepath_landcover_map, gdal.GA_ReadOnly)
landcover_map = dataset.GetRasterBand(1).ReadAsArray().astype('int64') # numpy array

# Read possible CN values --> cite source 
CN_values = np.loadtxt(filepath_CN_values_txt)  

# /////////// CAN TUNE CN VALUES HERE, MAKE EXTERNAL LOOP //////////// 

#==================#
#== MAKE CN MAP ===#  
#==================#
print('...making curve number map')
# Determine CN value per pixel based on soil hydrological class and land use (under standard moisture conditions, CN-II)

CN_map = np.zeros((len(precip['lat']), len(precip['lon'])), dtype=np.int64) # initialize array
soil_type = np.zeros((len(precip['lat']), len(precip['lon'])), dtype=np.int64) # initialize array
for i in range(len(precip['lat'])):
    for j in range(len(precip['lon'])):
        # determine hydrological soil type per pixel
        soil_type[i,j] = soilclass_hydro[soilclass_map[i,j] - 1] # check the minus one! 
        
        # if water body CN = 100 so then S (max soil retention) = 0, all of it becomes runoff
        if landcover_map[i,j] == 8 or soil_type[i,j] == 0:
            CN_map[i,j] = 100 # in old code set to 100 but was 0 in calculations I guess? Or should it be 100? 
        else: 
            #find CN value. this also depends on antecedent moisture condition (see settings file), default AMC is 2 (where?? in paper AMC=5)
            CN_map[i,j] = CN_values[landcover_map[i,j] - 1, soil_type[i,j] - 1] # check this ! I reduced index by 1 because of Python v. MATLAB indexing

if flag_plotinflow == 1:
    # --- Test --- 
    cmap = plt.get_cmap('viridis', 10)   # discrete colors
    cax = plt.imshow(CN_map, cmap = cmap)
    plt.colorbar(cax, orientation='vertical')
    plt.title('CN map, max {:.0f}, mean {:.0f} \n mean without 100 : {:.0f}'.format(np.nanmax(CN_map), np.nanmean(CN_map),  CN_map[CN_map != 100].mean()))
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)


# Inflow: Antecedent moisture
#==================#
#== SOLVE INFLOW ==# 
#==================#

if flag_plotinflow == 1:
    # Check: how to read precip correctly into numpy array (flip upside down if dims are time,lat,lon)
    precip_test = np.array(precip.mean('time'))     # the same as .values 
    precip_test = np.flipud(precip_test)            # flipud
    # --- Test --- 
    cax = plt.imshow(precip_test)                   # m/day
    plt.title('P mean m/day as array (to check spatial pattern)')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

# -------------------------------------------------------
# 1) Initialise runoff and antecedent moisture 3D arrays
# -------------------------------------------------------
Q_map = np.zeros((len(precip['time']), len(precip['lat']), len(precip['lon'])))     # numpy array 3D, runoff per pixel
AM_map = np.zeros((len(precip['time']), len(precip['lat']), len(precip['lon'])))    # numpy array 3D, antecedent moisture condition per pixel

# Check they are only zeros and there are no NaN 
check_zeroarray(Q_map, 'runoff')
check_zeroarray(AM_map, 'runoff')


# -------------------------------------------------------
# 2) Set the first 5 days of antecedent soil moisture to the initial AMC value
# -------------------------------------------------------
AM_map[0:amc_days] = amc_initialvalue

# mean AMC 0.0112 --> in Matlab this seems to be 0.0145, masking out the lake and the non-basin, here in Python seems to be 0.0167 \\ CHECK THIS \\

# Check initialisation is as desired, fxn: array to check, name of variable, desired value
check_valuesarray(AM_map[0:amc_days], 'antecedent moisture', amc_initialvalue)

if flag_plotinflow == 1:
    # --- Test ---                                                  
    cax = plt.imshow(AM_map[0])                                     
    plt.title('Antecedent moisture day 1')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)


# -------------------------------------------------------
# 3) Initialize empty data arrays of time-dependent curve numbers (modified by AMC), 
# soil water retention, and runoff per pixel
# -------------------------------------------------------
CN_AMC_da = np.zeros((len(precip['time']), len(precip['lat']), len(precip['lon']))) # numpy array 3D,curve number modified by AMC
S_da = np.zeros((len(precip['time']), len(precip['lat']), len(precip['lon'])))      # numpy array 3D, max soil moisture retention


# /////////// I FLIP HERE //////////// 

print('...calculating antecedent moisture')
# -------------------------------------------------------
# 4) Determine antecedent soil moisure for whole time period
# soil moisture (AMC) = the sum of the Ppn of the previous 5 days (P_5day)
# -------------------------------------------------------
for t in range(amc_days-1, len(DATEs)):                          # from 5 (6th item) to 8035-1 (last item)
    if t > amc_days-1:                                           # e.g. t = 5 (6th day)
        cumprecip = precip.isel(time=slice(t-amc_days,t))        # cumprecip, shape 5,130,130: precipitation in previous 5 days - e.g. from 0 to 4 (5 days)
        cummoisture = cumprecip.sum('time').values[::-1,:]        # cummoisture, shape 130,130: sum of the precip in the previous 5 days flipud on a 2d array flips along row axis, understand how to specify re 3d array                                                      # note. if you use precip_basin, becomes 0 outside of basin                                                                
        AM_map[t] = AM_map[t] + cummoisture                      # map of antecedent moisture, 3D [m.w.e.]

if flag_plotinflow == 1:        
    # --- Test --- 
    cax = plt.imshow(np.mean(AM_map, axis=0))                       # mean of this mean map = 0.0167 [m.w.e.] = 16.7 mm, including whole domain and lake -> check over only basin  
    plt.title('Antecedent moisture mean m.w.e.')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cax = plt.imshow(cummoisture)                                   # mean value of last cummoisture = 0.0163 m 
    plt.title('cumulative moisture day {} m.w.e'.format(t+1))
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # calc mean AM in basin and look at basin mask on CN 
    # --- Test --- 
    AM_basin = np.mean(AM_map,axis=0) * np.flipud(rmask_basinlake.values) 
    plt.imshow(AM_basin)
    plt.title('mean AM basinlake {:.6f}'.format(np.nanmean(AM_basin)))
    plt.colorbar(orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    AM_basin = np.mean(AM_map,axis=0) * np.flipud(rmask_basin.values) 
    plt.imshow(AM_basin)
    plt.title('mean AM basin {:.6f}'.format(np.nanmean(AM_basin)))  # 0.016285 m mean AM in basin
    plt.colorbar(orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cmap = plt.get_cmap('viridis', 10)   # discrete colors
    cax = plt.imshow(CN_map * np.flipud(rmask_basinlake.values), cmap = cmap)  
    plt.colorbar(cax, orientation='vertical')
    plt.title('CN map masked basinlake') 
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cmap = plt.get_cmap('viridis', 10)   # discrete colors
    cax = plt.imshow(CN_map * np.flipud(rmask_basin.values), cmap = cmap) 
    plt.colorbar(cax, orientation='vertical')
    plt.title('CN map maskedbasin ')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)  # --> outline of lake doesn't correspond, careful because this will affect inflow on the lake... 
                # I should change the values of these pixels to something other than 100, like take their nearest neighbor value... 
                # also check about the smaller lakes being treated as 100% runoff 

# OPTIONAL :             
# Expand the lake buffer a bit, get index numbers of non-nan values, compare with CN map. If any of the CN in those indices
# are equal to 100 set them to nearest non-100 value (i.e. turn them into land pixels, update CN_map)

# dilation: https://stackoverflow.com/questions/29027249/create-buffer-zone-within-a-numpy-array
# replace value by closest value below a threshold: https://gis.stackexchange.com/questions/207700/nearest-numpy-array-element-whose-value-is-less-than-the-current-element 

# Clear memory 
del cumprecip, cummoisture
gc.collect()

#%% Apply AM on CN
print('...applying antecedent moisture on curve numbers')

# -------------------------------------------------------------------------------------------------------------------------
# 5) Apply antecedent moisture condition on curve number : dry day (AMI), normal day (AMII), wet day (AMIII)
            # (Descheemaeker et al 2008) for AMI, AMII, AMIII conditions
            # (Ponce and Hawkins, 1996, Eq. 15-16) for CNI and CNIII formulas 
            # units AM: m 
# -------------------------------------------------------------------------------------------------------------------------
CN_AMC_da = np.where(AM_map < 0.0125, (CN_map / (2.281 - 0.01281 * CN_map))  , CN_map) # condition (dry day), where true (CN-I), where false (CN-II)
CN_AMC_da = np.where(AM_map >  0.0275,  (CN_map / (0.427 + 0.00573 * CN_map)) , CN_AMC_da) # condition (wet day), where true (CN-III), where false (whatever it was before, CN-I or CN-II)

# Each grid cell has 3 possible CN values depending on AM conditions

# -------------------------------------------------------------------------------------------------------------------------
# 6) Calculate S (max soil water retention), note where water: CN = 100 and S = 0 --> max CN, less soil retention, it all becomes runoff 
# (so I need to be careful to mask out the lake, as we already count the P_lake). 
# Where S=0, Q=P -->  So in fact P_lake is not strictly necessary, the CN method could calculate it by itself
            # (USDA, 2004 Eq 10-13)
            # units S: mm
# -------------------------------------------------------------------------------------------------------------------------
S_da = (25400 / CN_AMC_da) - 254 
             

if flag_plotinflow == 1:
    # --- Test --- 
    plt.set_cmap('viridis')
    cax = plt.imshow(CN_map) 
    plt.title('original CN map')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    
    # --- Test --- 
    cax = plt.imshow(np.nanmean(CN_AMC_da, axis=0))
    plt.title('time-dependent CN map mean')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)

if flag_plotCNdiff == 1:    
    # --- Test --- 
    CN_AMC_meandiff = np.nanmean(CN_AMC_da, axis=0) - CN_map
    cax = plt.imshow( CN_AMC_meandiff , vmin=-18, vmax=18)
    plt.title('Difference CN_AMC - original CN') # --> understand better the CNs and how AMC changes them and explain well in thesis 
    plt.set_cmap('bwr')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)           # --> CNs decrease over almost all of basin (i.e. S increases), they increase where there is most Ppn (S decreases there, makes sense I guess, the soil is saturated more often?)
    plt.set_cmap('viridis')
    
if flag_plotinflow == 1:    
    # --- Test --- 
    t_sel = 1
    cax = plt.imshow(CN_AMC_da[t_sel]) 
    plt.title('time-dependent CN map, day {}'.format(t_sel))
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    
    # --- Test --- 
    t_sel = 1
    cax = plt.imshow(S_da[t_sel]) 
    plt.title('max soil water retention S, day {} [mm]'.format(t_sel)) # Values of S look ok 
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    
    
    # --- Test --- 
    cax = plt.imshow(np.nanmean(S_da, axis=0)) 
    plt.title('max soil water retention S, mean [mm]')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)


#%% Calc Q per grid cell: New version using whole precip and not precipbasinlake and masking afterwards : Version 2 
# THIS IS VERY SLOW, IMPROVE

print('...calculating runoff Q_map')
# -------------------------------------------------------
# 7) Calculate outflow Q per grid cell (source: USDA 2004, Eq. 10-11) --> check if all is working, order of magnitude of Qin now is good but seasonality is not convincing 
    # condition: where P (mm) > 0.2*S
    # where true : calc Q
    # where false : Q = 0
# -------------------------------------------------------

# could use precip_raw instead but need to remove neg vals
precip_arr_mm = precip.values[:,::-1,:] * 1e3   

# Clear memory
del precip
gc.collect()

if flag_plotflip == 1:
    # --- Test --- 
    cax = plt.imshow(np.nanmean(precip_arr_mm, axis=0)) 
    plt.title('P mean mm/day as array flipped [mm]')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    
    # --- Test --- 
    cax = plt.imshow(np.nanmean(S_da, axis=0)) 
    plt.title('S mean mm/day as array [mm]')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)


# Calc Q part 2
Q_map = np.where( ( (precip_arr_mm) > (S_da * 0.2) ),  #  CHECK THIS CONDITION AND FORMULA MAKE IT INTO A FOR LOOP MAYBE... !!! removed : (S_da != 0) as this would get rid of rivers also (if they show in CN map, so I need to mask out lake surface), removed & (precip_arr_mm_test != 0) as this doesn't help re the divide by zero
                 ( (precip_arr_mm - 0.2*S_da)**2 / (precip_arr_mm + 0.8*S_da)  ), # I get warning: invalid value encountered in true_divide, check this, a divide by 0 error (if P and S are both = 0) but should give a nan and there are no nan, also this would lead to condition=False
                 0)                             # Q in mm l.w.e./day = kg m-2 day-1  

# \\\ CHECK FOR NaNs HERE \\\
    
# Clear memory
del precip_arr_mm
gc.collect()
    

print('...converting units of runoff Q_map')
# -------------------------------------------------------
# 8) Convert Q in mm/day per grid cell in basin to lake level equivalent m/day
        # kg m-2 day-1 (Q_map) * m2 (A_cell) = kg day-1
        # kg day-1 * 0.001 m3 kg-1 water = m3 / day (Q_map_m3) per grid cell
        # Masked over basin only, removing lake and outside of basin
        # Sum over all basin cells (nansum) = total m3 entering lake/day
        # Divide by lake area m2 = m lake level eq.
        
                # note.check alternate possibilities for A_lake
                # note. made this more compact for code running better
# -------------------------------------------------------
Q_map_masked_m3 = Q_map * inigeom.A_cell \
    * 1e-3 \
    * np.flipud(rmask_basin.values)                      
Qin_m = np.nansum(Q_map_masked_m3, axis=(1,2)) / inigeom.A_lake                             


# Plots inflow 
if flag_plotinflow == 1:
    # --- Test --- 
    cax = plt.imshow(np.flipud(rmask_basin.values)) 
    plt.title('check mask')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cax = plt.imshow(Q_map[10]) 
    plt.title('Qmap mm/day test')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cax = plt.imshow(Q_map_masked_m3[10]) 
    plt.title('Qmap m3/day test')
    plt.colorbar(cax, orientation='vertical')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cax = plt.imshow(np.nanmean(Q_map_masked_m3, axis=0)) # --> way more runoff in the water bodies (CN=100 all P becomes runoff)
    plt.title('Qmap m3/day mean (different scales)')
    plt.colorbar(cax, orientation='vertical')
    #c = plt.colorbar()
    plt.clim(0, 1e5) 
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
    # --- Test --- 
    cax = plt.imshow(np.nanmean(Q_map_masked_m3, axis=0)) 
    plt.title('Qmap m3/day mean (different scales)')
    plt.colorbar(cax, orientation='vertical')
    #c = plt.colorbar()
    plt.clim(0, 1e4) 
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)


#===============#
#== 1D: Q_in ===# 
#===============#

Qin_wb = Qin_m 
# 1D array

# --- Test --- 
plt.plot(Qin_wb)
plt.title('Qin_wb')
fig_showsave(plt, flag_savefig, fig_path, fig_n)
fig_n = update_fign(flag_savefig, fig_n)

# Clear memory
del Q_map_masked_m3 , Q_map, S_da, CN_AMC_da, AM_map
gc.collect()


#%% Calc Qout
print('...manipulating outflow')
#=====================#
#== 4) CALC OUTFLOW ==# 
#=====================#

if run_type == 0 and flag_obs_outflow == 0: # use multi-source outflow (not agreed curve) 
    # Cut to correct dates 
    outflow = outflow.loc[startDATE:endDATE] 
    
    # Check dates and Nan/negative values 
    if list(set(np.array(DATEs)) - set(np.array(outflow.index))): 
        print("model dates and outflow dates do not match")
    if np.any(outflow < 0) or np.any(outflow != outflow): 
        print("there are some negative or NaN values in your outflow \n ... please check them")
    
    
    # /////////// CAN CHOOSE LAKE AREA HERE (see inigeom) //////////// 
    
    #==============#
    #== 1D: Qout ==#   # dataframe 
    #==============#
    
    Qout_wb = outflow['outflow'] / inigeom.A_lake #try also A_lake_old // with A_lake mean = 0.001323 (higher), with A_lake_old mean = 0.001295 (lower), with A_shp mean = 0.001316 (somewhere inbetween) // difference about 1 mm I think? not a big deal
    
    # --- Test --- 
    plt.plot(Qout_wb)
    plt.title('Qout_wb')
    fig_showsave(plt, flag_savefig, fig_path, fig_n)
    fig_n = update_fign(flag_savefig, fig_n)
else:
    print('setting Qout to zero to calculate with agreed curve')
    Qout_wb = np.zeros(len(DATEs))

#==========================================#
#== EXTRA: COMPARE LAKE AREA ESTIMATES ====# 
#==========================================#

# Compare different A_lake estimates, this will affect Qin and Qout and differences in estimates suggest differences in how masking has been done
# See which one works best... Feed it two masked files so it can count the number of unmasked pixels and multiply by pixel area 

# compare_lake_areas(evap_lake[0], precip_lake[0]) # see function file to get information on different estimates and sources


#%% Calc LL
print('...calculating lake levels')
#=========================#
#== 5) CALC Lake Levels ==# 
#=========================#


# initialise empty lake level array 
L_wb = np.zeros(len(DATEs))

# Initialise first level in masl 
if run_type == 0: # OBS PPN, where I have lake levels
    L0 = lakelevels.loc[startDATE] 
    L_wb[0] = L0['water_level']
else: 
    lakelevels_slice = lakelevels.loc[lakelevels.index[0]:'1960-12-31'] 
    L0 = lakelevels_slice.mean()
    L_wb[0] = L0['water_level']

# Calculate lake levels time-stepping
for t in range(1,len(DATEs)):

    if run_type == 0 and flag_obs_outflow == 0: 
        # observational run, observational outflow  
        print('observed outflow being used to calculate lake levels')
        L_wb[t] = L_wb[t-1] - Qout_wb.values[t] + Qin_wb[t] - E_wb[t] + P_wb.values[t] # --> make these dtypes consistent 
    
    # For the test with the old NaN data
    # L_wb[t] = L_wb[t-1] - Qout_wb.values[t] + Qin_wb[t] - E_wb[t] + np.nan_to_num(P_wb.values)[t]

    else: 
        # Apply agreed curve formula each day to calculate Qout
        print('calculating agreed curve')
        L_insitu = L_wb[t-1] - 1123.319
        Qout_ac_cumecs = 66.3*((L_insitu - 7.96)**(2.01)) # formula from Vanderkelen et al 2018, from Sene 2000
        Qout_ac_cmd = Qout_ac_cumecs * inicon.sec_per_day  
        Qout_wb_ac = Qout_ac_cmd / inigeom.A_lake
        Qout_wb[t] = Qout_wb_ac
        if t == 1: 
            # fill first day
            Qout_wb[t-1] = Qout_wb[t]
        
        L_wb[t] = L_wb[t-1] - Qout_wb[t] + Qin_wb[t] - E_wb[t] + P_wb.values[t] # --> make these dtypes consistent 
            
    
print('...done !')

#%% Plot the levels 

plt.plot(L_wb)
plt.title('{} \n v0{}r0{}, DeltaL = {:.03f} m ({:.3f} m/yr)'.format(run_name, ver_n, run_n, (L_wb[len(DATEs)-1] - L_wb[0]), (L_wb[len(DATEs)-1] - L_wb[0]) /(endYEAR - startYEAR +1) ))
fig_showsave(plt, flag_savefig, fig_path, fig_n)
fig_n = update_fign(flag_savefig, fig_n)

#%% Save the output

data = {'date':  DATEs,
        'L_wb': L_wb,
        'P_lake': P_wb.values,
        'E_lake': E_wb,# ['E_lake']
        'Q_in': Qin_wb,
        'Q_out': Qout_wb,
        }

df = pd.DataFrame(data)
df = df.set_index('date')
df.to_csv(os.path.join(out_path, 'WBM_run_{}_v0{}r0{}_{}_{}_m.csv'.format(run_name, ver_n, run_n, startYEAR, endYEAR)))

data = { 'date':  DATEs, # maybe delete this extra column?? or keep as a check?
        'L_wb': L_wb,
        'P_lake': P_wb.values * 1e3,
        'E_lake': E_wb  * 1e3, #['E_lake']
        'Q_in': Qin_wb  * 1e3,
        'Q_out': Qout_wb  * 1e3,
        }

df_mm = pd.DataFrame(data)
df_mm = df_mm.set_index('date') # why was this not necessary before the AC ouflow??
df_mm.to_csv(os.path.join(out_path, 'WBM_run_{}_v0{}r0{}_{}_{}_mm.csv'.format(run_name, ver_n, run_n, startYEAR, endYEAR)), index=True)

# v01
# ----------
# run1 --> Qout was overwritten by Qin, looked quite good, rising lake levels but not super dramatic
# run2 --> Qin is much too high, something is wrong in the inflow calculations, fix this (it was a unit problem, fixed)
# run3 --> falling lake levels, by 6 meters, Qin is much too low, check this, get it approx to the magnitude of Qout
# run4 --> falling lake levels, by 8 meters, Qin is much too low, and I lowered P_lake by changing how I took the daily mean (I think P_lake is good now, we just need to get Qin up)
# run5 --> much better!! I rewrote the Qin calculation making it more explicit, but still dry
# run6 --> the same as run5 but initialised levels correctly, -75 mm/r
# run7 --> test 1983-2020, -111 mm/yr
# run6.5 --> the same result as run6, I changed the way of masking but no numerical effect apparent. 

# v02 (regionmasking and NaN data) :
# ----------------------
# run 1 (with old NaN data reordered to be read in correctly as numpy array as well) --> weirdly the same results as with previous masking..
# different NaN runs, coherence with matlab in P_lake and E_lake (when mean taken in same way)

# v03 observational: 
# ----------------------
        # Plake and Elake spatial mean is now taken over all cells, not first x then y
        # the lake and basin masks are made with regionmask in main script, always check in plots that they work!
        # Elake climatology is now correct 
# run 1 --> good, smaller bias (4.5 mm/month = 48 mm/year). P_lake is around 126 mm (ok), Evap has fallen a bit to 122 mm because of different way of taking spatial mean. Qin needs to be fixed. Qout ok 
# run 2 1983-2020 --> good ! 

# v04 observational: 
# ----------------------
# new outflow, AC earlier (from 2006)

# v03 historical:
# run 1 ISIMIP3b_hist_GFDL-ESM4 --> hadn't converted mm/s to mm/day, now done! it dries
# v02 hist-nat:
# it wets ! 

# version numbering starts again in version made for HPC, from v01 and each run is a different input data
# ----------------------


#%% Plot LL with dates 

df['L_wb'].plot()
plt.title('{} \n v0{}r0{}, DeltaL = {:.03f} m ({:.3f} m/yr)'.format(run_name, ver_n, run_n, (L_wb[len(DATEs)-1] - L_wb[0]), (L_wb[len(DATEs)-1] - L_wb[0]) /(endYEAR - startYEAR +1) ))
if flag_savefig == 1:
    plt.savefig(os.path.join(fig_path,'LL_{}_{}_{}_v0{}r0{}_dates.png'.format(run_name, startYEAR, endYEAR, ver_n, run_n)),dpi=300)


#%% Plot WB terms monthly acc

df_mm_resample = df_mm.resample('M', label='right').sum()
df_mm_resample.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"]) #
if flag_savefig == 1:
    plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_v0{}r0{}_M.png'.format(startYEAR, endYEAR, ver_n, run_n)),dpi=300)

#%% Plot WB terms yearly acc

fig, ax = plt.subplots(figsize=(7.5,4)) 

df_mm_resample = df_mm.iloc[:,1:5].resample('Y', label='right').sum()
df_mm_resample.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"], ax=ax)
df_residual = df_mm_resample["P_lake"] - df_mm_resample["E_lake"] + df_mm_resample["Q_in"] - df_mm_resample["Q_out"]
df_residual.plot(color='lightgray', ax=ax, label = 'residual' )
plt.axhline(y=0, color='gray', ls='--')
# plot formatting
plt.legend(loc='lower left')
ax.grid(True)
ax.set_ylabel('mm/year')
# add text
pos_x = 1.02
pos_y = 0.9
ax.text(pos_x, pos_y, 'yearly mean \n(mm/year):', transform=ax.transAxes)
summarystats = pd.DataFrame(round(df_mm_resample.mean(),1))
summarystats.columns = ['']
pos_y = 0.88
ax.text(pos_x, pos_y, str(summarystats),
          horizontalalignment='left',
          verticalalignment='top', transform=ax.transAxes)

mean_bias = df_mm_resample.mean()[0] - df_mm_resample.mean()[1] + df_mm_resample.mean()[2] - df_mm_resample.mean()[3]
pos_y = 0.02
ax.text(pos_x,pos_y,'{} \n({}-{}) v0{}r0{} \n \nmean residual {:.2f} mm/yr'.
        format(run_name, startYEAR, endYEAR,ver_n, run_n, mean_bias), transform=ax.transAxes)

fig.tight_layout()
if flag_savefig == 1:
    plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_{}_v0{}r0{}_Y_res.png'.
                             format(run_name, startYEAR, endYEAR, ver_n, run_n)),dpi=300)


#%% WBM terms climatology

fig, ax = plt.subplots(figsize=(7.5,4)) 

df_mm_resample = df_mm.iloc[:,1:5].resample('M', label='right').sum()
df_mm_climatology = df_mm_resample.groupby(df_mm_resample.index.month).mean() # climatology
df_mm_climatology.plot(y=["P_lake", "E_lake", "Q_in", "Q_out"], ax=ax)
ax.grid(True)
# add summary stats
pos_x = 1.02
pos_y = 0.9
ax.text(pos_x, pos_y, 'monthly mean climatology \n(mm/month):', transform=ax.transAxes)
summarystats = pd.DataFrame(round(df_mm_climatology.mean(),1))
summarystats.columns = ['']
pos_y = 0.88
ax.text(pos_x, pos_y, str(summarystats),
         horizontalalignment='left',
         verticalalignment='top',
         transform=ax.transAxes)
mean_bias = df_mm_climatology.mean()[0] - df_mm_climatology.mean()[1] + df_mm_climatology.mean()[2] - df_mm_climatology.mean()[3]
pos_y = 0.02
ax.text(pos_x, pos_y,'{} \n({}-{}) WBM v0{}r0{} \n \nmean residual {:.2f} mm/mo'.
         format(run_name, startYEAR, endYEAR,ver_n, run_n, mean_bias),
         transform=ax.transAxes)
fig.tight_layout()

if flag_savefig == 1:
    plt.savefig(os.path.join(fig_path,'WBM_terms_{}_{}_{}_v0{}r0{}_climatology.png'.format(run_name,startYEAR, endYEAR, ver_n, run_n)),dpi=200, bbox_inches = "tight")
