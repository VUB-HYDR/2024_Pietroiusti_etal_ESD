#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:47:50 2022

@author: rosapietroiusti

WBM: initialise geometry
"""

# import modules
import sys, os , glob
import geopandas as gpd # check consistency of name import in different files 
import numpy as np
import regionmask
import xarray as xr
import matplotlib.pyplot as plt

# link to other scripts
WBM_path = os.path.dirname(__file__) #'/Users/rosapietroiusti/Documents/Thesis/WBM_pythonversion'
sys.path.append(WBM_path)
import WBM_inicon as inicon 
import WBM_settings as settings 


# initialise Lake Victoria boundaries - do I need this ? 
# lat_min_Vict = -3.1
# lat_max_Vict = 0.6
# lon_min_Vict = 31.4
# lon_max_Vict = 35
# bounds_Vict  = [lat_min_Vict, lat_max_Vict, lon_min_Vict, lon_max_Vict]


# define size of a grid point (grid resolution)
res_grid = 0.065 # pixel length [°]
res_m = res_grid*inicon.c_earth/360 # pixel length [m] - at the equator?
A_cell = res_m**2 # pixel area - note, using this is inconsistent with A_lake being calculated by projecting shapefile and then calculating area


# Calculate area of lake - projecting shapefiles (compare these values with literature)
lake_shp = gpd.read_file(settings.filepath_shp_lake, crs="epsg:4326") # check this is correct crs, check what file this is in
basin_shp = gpd.read_file(settings.filepath_shp_basin, crs="epsg:4326")
lake_proj = lake_shp.to_crs("EPSG:32736") # UTM zone 36S
basin_proj = basin_shp.to_crs("EPSG:32736")
A_lake = sum(lake_proj.area) # m^2
A_basin = sum(basin_proj.area) # m^2


# Alternative A_lake values 
A_lake_old = 68272645811.8022 # m^2 # (Vanderkelen 2018) was calculated as lake pixel number x A_cell # note pixel number of lake mask are slightly different in python version here cfr to origianl MATLAB version
A_lake_shp = lake_shp['Lake_area'].iloc[0] * 10**6 # from km2 to m2 # shapefile metadata


# Get grid coordinates from Ppn file --> might not need this if masks are made in main script
# var = 'precipitation'
# with xr.open_dataset(settings.filepath_precip, decode_coords="all") as ds:
#     lats_precip = np.array(ds.get('lat'))
#     lons_precip = np.array(ds.get('lon'))


# Make lake and basin masks  --> I moved this to main script, it was making problems here (small difs in lat/lon probably)
#rmask_lake = regionmask.mask_geopandas(lake_shp, lons_precip, lats_precip) + 1       # only lake set = 1
#rmask_basinlake = regionmask.mask_geopandas(basin_shp, lons_precip, lats_precip) + 1 # lake+ basin set = 1
#rmask_basin = rmask_basinlake * np.where(np.isnan(rmask_lake),1,np.nan) # only basin set = 1




