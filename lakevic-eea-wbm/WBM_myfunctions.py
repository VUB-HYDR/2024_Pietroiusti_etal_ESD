#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 16:11:31 2022

@author: rosapietroiusti
"""
import numpy as np

import WBM_inigeom as inigeom # why don't I need WBM_path here? 
import WBM_inicon as inicon
import os

#====================#
#== MY FUNCTIONS ====# 
#====================#

def compare_lake_areas(masked_ncfile1, masked_ncfile2):
    # Compare different A_lake estimates, this will affect P_lake and differences suggest differences in how masking has been done
    print('\n A_lake estimates [km2]')
    print('calculated here projecting shp', inigeom.A_lake / 10**6) # from projected shapefile projected and calculated with geopandas
    print('original area', inigeom.A_lake_old / 10**6) # from Inne's script (number of lake cells in MATLAB * area of cell)
    print('calculated here with original method', np.nansum(masked_ncfile1 == masked_ncfile1) * inigeom.A_cell / 10**6, '(file1)') # from number of masked grid cells (here) using rioxarray * area of cell from Inne's script
    print('from shapefile metadata', inigeom.A_lake_shp)

    print('\n number of lake pixels')
    print('here', np.nansum(masked_ncfile1 == masked_ncfile1), '(file2) or', np.nansum(masked_ncfile2 == masked_ncfile2), '(file2)'  )
    print('here if all_touched in rio.clip is set to True', 1549)
    print('WBM islake', 1304)
    print('WBM islake_intp', 1392)

    print('\n basin area here [km2] ', inigeom.A_basin / 10**6)
    return

def check_zeroarray(array, varname):
    # Check they are only zeros and there are no NaN 
    if np.any(array != 0): 
        print("there are some non-zero values in your empty initialised", varname, " \n... check them")
    if np.any(array != array): 
        print("there are some Nan values in your empty initialised", varname, " \n... check them")
    return 

def check_valuesarray(array, varname, value):
    # Check they are only desired value and there are no NaN in your array
    if np.any(array != value): 
        print("there are some non-desired values in your empty initialised", varname, " \n... check them")
    if np.any(array != array): 
        print("there are some Nan values in your empty initialised", varname, " \n... check them")
    return 

def check_dimorders(ds):
    # check your netcdf file (esp for precipitation) is in the correct dimension order to be correctly read as a numpy array
    if list(ds.dims) != ['time', 'lat', 'lon']:
        print("\n                  ***WARNING*** \nThe order of your dimensions in the file {} is {} not ['time', 'lat', 'lon'] \nif this is your PRECIPITATION file, this will create errors in the model \n".format( list(ds.keys()), list(ds.dims)   ))
        if list(ds.dims) == ['time', 'lon', 'lat']:
            flag_orderdims = 1 

def check_punits(da): 
    # Checks units of the precip file and if they are labelled as kg m-2 s-1 (mm/s) it modified them to kg m-2 d-1 (mm/d), otherwise it leaves the array as it is and tells you what units you have
    if da.units != 'kg m-2 d-1':
        print("... check your precip units")
        if da.units == 'kg m-2 s-1': 
            print('units are kg m-2 s-1, converting them to kg m-2 d-1')
            return da * inicon.sec_per_day
        else:
            print('units are {}'.format(da.units))
            return da 

def fig_showsave(plt, flag_savefig, fig_path, fig_n):
    if flag_savefig == 1:
        plt.savefig(os.path.join(fig_path, "fig_{:02d}".format(fig_n)))
        plt.clf()
    else: 
        plt.show() # check this is ok 

def update_fign(flag_savefig,fig_n):
    if flag_savefig == 1:
        return fig_n + 1
    else:
        return fig_n
        
