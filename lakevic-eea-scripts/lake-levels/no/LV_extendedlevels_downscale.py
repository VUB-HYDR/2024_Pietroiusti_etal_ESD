# -*- coding: utf-8 -*-
"""
Created on Thu May  4 13:06:18 2023

@author: rpietroi

Downscale the extended lakelevels to check if the positive trend was due to different resolution of timeseries 
"""



#=============#
#== IMPORTS ==#
#=============#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os , glob
import math
import pandas as pd
from matplotlib import transforms
from matplotlib.ticker import MaxNLocator
import scipy as scipy
from scipy import stats
from matplotlib.pyplot import cm

#===========#
#== DATA  ==#
#===========#

inDIR_Obs = r'C:\Users\rpietroi\OneDrive - Vrije Universiteit Brussel\repos_cloud\lakevic-eea\lakevic-eea-scripts\data\data-modified\lakelevels'

filepath = os.path.join(inDIR_Obs, 'lakelevel_ext_intr_HCDH22.csv')

df = pd.read_csv(filepath,
                     index_col = 'date',
                     parse_dates=True)

#%%


df_firstday_intr = df[df.index.day == 1].resample('D').asfreq().interpolate(method='linear').round(3) 

