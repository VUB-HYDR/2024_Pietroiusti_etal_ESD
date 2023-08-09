#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:41:48 2022

@author: rosapietroiusti

WBM: Initialise Physical constants
"""

# specific gas constant for water vapour (Van Lipzig, 2009 (W&K))
Rv = 461.5 # [J*kg^-1*K^-1]

# Specific gas constant for dry air (Van Lipzig, 2009 (W&K))
Rd = 287.058 # [J*kg^-1*K^-1]

#Water vapour pressure at freezing point of water (611Pa) (Van Lipzig, 2009 (W&K));
es0 = 611 # [Pa] !dit is op basis van de formule van de verzadigde dampspanning berekend!

# Temperature at freezing point of water (611Pa) (Van Lipzig, 2009)
T0 = 273.16 # [K]

# Standard reference pressure (Akkermans et al., 2014)
P0 = 100000 # [Pa]

# Latent heat of vaporization, assumed constant
Lvap = 2.50E6 # [J*kg^-1]

# Stephan-Boltzmann constant
sigma = 5.67E-8 # [W*m^-2*K^-4]

# Specific heat capacity of water at constant pressure
cp  = 4.1813E3 # [J kg^-1 K^-1]

# Specific heat capacity of air at constant pressure (Akkermans et al., 2014)
cp_air  = 1005 # [J kg^-1 K^-1]

# Latent heat of vaporization, assumed constant
rhow = 1E3 # [J*kg^-1]

# earth circumference (wikipedia)
c_earth = 40075017 # [m]

# earth gravitation (wikipedia)
g_earth = 9.81 # [m s^-1]

# seconds per day
sec_per_day = 60 * 60 * 24 # [s d^-1]
