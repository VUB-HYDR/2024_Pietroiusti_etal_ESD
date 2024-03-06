# Pietroiusti et al. 2024 ESD

This repository contains the water balance model and analysis scripts used in Pietroiusti et al. 2024 "Possible role of anthropogenic climate change in the record-breaking 2020 Lake Victoria levels and floods" (accepted). 

The water balance model simulates lake levels for Lake Victoria based on the following inputs: precipitation, evaporation, outflow and information on soil types and land use in the basin. The analysis scripts reproduce all results in the paper including analysing the main drivers of the high 2020 lake levels and applying a probabilistic extreme event attribution methodology to estimate the role of anthropogenic climate change in the 2020 high lake levels. 

<img src=/lakevic-eea-wbm/input_data/shapefiles/fig01.png alt="drawing" width="400" ALIGN=”left” />

Figure: Vanderkelen et al., 2018, HESS

## Water balance model

`lakevic-eea-wbm` contains the water balance model used in the study, and can be run from the main `WBM_model.py` script. 

This model takes as input:
1. Precipitation data, remapped using the provided grid specifications (`mygrid.txt`). 
    - observational run: PERSIANN-CDR 1983-2020
    - attribution runs: CMIP6 GCM output historical and hist-nat bias adjusted with ISIMI3BASD from ISIMIP3b
2. Evaporation data (provided) 
3. Information on soil classes, land use, basin and lake shapefiles (provided)
4. Outflow
    - observational run: semi-observational outflow (provided)
    - attribution runs: calculated using the Agreed Curve. 

`WBM_settings.py` acts as a namelist and should be modified to run the model under different scenarios and forcings and to specify input data paths, `WBM_inicon.py` initiates constants in the model, `WBM_inigeom.py` imports shapefiles and geoinformation about Lake Victoria and the basin.

The water balance model is based on the MATLAB version developed in Vanderkelen et al. 2018, HESS, "Modelling the water balance of Lake Victoria".

All input data necessary to run the model is in the `lakevic-eea-wbm/input_data` and `lakevic-eea-wbm/modified_data` directories (see 'Data Availability' below for details on where this data comes from), except for the precipitation data (see 'Data availability' below for details on where to get this data). 

## Analysis of 2020 flood event

`lakevic-eea-scripts` contains the analysis scripts used inthe attribution study. 

Scripts are divided into sections which correspond to the different sections of the analysis in the paper: 0 preprocessing, 1 analysis of observational precipitation data, 2 lake levels analysis and interpolation, 3 event definition, 4 evaluation of the WBM, 5 analysis of WB drivers, 6 GCM evaluation, 7 statistical fits for attribution and synthesis. 

## Data availability 

### Water balance model 

All data necessary to run the WBM is provided in this GitHub repository (`lakevic-eea-wbm/input-data/`), except for precipitation NetCDF files, which are available publicly, here:
1. PERSIANN-CDR: https://www.ncei.noaa.gov/products/climate-data-records/precipitation-persiann 
2. ISIMIP3b CMIP6 GCMs https://data.isimip.org/

The location of precipitation input files should be specified in `WBM_settings.py`. 

### Extreme event analysis 

All analysis data that can be publicly shared is available in this repository (input data in `lakevic-eea-scripts/data/input-data/` and outputs of the WBM runs that are behind the paper in `lakevic-eea-wbm/output/2022`).

All data included in this repository is also archived on  [Zenodo](https://zenodo.org/record/8233523).

## Version
Version November 2023

## Author
Rosa Pietroiusti
