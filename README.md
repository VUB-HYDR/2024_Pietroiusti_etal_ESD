# Pietroiusti et al. 2024 (ESD) Possible role of anthropogenic climate change in the record-breaking 2020 Lake Victoria levels and floods

This repository contains the water balance model and analysis scripts used in Pietroiusti et al. 2024 "Possible role of anthropogenic climate change in the record-breaking 2020 Lake Victoria levels and floods" (accepted). 

The water balance model `lakevic-eea-wbm` simulates lake levels for Lake Victoria based on the following inputs: precipitation, evaporation, outflow and information on soil types and land use in the basin. The analysis scripts in `lakevic-eea-analysis` reproduce all results in the paper including analysing the main drivers of the high 2020 lake levels and applying a probabilistic extreme event attribution methodology to estimate the role of anthropogenic climate change in the 2020 high lake levels. 

<img src=/lakevic-eea-wbm/input_data/shapefiles/fig01.png alt="drawing" width="400" ALIGN=”center” />

Figure: Vanderkelen et al., 2018, HESS

## Water balance model

`lakevic-eea-wbm` contains the water balance model used in the study, it simulates daily lake levels for Lake Victoria.

The model can be run from the main `WBM_model.py` script. `WBM_settings.py` acts as a namelist and should be modified to run the model under different scenarios and forcings and to specify input data paths, `WBM_inicon.py` initiates constants in the model, `WBM_inigeom.py` imports shapefiles and geoinformation about Lake Victoria and the basin.

All input data necessary to run the model is in `lakevic-eea-wbm/input_data` and `lakevic-eea-wbm/modified_data`, except for precipitation (see 'Data availability'). 

The model takes as inputs:
1. Precipitation data, which should be remapped using the grid specifications in `mygrid.txt`. 
    - observational run: PERSIANN-CDR 1983-2020
    - attribution runs: CMIP6 GCM output historical and hist-nat bias adjusted with ISIMI3BASD from ISIMIP3b
2. Evaporation data `input_data/evap`
3. Information on soil classes, land use, basin and lake shapefiles `input_data/CN` and `input_data/shapefiles`
4. Outflow
    - observational run: semi-observational outflow `input_data/outflow` 
    - attribution runs: calculated by the model using the Agreed Curve.
  
The water balance model is based on the MATLAB version developed in [Vanderkelen et al. 2018, HESS, "Modelling the water balance of Lake Victoria"](https://hess.copernicus.org/articles/22/5509/2018/).

## Analysis of 2020 flood event

`lakevic-eea-scripts` contains the analysis scripts used in the attribution study. 

Scripts are divided into sections which correspond to the different sections of the analysis in the paper: 0 preprocessing, 1 analysis of observational precipitation data, 2 lake levels analysis and interpolation, 3 event definition, 4 evaluation of the WBM, 5 analysis of WB drivers, 6 GCM evaluation, 7 statistical fits for attribution and synthesis. 

## Data availability 

### Water balance model 

All data necessary to run the WBM, except for precipitation data, is provided in this GitHub repository (`lakevic-eea-wbm/input-data/`). 

This includes: 
1. Evaporation data `input_data/evap` and information on soil classes, land use, basin and lake shapefiles `input_data/CN` and `input_data/shapefiles` from [Vanderkelen et al. 2018, HESS](https://hess.copernicus.org/articles/22/5509/2018/).
2. Outflow `input_data/outflow`, which is a combination of the version in [Vanderkelen et al. 2018, HESS](https://hess.copernicus.org/articles/22/5509/2018/) and new data. 

Precipitation NetCDF files used in the paper are available publicly:
1. PERSIANN-CDR: https://www.ncei.noaa.gov/products/climate-data-records/precipitation-persiann 
2. ISIMIP3b CMIP6 GCMs https://data.isimip.org/
3. 
The location of precipitation input files should be specified in `WBM_settings.py`. 

### Extreme event analysis 

All data used in the analysis that can be publicly shared is available in this repository (input data in `lakevic-eea-scripts/data/input-data/` and outputs of the WBM runs that are behind the paper in `lakevic-eea-wbm/output/2022`).

`lakevic-eea-scripts/data/input-data/` contains the following: 
1.	IOD timeseries from NOAA, 
2.	Annual GMST time series from ISIMIP3b bias-adjusted CMIP6 data.
3.	GMST time series from GISTEMP, 
4.	Lake level timeseries based on a combination of DAHITI data and in situ measurements from *Vanderkelen et al. 2018, HESS*

All data included in this repository is also archived on  [Zenodo](https://zenodo.org/record/8233523).

## Version
Version November 2023

## Author
Rosa Pietroiusti

## References 

If you use this 
