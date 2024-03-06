# Extreme event attribution of 2020 Lake Victoria floods

This repository contains the water balance model and analysis scripts used in [**Pietroiusti et al. 2024 "Possible role of anthropogenic climate change in the record-breaking 2020 Lake Victoria levels and floods"**](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-1827/) (accepted). 

The water balance model `lakevic-eea-wbm` simulates lake levels for Lake Victoria based on the following inputs: precipitation, evaporation, outflow and soil types and land cover in the basin. The analysis scripts in `lakevic-eea-analysis` reproduce all results in the paper including analysing the main drivers of the high 2020 lake levels and applying a probabilistic extreme event attribution methodology to estimate the role of anthropogenic climate change. 

<img src=/lakevic-eea-wbm/input_data/shapefiles/fig01.png alt="drawing" width="400" ALIGN=”center” />
Figure source: Vanderkelen et al. 2018

## Water balance model

`lakevic-eea-wbm` contains the water balance model used in the study, it simulates daily lake levels for Lake Victoria.

The model can be run from the main `WBM_model.py` script. `WBM_settings.py` acts as a namelist and should be modified to run the model under different scenarios and forcings and to specify input data paths, `WBM_inicon.py` initiates constants in the model, `WBM_inigeom.py` imports shapefiles and geoinformation about Lake Victoria and the basin.

All input data necessary to run the model is in `lakevic-eea-wbm/input_data` and `lakevic-eea-wbm/modified_data`, except for precipitation (see 'Data availability'). The model takes as inputs:
1. Precipitation data, which should be remapped using the grid specifications in `mygrid.txt`. 
    - observational run: PERSIANN-CDR 1983-2020
    - attribution runs: CMIP6 GCM output bias-adjusted within ISIMIP3b
2. Evaporation data `input_data/evap`
3. Information on soil classes, land use, basin and lake shapefiles `input_data/CN` and `input_data/shapefiles`
4. Outflow
    - observational run: semi-observational outflow `input_data/outflow` 
    - attribution runs: calculated by the model using the Agreed Curve.
  
The model is based on the MATLAB version from [Vanderkelen et al. 2018, HESS, "Modelling the water balance of Lake Victoria"](https://hess.copernicus.org/articles/22/5509/2018/).

## Extreme event analysis

`lakevic-eea-analysis` contains the analysis scripts used in the attribution study. 

Scripts are divided into sections which correspond to the different sections of the analysis in the paper: 

0. preprocessing
1. analysis of observational precipitation data
2. lake levels analysis and interpolation
3. event definition
4. evaluation of the WBM
5. analysis of WB drivers
6. GCM evaluation
7. statistical fits for attribution and synthesis

## Data availability & sources

### Data for water balance model 

All data necessary to run the WBM, except for precipitation data, are provided in this GitHub repository `lakevic-eea-wbm/input-data/`. This includes: 
1. Evaporation data `input_data/evap` from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/), originally from [Thiery et al. 2015](https://journals.ametsoc.org/view/journals/clim/28/10/jcli-d-14-00565.1.xml)
2. Information on soil classes and land cover `input_data/CN` from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/), originally from [JRC](https://publications.jrc.ec.europa.eu/repository/handle/JRC24914) and [Dewitte et al. 2013](https://www.sciencedirect.com/science/article/abs/pii/S0016706113002401?via%3Dihub) 
3. basin and lake shapefiles `input_data/shapefiles` from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/) and [Harvard dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/PWFW26)
4. Outflow `input_data/outflow` which is a combination of data from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/) and new in situ data. 

Precipitation NetCDF files used in the paper are available publicly:
1. PERSIANN-CDR from [NOAA](https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ncdc:C00854/html#)
2. ISIMIP3b bias-adjusted CMIP6 GCM output from [ISIMIP](https://doi.org/10.48364/ISIMIP.842396.1)

### Data for extreme event analysis 

All data used in the analysis that can be publicly shared is available in this repository. Input data in `lakevic-eea-analysis/data/input-data/` and outputs of the WBM runs that are behind the paper in `lakevic-eea-wbm/output/2022`. 

`lakevic-eea-analysis/data/input-data/` contains the following: 
1.	`IOD` timeseries from [NOAA](https://psl.noaa.gov/gcos_wgsp/Timeseries/DMI/)
2.	`gmst-models`: annual GMST time series from ISIMIP3b bias-adjusted CMIP6 data from [ISIMIP](https://doi.org/10.48364/ISIMIP.842396.1)
3.	`gmst-obs`: GMST time series from [GISTEMP](https://data.giss.nasa.gov/gistemp/)
4.	`lakelevels`: lake levels from [DAHITI](https://dahiti.dgfi.tum.de/en/products/water-level-altimetry/) and in situ measurements from [Vanderkelen et al. 2018, HESS](https://hess.copernicus.org/articles/22/5509/2018/)

### Archive

All data included in this repository is also archived on  [Zenodo](https://zenodo.org/record/8233523).

## Version
Version November 2023

## Author
Rosa Pietroiusti

## Citing  
If you use this code or data created by this study please cite [**Pietroiusti et al. 2024**](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-1827/). If you use data in this repository coming from other sources, please cite them, all DOI/links are present in this documentation. 
