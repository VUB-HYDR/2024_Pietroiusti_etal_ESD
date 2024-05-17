# Extreme event attribution of 2020 Lake Victoria floods

[![DOI](https://zenodo.org/badge/676547448.svg)](https://zenodo.org/doi/10.5281/zenodo.10794481)


This repository contains the water balance model and analysis scripts used in [**Pietroiusti et al. 2024 "Possible role of anthropogenic climate change in the record-breaking 2020 Lake Victoria levels and floods"**](https://esd.copernicus.org/articles/15/225/2024/) (Earth System Dynamics)

The water balance model `lakevic-eea-wbm` simulates lake levels for Lake Victoria based on the following inputs: precipitation, evaporation, outflow and soil types and land cover in the basin. 

The analysis scripts in `lakevic-eea-analysis` reproduce all results in the paper including analysing the main drivers of the high 2020 lake levels and applying a probabilistic extreme event attribution methodology to estimate the role of anthropogenic climate change in the event. 

<img src=/lakevic-eea-wbm/input_data/shapefiles/fig01.png alt="drawing" width="400" ALIGN=”center” />
Figure source: Vanderkelen et al. 2018

## Water balance model

`lakevic-eea-wbm` contains the water balance model used in the study, it simulates daily lake levels for Lake Victoria.

The model can be run from the main `WBM_model.py` script. `WBM_settings.py` acts as a namelist and should be modified to run the model under different scenarios and forcings and to specify input data paths, `WBM_inicon.py` initiates constants in the model, `WBM_inigeom.py` imports shapefiles and geoinformation about Lake Victoria and the basin.

The model takes as inputs:
1. Precipitation data, which should be remapped using the grid specifications in `mygrid.txt`. 
2. Evaporation data `input_data/evap`
3. Information on soil classes, land use, basin and lake shapefiles `input_data/CN` and `input_data/shapefiles`
4. Outflow
    - observational run: semi-observational outflow `input_data/outflow` 
    - attribution runs: calculated by the model using the Agreed Curve.

The model is based on the MATLAB version from [Vanderkelen et al. 2018, HESS, "Modelling the water balance of Lake Victoria"](https://hess.copernicus.org/articles/22/5509/2018/).


## Data for water balance model 

For ease of reproducibility, all data necessary to run the WBM, except for precipitation data, are provided in this GitHub repository.  

`lakevic-eea-wbm/input-data/` and `lakevic-eea-wbm/modified_data` contain:
1. Evaporation data `input_data/evap` from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/), originally from [Thiery et al. 2015](https://journals.ametsoc.org/view/journals/clim/28/10/jcli-d-14-00565.1.xml)
2. Information on soil classes and land cover `input_data/CN` from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/), originally from [JRC](https://publications.jrc.ec.europa.eu/repository/handle/JRC24914) and [Dewitte et al. 2013](https://www.sciencedirect.com/science/article/abs/pii/S0016706113002401?via%3Dihub) 
3. basin and lake shapefiles `input_data/shapefiles` from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/) and [Harvard dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/PWFW26)
4. Outflow `input_data/outflow` which is a combination of data from [Vanderkelen et al. 2018](https://hess.copernicus.org/articles/22/5509/2018/) and new in situ data.

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

## Data availability

All data necessary to reproduce analyses are archived on Zenodo ([here](https://zenodo.org/record/8233523), DOI: 10.5281/zenodo.8229505).

For ease of reproducibility, data necessary to run the water balance model is also provided in this repo. 

## Version
Version November 2023

## Author
Rosa Pietroiusti

## References  
Pietroiusti, R., Vanderkelen, I., Otto, F. E. L., Barnes, C., Temple, L., Akurut, M., Bally, P., van Lipzig, N. P. M., and Thiery, W.: Possible role of anthropogenic climate change in the record-breaking 2020 Lake Victoria levels and floods, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1827, 2023. 

