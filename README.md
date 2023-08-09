# 2023_Pietroiusti_etal_ESD

Water balance model and analysis scripts used in Pietroiusti et al. 2023 Lake Victoria flood attribution study. Data used in the analysis is available on Zenodo: [INSERT LINK] or by request.

## Description 

### lakevic-eea-wbm

lakevic-eea-wbm contains the water balance model component of the analysis, and can be run from the main WBM_model.py script. 

WBM_settings.py acts as a namelist and should be modified to run the model under different scenarios and forcings, WBM_inicon.py initiates constants in the model, WBM_inigeom.py imports shapefiles and geoinformation about Lake Victoria and the basin.

The water balance model is based on the MATLAB version developed in Vanderkelen et al. 2018, HESS, "Modelling the water balance of Lake Victoria".

### lakevic-eea-scripts

lakevic-eea-scripts contains the analysis scripts. 

## Version
Version August 2023
