#!/bin/bash

#==============================================================================
#PERSIANN REPROJECT SCRIPT
#Adapted from Inne's script, reprojects AGL PERSIANN CDR files for all years to 
#the Lake Victoria WBM model grid
#Created: 25 November 2021
#By: Rosa Pietroiusti
#==============================================================================

# load modules
module purge && module load CDO/1.9.10-intel-2019b

#==============================================================================
# initialisation
#==============================================================================

# set working directory
cd $VSC_DATA/Thesis/PERSIANN/output

#==============================================================================
# manipulations
#==============================================================================

cat > mygrid << EOF
gridtype = lonlat
xsize    = 130
ysize    = 130
xfirst   = 28
xinc     = 0.065
yfirst   = -5.5
yinc     = 0.065
EOF

startYEAR=1983
endYEAR=2020
# for Inne this was 1983-2014

cdo -L -remapcon2,mygrid -sellonlatbox,27.5,38.5,-6,3 PERSIANN-CDR_v01r01_${startYEAR}_${endYEAR}_AGL.nc PERSIANN-CDR_v01r01_${startYEAR}_${endYEAR}_remapped_owngrid.nc

