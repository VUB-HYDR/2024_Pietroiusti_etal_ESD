#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --ntasks=3
#SBATCH --mail-type=ALL

#==============================================================================
# PERSIANN REORDER DIMENSIONS FOR WBM_Py to time, lat, lon (t,r,c)
# For easy reading into numpy arrays from xarrays 
# from: https://code.mpimet.mpg.de/boards/2/topics/11747 

# Created: 13 May 2022
#==============================================================================

# tic
START=$(date +%s.%N)

# load modules
module purge && module load NCO

#==============================================================================
# initialisation
#==============================================================================


# set work directory
workDIR=$VSC_DATA_VO_USER/Thesis/PERSIANN/modifyPERSIANN_old

infile=PERSIANN_old_1992_2014_owngrid_resample.nc #file from Inne, with NaN instead of zero
fileNAME=$(basename -s .nc ${infile}) 

varNAME=precipitation

#==============================================================================
# manipulations
#==============================================================================

# reorder the dimensions of the variables in the file 

ncpdq --rdr=time,lat,lon $workDIR/$infile $workDIR/${fileNAME}_temp.nc

# reorder the dimensions of the file to time,lat,lon and append variables 

ncks -A -v time $workDIR/${fileNAME}_temp.nc workDIR/${fileNAME}_reorder_timelatlon.nc
ncks -A -v lat $workDIR/${fileNAME}_temp.nc workDIR/${fileNAME}_reorder_timelatlon.nc
ncks -A -v lon $workDIR/${fileNAME}_temp.nc workDIR/${fileNAME}_reorder_timelatlon.nc
ncks -A -v $varNAME $workDIR/${fileNAME}_temp.nc workDIR/${fileNAME}_reorder_timelatlon.nc

# for some reason these last commands didn't work if i ran the script in HPC as "bash reorderdims....sh" so I did them by hand
# it gave a permission denied to create the outfile error, understand why 
# i think i had forgotten $ behind "workDIR"
