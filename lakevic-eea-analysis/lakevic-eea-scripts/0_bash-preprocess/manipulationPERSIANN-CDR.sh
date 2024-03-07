#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --ntasks=3
#SBATCH --mail-type=ALL

#==============================================================================
# PERSIANN MERGE AND CROP SCRIPT
# Adapted from Inne's script, merges daily PERSIANN-CDR data and crops over AGL domain
# Created: 20 November 2021 
#==============================================================================

# tic
START=$(date +%s.%N)

# load modules
module purge && module load CDO/1.9.10-intel-2019b

#==============================================================================
# initialisation
#==============================================================================


# set output directory
outDIR=/data/brussel/104/vsc10419/Thesis/PERSIANN/output

# set input directory TEST
#inDIR=/data/brussel/104/vsc10419/Thesis/PERSIANN/test

# PERSIANN data folder REAL
inDIR=/data/brussel/vo/000/bvo00012/data/dataset/PERSIANN-CDR/

# set working directory for intermediate files 
cd $VSC_SCRATCH/Thesis/PERSIANN

#==============================================================================
# manipulations
#==============================================================================

# merge daily files to create a yearly file for each year and copy it in my scratch folder
echo "**MERGING DAILY FILES**"

startYEAR=1983
endYEAR=2020
for ((YEAR=startYEAR;YEAR<=endYEAR;YEAR++));
	do
		echo "processing $YEAR"
		cdo mergetime $inDIR/${YEAR}/PERSIANN-CDR_v01r01*.nc PERSIANN-CDR_v01r01_${YEAR}.nc
	done


# cut each yearly file to AGL 
echo "**CROPPING YEARLY FILES TO AGL**"

for ((YEAR=startYEAR;YEAR<=endYEAR;YEAR++));
        do
		echo "processing $YEAR"
		cdo sellonlatbox,20,40,-15,10 PERSIANN-CDR_v01r01_${YEAR}.nc AGL_${YEAR}.nc
	done
 

# merge all yearly AGL files 
echo "**MERGING ALL YEARLY AGL FILES AND COPYING OUTPUT**"

cdo mergetime AGL_????.nc PERSIANN-CDR_v01r01_${startYEAR}_${endYEAR}_AGL.nc

cp PERSIANN-CDR_v01r01_${startYEAR}_${endYEAR}_AGL.nc $outDIR/PERSIANN-CDR_v01r01_${startYEAR}_${endYEAR}_AGL.nc

# toc
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo
