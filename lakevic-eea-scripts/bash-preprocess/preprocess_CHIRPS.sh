#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --ntasks=3
#SBATCH --mail-type=ALL

#==============================================================================
# PREPROCESS CHIRPS data 

# crops daily CHIRPS precipitation data to AGL, merges through time and reprojects to own grid.

# 1 day temporal resolution, 0.05 deg spatial resolution (approx 5.5km)
# CHIRPS version 2.0 currently on Hydra years 1981-2015

# Created: 26 April 2022. 
# By: Rosa Pietroiusti
#==============================================================================

# tic
START=$(date +%s.%N)

# load modules
module purge && module load CDO/1.9.10-intel-2019b


#==============================================================================
# settings
#==============================================================================

# run settings
flags_exper=1
              # 1: CHIRPS 1981-2015 
              # 2: fill when new data available


# set input directories and start and end year of simulations
if [ $flags_exper -eq 1 ]
then 
	inDIR=/data/brussel/vo/000/bvo00012/data/dataset/chirps/v2/0.05deg_lat-lon_1d/original
	startYEAR=1981
	endYEAR=2015
elif [ $flags_exper -eq 2 ]
then
	inDIR=
	startYEAR=
	endYEAR=
fi

echo **inDIR is $inDIR, period is $startYEAR to $endYEAR

# set output directory
outDIR=$VSC_DATA_VO_USER/Thesis/CHIRPS

# set working directory for intermediate files 
workDIR=$VSC_DATA_VO_USER/Thesis/CHIRPS/temp

if [ ! -e $workDIR ]
then
	mkdir $workDIR
fi
cd $workDIR

echo **workDIR is $workDIR

#==============================================================================
# initialisation
#==============================================================================

# set data name
dNAME=chirps
echo **preprocessing ${data} data

#==============================================================================
# grid for reprojection
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


# ==============================================================================
# processing : crop all files, merge all timesteps, reproject files
# ==============================================================================


# loop through files to crop files, merge files, reproject files

echo **cropping $dNAME files
for FILE in $inDIR/${dNAME}-*.????.*.nc; do
	FILENAME=$(basename -s .nc ${FILE})
	cdo -O sellonlatbox,20,40,-15,10 $FILE $workDIR/${FILENAME}_cropped.nc
done


echo **merging $dNAME files
cdo -O mergetime $workDIR/*_cropped.nc $workDIR/${dNAME}_v2.2_pr_AGL_${startYEAR}_${endYEAR}.nc 


echo **reprojecting $dNAME files
cdo -O -L -remapcon2,mygrid -sellonlatbox,27.5,38.5,-6,3 $workDIR/${dNAME}_v2.2_pr_AGL_${startYEAR}_${endYEAR}.nc $outDIR/${dNAME}_pr_owngrid_${startYEAR}_${endYEAR}.nc



#==============================================================================
# end
#==============================================================================

# toc
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo
