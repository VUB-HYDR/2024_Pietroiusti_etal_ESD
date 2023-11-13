#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --ntasks=3
#SBATCH --mail-type=ALL

#==============================================================================
# PREPROCESS ISIMIP3a data 
# crops daily ISIMIP3a precipitation data to AGL, merges through time (1850-2020) and reprojects to own grid
# Created: 6 April 2022. Modified: 11 May 2022. Adapted from ISIMIP3b script 20 June 2023. 
# By: Rosa Pietroiusti rosa.pietroiusti@vub.be
#==============================================================================

# tic
START=$(date +%s.%N)

# load modules
module purge && module load CDO/2.0.6-gompi-2022a #CDO/1.9.10-intel-2019b


#==============================================================================
# settings
#==============================================================================

# run settings
flags_exper=2
              # 1: obsclim
              # 2: counterclim
flags_data=1
	      # 1: 20CRv3-ERA5 
          # other data when available


# set input directories, start and end year of simulations and name to save
if [ $flags_data -eq 1 ]
then 
	startYEAR=1901
	endYEAR=2020
    data_name=('20CRv3-ERA5')
	if [ $flags_exper -eq 1 ]
	then
		inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3a/InputData/climate/atmosphere/obsclim/global/daily/historical/20CRv3-ERA5/ 
	elif [ $flags_exper -eq 2 ]
	then
		inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3a/InputData/climate/atmosphere/counterclim/global/daily/historical/20CRv3-ERA5/
	fi
fi

echo **inDIR is $inDIR, period is $startYEAR to $endYEAR

# set working directory for intermediate files 
workDIR=$VSC_DATA_VO_USER/Thesis/ISIMIP3a  

echo **workDIR is $workDIR

#==============================================================================
# initialisation
#==============================================================================

# set experiment name and experiment years based on flag_exper
exper=('empty' 'obsclim' 'counterclim')

# set variable names
VARs=('pr')


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

mygrid=$workDIR/mygrid

# ==============================================================================
# make directories
# ==============================================================================

# choose experiment
for flag_exper in "${flags_exper[@]}"; do
echo **making directories for ${exper[$flag_exper]} simulations

# make and set experiment output directory and GCM directories if they don't already exist 
EXP=${exper[$flag_exper]} # e.g. hist-nat

workDIREXP=$workDIR/$EXP  
if [ ! -e $workDIREXP ]
then
	mkdir $workDIREXP
fi
echo **work directory is $workDIREXP

echo **data is
	for data in "${data_name[@]}"; do
	echo $data

	if [ ! -e $workDIREXP/$data ] 
	then
		mkdir $workDIREXP/$data
	fi

	done
done

# ==============================================================================
# processing : crop all files, merge all timesteps, reproject files
# ==============================================================================

cd $workDIREXP

# choose experiment, loop through files to crop files, merge files, reproject files
for flag_exper in "${flags_exper[@]}"; do
echo **preprocessing ${exper[$flag_exper]} simulations

	for data in "${data_name[@]}"; do    
	
	echo **cropping $data_name files
	for FILE in $inDIR/*_${VARs}_*.nc; do
		FILENAME=$(basename -s .nc ${FILE})
		cdo -O sellonlatbox,20,40,-15,10 $FILE $workDIREXP/${FILENAME}_cropped.nc
	done

	echo **merging $data_name files and reprojecting 
    cdo -O -L -remapcon2,${mygrid} -sellonlatbox,27.5,38.5,-6,3 -mergetime $workDIREXP/*_cropped.nc $workDIREXP/${data}/${data}_${flag_exper}_pr_owngrid_${startYEAR}_${endYEAR}.nc # change name so exper is in the final filename

	find $workDIREXP/ -type f -name '*_cropped.nc' -exec mv {} ${VSC_DATA_VO_USER}/trash/ \;

	done
done



#==============================================================================
# end
#==============================================================================

# toc
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo
