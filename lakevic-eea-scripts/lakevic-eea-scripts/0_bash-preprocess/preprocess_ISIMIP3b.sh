#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --ntasks=3
#SBATCH --mail-type=ALL

#==============================================================================
# PREPROCESS ISIMIP3b data 
# crops daily ISIMIP3b precipitation data to AGL, merges through time (1850-2020) and reprojects to own grid
# Created: 6 April 2022. Modified: 11 May 2022. 
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
              # 1: hist
              # 2: ssp3-rcp70
              # 3: hist-nat
flags_data=2
	      # 1: primary input data
	      # 2: secondary input data


# set input directories and start and end year of simulations
if [ $flags_exper -eq 1 ]
then 
	startYEAR=1850
	endYEAR=2014
	if [ $flags_data -eq 1 ]
	then
		inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily/historical/
	elif [ $flags_data -eq 2 ]
	then
		inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/SecondaryInputData/climate/atmosphere/bias-adjusted/global/daily/historical/
	fi

elif [ $flags_exper -eq 2 ]
then
	startYEAR=2015
	endYEAR=2020
	if [ $flags_data -eq 1 ]
	then 
		inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily/ssp370/
	elif [ $flags_data -eq 2 ]
	then
		inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/SecondaryInputData/climate/atmosphere/bias-adjusted/global/daily/ssp370/
	fi

elif [ $flags_exper -eq 3 ]
then
	inDIR=/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/SecondaryInputData/climate/atmosphere/bias-adjusted/global/daily/hist-nat/
	startYEAR=1850
	endYEAR=2020
fi

echo **inDIR is $inDIR, period is $startYEAR to $endYEAR

# set output directory
outDIR=/data/brussel/104/vsc10419/Thesis/ISIMIP3b

# set working directory for intermediate files 
workDIR=$VSC_DATA_VO_USER/Thesis/ISIMIP3b

echo **workDIR is $workDIR

#==============================================================================
# initialisation
#==============================================================================

# set experiment name and experiment years based on flag_exper
exper=('empty' 'hist' 'rcp70' 'hist-nat')

# set variable names
VARs=('pr')

# set GCM names
if [ $flags_exper -eq 3 ]
then
	GCMs=('CanESM5' 'CNRM-CM6-1' 'GFDL-ESM4' 'IPSL-CM6A-LR' 'MIROC6' 'MRI-ESM2-0') # hist-nat
else
	GCMs=('CanESM5' 'CNRM-CM6-1' 'MIROC6') # GCMs=('GFDL-ESM4' 'IPSL-CM6A-LR' 'MRI-ESM2-0') # hist and rcp70 (1) Secondary Input Data (2) Input Data
fi

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

mygrid=$outDIR/mygrid

# ==============================================================================
# make directories
# ==============================================================================

# choose experiment
for flag_exper in "${flags_exper[@]}"; do
echo **making directories for ${exper[$flag_exper]} simulations

# make and set experiment output directory and GCM directories if they don't already exist 
EXP=${exper[$flag_exper]} # e.g. hist-nat
outDIREXP=$outDIR/$EXP
if [ ! -e $outDIREXP ]
then
	mkdir $outDIREXP
fi
echo **output directory is $outDIREXP

workDIREXP=$workDIR/$EXP
if [ ! -e $workDIREXP ]
then
	mkdir $workDIREXP
fi
echo **work directory is $workDIREXP

echo **GCMs are
	for GCM in "${GCMs[@]}"; do
	echo $GCM
	if [ ! -e $outDIREXP/$GCM ]
	then
		mkdir $outDIREXP/$GCM
	fi

	if [ ! -e $workDIREXP/$GCM ]
	then
		mkdir $workDIREXP/$GCM
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

	for GCM in "${GCMs[@]}"; do    #!! REPLACE 0 WITH @ TO RUN ALL GCMs and vice versa for test !!!
	
	echo **cropping $GCM files
	for FILE in $inDIR/$GCM/*_${VARs}_*.nc; do
		FILENAME=$(basename -s .nc ${FILE})
		#echo $FILENAME
		#echo $workDIREXP/$GCM/${FILENAME}_cropped.nc  
		cdo -O sellonlatbox,20,40,-15,10 $FILE $workDIREXP/$GCM/${FILENAME}_cropped.nc
	done


	echo **merging $GCM files
	if [ $flags_exper -eq 2 ] #rcp70 only take first file, up to end 2020
	then
	cp $workDIREXP/$GCM/*${startYEAR}_${endYEAR}_cropped.nc $workDIREXP/$GCM/${GCM}_pr_AGL_${startYEAR}_${endYEAR}.nc 
	else
	cdo -O mergetime $workDIREXP/$GCM/*_cropped.nc $workDIREXP/$GCM/${GCM}_pr_AGL_${startYEAR}_${endYEAR}.nc 
	fi

	echo **reprojecting $GCM files
	cdo -O -L -remapcon2,${mygrid} -sellonlatbox,27.5,38.5,-6,3 $workDIREXP/$GCM/${GCM}_pr_AGL_${startYEAR}_${endYEAR}.nc $workDIREXP/$GCM/${GCM}_pr_owngrid_${startYEAR}_${endYEAR}.nc

#	echo **copying $GCM files to outDIR
#	cp $workDIREXP/$GCM/${GCM}_pr_owngrid_${startYEAR}_${endYEAR}.nc $outDIREXP/$GCM/${GCM}_pr_owngrid_${startYEAR}_${endYEAR}.nc

	find $workDIREXP/$GCM/ -type f -name '*_cropped.nc' -exec mv {} $workDIR/trash-intermediatefiles-2/ \;

	done
done



#==============================================================================
# end
#==============================================================================

# toc
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo
