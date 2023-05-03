#!/bin/bash

# run:
# ./run_trajectory.sh year month day
#
# Ex:
# ./run_trajectory.sh 2023 05 03

# Path to P-files/foreacst directory
Pfile_dir=/scratch/sr9/testdata/ec.oper/${1}${2}${3}_00/cdf

# Check if forecast data exists
if [ ! -e $Pfile_dir ]; then echo 'Forecast not ready for '${1} ${2} ${3} ; exit; fi

# Link input file and Oden location at forecast time 
ln -s ${Pfile_dir}/P* .
ln -s ${Pfile_dir}/ml_cst .
ln -s ${Pfile_dir}/Odenloc.txt .

# Run lagranto and plot
python3 RUN_forecastTrajs_ARTofMELT_EFI_combi.py ${1} ${2} ${3} 

# Move to output directory on scratch
out_dir=${SCRATCH}/artofmelt/${1}${2}${3}_00
echo 'Move plots and trajectories to '$out_dir
if [ ! -e ${out_dir} ]; then mkdir -p ${out_dir}; fi
mv PLOT* ${out_dir}
mv traj* ${out_dir}

rm start*

#Copy to box
rclone copy --ignore-existing ${out_dir} box:artofmelt/${1}${2}${3}_00/trajectory_extras
rclone copyto EFIloc.txt box:artofmelt/Extra_trajectory_${USER}.txt 

#remove linked files
rm P*
rm Odenloc.txt
rm ml_cst


