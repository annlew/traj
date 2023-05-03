#!/bin/bash


#SBATCH -n 1
#SBATCH -J 'oden'
#SBATCH -t 1:00:00


#####################################################################
# forecast initial time (or basetime, BT)
module load python3
module load rclone


# Some forecast date for testing
#BT="20230502_00"

# BT from ecaccess

#MSJ_BASETIME="00"
#MSJ_STEP="144"
#MSJ_YEAR="2023"
#MSJ_MONTH="05"
#MSJ_DAY="02"
#MSJ_EXPVER="0001"
#MSJ_EVENT="fc00h144"

# Fetch ECaccess event variables for foreacts date
BT=${MSJ_YEAR}${MSJ_MONTH}${MSJ_DAY}_${MSJ_BASETIME}

echo $BT


traj_home=/home/sr9/artofmelt/TrajectoryScripts/ecaccess_trajectories
echo $traj_home

# forecast lead times
step_start=0
step_end=240
step_inc=6 # increment

# Forecast steps for mars command (needed in special format)
steps_mars=$step_start/to/$step_end/by/$step_inc

# Forecast steps for data renaming
steps=($(seq $step_start $step_inc $step_end))

# time interval (in h) for how long the trajectories should be computed: e.g. 48h
trajectory_duration=48

# area for data retrieval
lonmin="-180"
lonmax="180"
latmin="0"
latmax="90"
area=${latmax}/${lonmin}/${latmin}/${lonmax} 

# resolution for data retrieval
res=1.0

# startfile command
#for help, just type create_startf into the command line
# this determines the starting points of the trajectories (e.g. in a rect. box with equidistant grid spacing of 100 km and at pressure levels between 1000 and 700 hpa)
# MODIFY THIS
# not needed
#startfile_command="box.eqd(-130,80,15,60,100) @ profile(1000,700,7) @ hPa"

# selection criterion
# filter the computed trajectories by criteria
# e.g. WCB-crieterion: ascent of at least 600 hPa in 48 hours
# not needed
#selection_criterion="SPECIAL:WCB:600,0,48"

### Define directories
# Base directory
#basedir=/perm/nemp/testdata/ec.oper
basedir=/scratch/sr9/testdata/ec.oper
#lagrantodir="/path/to/varf_api" # for varf_api
lagrantodir="/home/sr9/lagranto_oden" # for varf_api



# Directory where data is downloaded to
datadir_grb=${basedir}/${BT}/grb
if [ ! -d ${datadir_grb} ] ; then mkdir -p ${datadir_grb} ; fi


### Data retrieval
# model level
mars <<EOF
RETRIEVE,
date=${BT:0:8},
time=${BT:9:2},
class=od,
expver=1,
type=fc,
stream=oper,
param=130/131/132/133/135,
levtype=ml,
levelist=46/to/137,
area=$area,
grid=$res/$res,
step=$steps_mars,
target="$datadir_grb/hres_ml_${BT:0:8}${BT:9:2}_[step]"

EOF

# surface level
mars <<EOF
RETRIEVE,
date=${BT:0:8},
time=${BT:9:2},
class=od,
expver=1,
type=fc,
stream=oper,
param=134,
levtype=sfc,
area=$area,
grid=$res/$res,
step=$steps_mars,
target="$datadir_grb/hres_sfc_${BT:0:8}${BT:9:2}_[step]"

EOF

# Rename data (including leading zeres)
for step in ${steps[@]};
do	step_long=$(printf "%03d" ${step})
	mv $datadir_grb/hres_ml_${BT:0:8}${BT:9:2}_${step} $datadir_grb/hres_ml_${BT:0:8}${BT:9:2}_${step_long}
	mv $datadir_grb/hres_sfc_${BT:0:8}${BT:9:2}_${step} $datadir_grb/hres_sfc_${BT:0:8}${BT:9:2}_${step_long}
done


### Data preprocessing
echo "Calculating lagranto $step_inc hourly until lead time $step_end"

# Depending on how long the trajectories are computed for (e.g. 48 hour forward trajectories), the max. time step at which trajectories can be started is the max lead time minus the trajectory time period

# create loop variables for forecast lead time
let step_end_lagranto=step_end-trajectory_duration
steps=($(seq -w ${step_start} ${step_inc} ${step_end}))
steps_lagranto=($(seq -w ${step_start} ${step_inc} ${step_end_lagranto}))

n_steps=${#steps_lagranto[@]}

echo "###"
echo "Data preprocessing"
datadir_cdf="${basedir}/${BT}/cdf"
if [ ! -e $datadir_cdf ]; then mkdir -p $datadir_cdf; fi

## MOVE SONJA'S stuff

echo ${traj_home}
cp ${traj_home}/runOden_48h.sh ${datadir_cdf}
cp ${traj_home}/runOden_96h.sh ${datadir_cdf}
cp ${traj_home}/runLat70N_0h.sh ${datadir_cdf}
cp ${traj_home}/runLat70N_48h.sh ${datadir_cdf}
cp ${traj_home}/plot_Oden_48h_forecastTrajs_ARTofMELT.py ${datadir_cdf}
cp ${traj_home}/plot_Oden_96h_forecastTrajs_ARTofMELT.py ${datadir_cdf}
cp ${traj_home}/plot_70Nlatband_0h_forecastTrajs_ARTofMELT.py ${datadir_cdf}
cp ${traj_home}/plot_70Nlatband_48h_forecastTrajs_ARTofMELT.py ${datadir_cdf}
cp ${traj_home}/tracevars ${datadir_cdf} 

cd ${datadir_cdf}



# Fetch Oden location
if [ ! -e Odenloc.txt ]; then

echo 'Fetch new Odenloc.txt' 
ftp bolftp.ecmwf.int<< EOF
cd artofmelt
get Odenloc.txt
quit
EOF

else
  echo 'File Odenloc.txt already exists' 
fi



# loop over forecast lead times
for step in ${steps[@]}
do
	# vt: valid time
	vt=$(newtime ${BT} ${step})

	if [ ! -e P$vt ]; then
		# Model level file
		grb2=hres_ml_${BT:0:8}${BT:9:2}_${step}
		# surface file
		grb1=hres_sfc_${BT:0:8}${BT:9:2}_${step}
		
		# Link grib data files to directory where lagranto will be exectued (no long path names!).
		if [ ! -f ${datadir_cdf}/$grb2 ]
		then
			ln -s ${datadir_grb}/${grb2} ${datadir_cdf}/.
		fi

                if [ ! -f ${datadir_cdf}/$grb1 ]
                then
                        ln -s ${datadir_grb}/${grb1} ${datadir_cdf}/.
                fi

		echo "*** Convert grib to netcdf in ive format"
	        fgrb2cdf_api -c ml_cst  -v $lagrantodir/varf_api -3 U V OMEGA T Q -2 PS -o P$vt -g $grb2 $grb1

	else
        	echo 'File 'P$vt' already exists' 
	fi
done
rm ${datadir_cdf}/hres_ml* ${datadir_cdf}/hres_sfc*


# prepare trajectory computation


IFS=' ' read clon clat < <(sed 2p Odenloc.txt)
lon=$clon
lat=$clat
echo $lon
echo $lat

${datadir_cdf}/runOden_48h.sh $BT $lon $lat
${datadir_cdf}/runOden_96h.sh $BT $lon $lat
${datadir_cdf}/runLat70N_0h.sh $BT $lon $lat
${datadir_cdf}/runLat70N_48h.sh $BT $lon $lat
python3 ${datadir_cdf}/plot_Oden_48h_forecastTrajs_ARTofMELT.py ${MSJ_YEAR} ${MSJ_MONTH} ${MSJ_DAY}
python3 ${datadir_cdf}/plot_Oden_96h_forecastTrajs_ARTofMELT.py ${MSJ_YEAR} ${MSJ_MONTH} ${MSJ_DAY}
python3 ${datadir_cdf}/plot_70Nlatband_0h_forecastTrajs_ARTofMELT.py ${MSJ_YEAR} ${MSJ_MONTH} ${MSJ_DAY}
python3 ${datadir_cdf}/plot_70Nlatband_48h_forecastTrajs_ARTofMELT.py ${MSJ_YEAR} ${MSJ_MONTH} ${MSJ_DAY}

if [ ! -e plots ]; then mkdir -p plots; fi
mv PLOT* plots
mv traj* plots
cp Odenlox.txt plots

rclone copy ${datadir_cdf}/plots box:artofmelt/${BT}/trajectory


echo '###############################################################' 
echo  'LAGRANTO  done' 
echo '###############################################################' 





