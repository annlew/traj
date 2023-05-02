#!/bin/bash

#this script assumes that the Pfiles for the forecast data already exists in the current rundirectory
# this script creates 4day backward trajectories initialized at 96h into the forecast (ie duration of -96h)

#Arguments for this script to work are following:
#forecast date (given as 'YYYYmmdd_hh')
#start longitude, eg. 5 or 5.63 or -4
#start latitude, eg. 84.53
#E.g., ./runOden_96h.sh 20200417_00 -5 15
#note that the lat and lon positions if given in decimal degrees are separated with a "dot" and not "comma"! 

# -------------------------------------------------------
# time references
# -------------------------------------------------------
forecastdate=$1
startdate=`newtime ${forecastdate} 96` #start 96 hours into the forecast
dur_h=-96

#get the end date
enddate=`newtime ${startdate} ${dur_h}`
echo ${startdate} ${enddate}

# -------------------------------------------------------
# Trajectory calculations
# -------------------------------------------------------
# create startfiles
#Oden location (read from the textfile Odenloc.txt)
lon=$2
lat=$3
# name of startfile
STARTFILE=startf_ARTofMELT_forecast_${forecastdate}_fromOden

# check if file exists
if [ -f "$STARTFILE" ]; then
    echo "$STARTFILE exist"
else
    echo "$STARTFILE does not exist"
    # create file if it does not exist
    echo "lets create the startf file"
    # Radius of circle:180km, distance between startpoint within the circle: 50km, origo of circle = Oden
    create_startf ${startdate} ${STARTFILE} 'circle.eqd('${lon}','${lat}',180,50) @profile(50,300,6) @hPa,agl'

   echo "Startfile at $startdate done."
fi

# compute trajectories
# name of trajectory file
TRAJFILE=traj_ARTofMELT_forecast_${forecastdate}_fromOden_96h_4dbw

# check if file exists
if [ -f "$TRAJFILE" ]; then
    echo "$TRAJFILE exist"
else
    echo "$TRAJFILE does not exist"

    echo "lets create the traj file"

    # hourly output
    caltra ${startdate} ${enddate} ${STARTFILE} ${TRAJFILE} -j -o 60
    echo "Backward Trajectory calculations from Oden on $startdate done."
fi
# trace T, Q and TH along the trajectories
trace ${TRAJFILE} ${TRAJFILE}


