#!/bin/bash

#this script assumes that the Pfiles for the forecast data already exists in the current rundirectory
# this script creates 2day backward trajectories initialized at 48h into the forecast (ie duration of -48h)

#Arguments for this script to work are following:
#forecast date (given as 'YYYYmmdd_hh')
#startdate for the bw trajectories (given as 'YYYYmmdd_hh') (here +48h)
#start longitude, eg. 5 or 5.63 or -4
#start latitude, eg. 84.53
#E.g., ./runOden_48h.sh 20200417_00 20200419_00 -5 15
#note that the lat and lon positions if given in decimal degrees are separated with a "dot" and not "comma"! 

# -------------------------------------------------------
# time references
# -------------------------------------------------------
forecastdate=$1
startdate=$2
dur_h=-48

#get the end date
enddate=`newtime ${startdate} ${dur_h}`
echo ${enddate}

# -------------------------------------------------------
# Trajectory calculations
# -------------------------------------------------------
# create startfiles
lon=$3
lat=$4
#Radius of circle:180km, distance between startpoint within the circle: 50km
STARTFILE=startf_ARTofMELT_forecast_${forecastdate}_fromOden

if [ -f "$STARTFILE" ]; then
    echo "$STARTFILE exist"
else
    echo "$STARTFILE does not exist"

    echo "lets create the startf file"

    create_startf ${startdate} ${STARTFILE} 'circle.eqd('${lon}','${lat}',180,50) @profile(50,300,6) @hPa,agl'

   echo "Startfile at $startdate done."
fi

# compute trajectories
TRAJFILE=traj_ARTofMELT_forecast_${forecastdate}_fromOden_48h_2dbw
if [ -f "$TRAJFILE" ]; then
    echo "$TRAJFILE exist"
else
    echo "$TRAJFILE does not exist"

    echo "lets create the traj file"

    #hourly output
    caltra ${startdate} ${enddate} ${STARTFILE} ${TRAJFILE} -j -o 60
    echo "Backward Trajectory calculations from Oden on $startdate done."
fi
#trace T and q
trace ${TRAJFILE} ${TRAJFILE}


