#!/bin/bash

#this script assumes that the Pfiles for the forecast data already exists in the current rundirectory
# this script creates 4day forward trajectories initialized at the date of the forecast along 70N for 120 equally spaced dots (delta lon = 3deg)

#Arguments for this script to work are following:
#forecast date (given as 'YYYYmmdd_hh')
#(here startdate for the fw trajectories (given as 'YYYYmmdd_hh') at +0h)
#E.g., ./runLat70N_0h.sh 20200417_00
#note that the lat and lon positions if given in decimal degrees are separated with a "dot" and not "comma"! 

# -------------------------------------------------------
# time references
# -------------------------------------------------------
forecastdate=$1
startdate=${forecastdate}
dur_h=96

#get the end date
enddate=`newtime ${startdate} ${dur_h}`
echo ${startdate} ${enddate}

# -------------------------------------------------------
# Trajectory calculations
# -------------------------------------------------------
# create startfiles
#name of startfile
STARTFILE=startf_ARTofMELT_forecast_${forecastdate}_70Nlatband

#check if file exists
if [ -f "$STARTFILE" ]; then
    echo "$STARTFILE exist"
else
    echo "$STARTFILE does not exist"

    echo "lets create the startf file"
    #120 points from along 70N, equally distributed
    create_startf ${startdate} ${STARTFILE} 'line(-180,177,70,70,120) @profile(50,300,6) @hPa,agl'

   echo "Startfile at $startdate done."
fi

# compute trajectories
#name of trajectory file
TRAJFILE=traj_ARTofMELT_forecast_${forecastdate}_70Nlatband_0h_4dfw

# check if file exists
if [ -f "$TRAJFILE" ]; then
    echo "$TRAJFILE exist"
else
    echo "$TRAJFILE does not exist"

    echo "lets create the traj file"

    #hourly output
    caltra ${startdate} ${enddate} ${STARTFILE} ${TRAJFILE} -j -o 60
    echo "Forward Trajectory calculations from 70N on $startdate done."
fi
# trace T, Q and TH along the trajectories
trace ${TRAJFILE} ${TRAJFILE}


