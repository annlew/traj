#!/bin/bash

#this script assumes that the Pfiles for the forecast data already exists in the current rundirectory
# this script creates 4day forward trajectories initialized at Xh wrt time of forecast initialization (duration of trajectories is 96h) started at a user specified location

#Arguments for this script to work is following:
#forecast date (given as 'YYYYmmdd_hh')
#E.g., ./runEFI_0h.sh 20200417_00
#note that the lat and lon positions if given in decimal degrees are separated with a "dot" and not "comma"!

# -------------------------------------------------------
# User specified lon,lat location
# -------------------------------------------------------

# ask the user for lon and lat position
echo "give longitude (in decimal degrees, separated by dot: e.g. 5.2)"
read lon
lon0=${lon}
echo "give latitude (in decimal degrees, separated by dot: e.g. 87.4)"
read lat
lat0=${lat}


#save the lon,lat as integers for filenames
int_lon0=$( LC_NUMERIC="en_US.UTF-8" printf "%.0f" $lon0 )
int_lat0=$( LC_NUMERIC="en_US.UTF-8" printf "%.0f" $lat0 )


# -------------------------------------------------------
# time references
# -------------------------------------------------------

# ask the user for hours into the forecast when the calculations should start (specifying the start time. If 0, then startdate = forecastdate)
echo "how many hours into the forecast shall the trajectories start? (give hours, eg 0 or 24; max 144 (equals 6d))"
read hours

echo "Forward trajectories are computed from a circle with origo at ($lon,$lat), started at +$hours into the forecast"

hours0=${hours}
forecastdate=$1
startdate=`newtime ${forecastdate} ${hours0}`
echo ${startdate}
dur_h=96 #4d forward

#get the end date
enddate=`newtime ${startdate} ${dur_h}`
echo ${enddate}

# -------------------------------------------------------
# Trajectory calculations
# -------------------------------------------------------


#save chosen startlocation in a file (append if many times chosen)
#format: forecastdate hours(into forecast for trajstart) lon(org) lat(org) lon(rounded) lat(rounded (per row)
file=EFIloc.txt
printf '%s\n' "${forecastdate} ${hours0} $lon $lat ${int_lon0} ${int_lat0}" >> "$file"

# create startfiles
#Radius of circle:180km, distance between startpoint within the circle: 50km
STARTFILE=startf_ARTofMELT_forecast_${forecastdate}_${int_lon0}_${int_lat0}_EFI

if [ -f "$STARTFILE" ]; then
    echo "$STARTFILE exist"
else
    echo "$STARTFILE does not exist"

    echo "lets create the startf file"

    create_startf ${startdate} ${STARTFILE} 'circle.eqd('${lon0}','${lat0}',180,50) @profile(50,300,6) @hPa,agl'

   echo "Startfile at $startdate done."
fi

# compute trajectories
TRAJFILE=traj_ARTofMELT_forecast_${forecastdate}_EFI_${int_lon0}_${int_lat0}_${hours0}h_4dfw
if [ -f "$TRAJFILE" ]; then
    echo "$TRAJFILE exist"
else
    echo "$TRAJFILE does not exist"

    echo "lets create the traj file"

    #hourly output
    caltra ${startdate} ${enddate} ${STARTFILE} ${TRAJFILE} -j -o 60
    echo "Forward trajectory calculations from $lon0 , $lat0 on $startdate done."
fi
#trace T and q
trace ${TRAJFILE} ${TRAJFILE}


