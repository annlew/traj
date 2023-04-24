#!/usr/bin/python
#!/usr/bin/env python

#import modules
import subprocess
import shlex
import sys,glob,os
import numpy as np
import matplotlib.pyplot as plt
import time,datetime
import string
import pandas as pd
from pathlib import Path

#NOTE! This function pre-requirements are that P-files are created in a specific folder and the Odenloc.txt file is updated/replaced with the new location!
#define subprocess functions

def calcOden48h(filename,forecastdate='20200417_00',startdate='20200419_00',lon0=5,lat0=84):
    '''
    This function is only for calculating the trajectories. P files required in the same folder as the rundirectory!
    Remember to check that the folder directories in the files used are correct!
    filename=runOden_48h.sh
    forecastdate= string of the forecast initialization date
    startdate = string of the startdate for trajectories
    lon0, lat0 = location of Oden extracted from the sounding BUFR file at 00 UTC
    Duration (in hours) for trajectory computations. Observe that "-" means backward trajectories, eg. -48 for -2d bw
    Variabels along the trajectory will be traced for T and q
    '''
    subprocess.call(shlex.split('./%s %s %s %i %i'%(filename,forecastdate,startdate,lon0,lat0)))

def calcOden96h(filename,forecastdate='20200417_00',startdate='20200419_00',lon0=5,lat0=84):
    '''
    This function is only for calculating the trajectories. P files required in the same folder as the rundirectory!
    Remember to check that the folder directories in the files used are correct!
    filename=runOden_96h.sh
    forecastdate= string of the forecast initialization date
    startdate = string of the startdate for trajectories
    lon0, lat0 = location of Oden extracted from the sounding BUFR file at 00 UTC
    Duration (in hours) for trajectory computations. Observe that "-" means backward trajectories, eg. -96 for -4d bw
    Variabels along the trajectory will be traced for T and q
    '''
    subprocess.call(shlex.split('./%s %s %s %i %i'%(filename,forecastdate,startdate,lon0,lat0)))

#lat band, plotting functions also included here
# names: ARTofMELT_forecast_YYYYmmdd_00_fromOden_48h_2dbw, ARTofMELT_forecast_YYYYmmdd_00_fromOden_96h_4dbw, ARTofMELT_forecast_YYYYmmdd_00_70Nlatband_4dfw
#ARTofMELT_forecast_YYYYmmdd_00_EFI_4dfw

def linkPfiles():
    '''
    This function links the forecast Pfiles from a directory defined in the bash script to the current runninf directory 
    '''
    subprocess.call('./linkPfiles.sh')


def remlinkPfiles():
    '''
    This function removes the link to the forecast Pfiles from the current rundirectory  
    '''
    subprocess.call('./removelinkedPfiles.sh')


#define python funtions

def strtodatetime_newtime(datestr='20230508_00',format_st='%Y%m%d_%H',hours_dif=-1):
    '''
    Input: 
    datestr : date of forecast at 00 UTC, as a string: 'YYYYmmdd_hh', where hh = 00
    format_str : format of the given input date
    hours_dif = hours plus/ minus (-) compared to the forecast initialisation time
    
    Output: date +/-X hours, given in same format: 'YYYYmmdd_hh',
    eg. for a sounding launched at 00UTC, we want the sounding file by name one hour earlier, eg. hours_dif= -1, so that hh = 23 and dd = dd(in) - 1day,
    '''
    Date_forecast=datetime.datetime.strptime(datestr, format_st)
    Date_new=Date_forecast + datetime.timedelta(hours=hours_dif)
    Date_tuple=Date_new.year,Date_new.month,Date_new.day,Date_new.hour
    Strdate_new=str(Date_tuple[0])+"{:0>2}".format(str(Date_tuple[1]))+"{:0>2}".format(str(Date_tuple[2]))+'_'+"{:0>2}".format(str(Date_tuple[3]))
    return(Strdate_new)


if __name__ == "__main__":
   
    ########################Arguments for the run##################

    Y=int(sys.argv[1]) #same year 2023
    M=int(sys.argv[2]) #month of forecast initialisation
    d=int(sys.argv[3]) # day of forecast initialisation
    h=00 #forecast initialized at 00 UTC default time

    ####################### Define start dates ###################
    Date_forecast=str(Y) + str(M).zfill(2) + str(d).zfill(2) + '_' + str(h).zfill(2)
    Date_48h=strtodatetime_newtime(Date_forecast, hours_dif=48)
    Date_96h=strtodatetime_newtime(Date_forecast, hours_dif=96)
    print(Date_forecast, Date_48h,Date_96h)

    ##################### Read location of Oden from text file extracted from sounding Bufr files #####
    
    # Start target location of Oden (default location):
    #lon=0.
    #lat=84.
    #Read a NEW/updated file every day!
    loc=pd.read_csv('Odenloc.txt', sep=' ',header=None)
    loc.columns=['lon','lat']
    lon0=loc['lon'][0]
    lat0=loc['lat'][0]
    print(lon0, lat0)
    
    ############################# Calculate trajectories with LAGRANTO ###################
    
    # LINK P files to the run directory for trajectory computations
    linkPfiles()
    print('P files linked to the run directory for the %s forecast'%Date_forecast)    
    # backward from Oden, initialized from 6 vertical levels 50 - 200 hPa AGL ,dP = 50hPa, with a circle of 180km radius, start. points eq 50km
    calcOden48h('runOden_48h.sh',Date_forecast,Date_48h,lon0,lat0)
    calcOden96h('runOden_96h.sh',Date_forecast,Date_96h,lon0,lat0)
    print('calculation done for trajectories initialized at Oden (%s,%s)'%(str(lon0),str(lat0)))

    #calcLatband70N()
    #print('calculation done for trajecectories initialized from the 70N latitude band')

    # remove the linked P files from the current run directory
    remlinkPfiles()

    ########################### PLOT trajectories ####################    
