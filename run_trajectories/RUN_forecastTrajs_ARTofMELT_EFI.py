#!/usr/bin/python
#!/usr/bin/env python

#import modules
import subprocess
import shlex
import sys,glob,os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time,datetime
import string
from pathlib import Path

# The idea of this function is to create a startfile for 4d forward trajectories initialized at +Xh wrt forecast date initialization
# startlocations are based on a circle with lon/lat given by the user (lon, lat in decimal degrees, separated by "dot", -180 to 180 lon and -90 to 90 lat
#Running this script requires that the P files and the'Odenloc.txt' - file is updated and done
#define subprocess functions

def calcEFI(filename,forecastdate='20200417_00'):
    '''
    This function is only for calculating the trajectories. P files linked to the same folder as the rundirectory!
    Remember to check that the directory for the files is correct
    filename=runEFI_Xh.sh
    Input:
    forecastdate=string of the forecast initialization date.
    The script asks the user to give lon, lat and time reference (in h) for starting location of forward trajs:
    location for the origo of the circle given manually by the user (lon0, lat0), and hours (how many hours into the forecast is wanted to use as starting time of trajectories (if 0, then the forecast date is used as initialization time for the trajs), if 24, the trajectories are initialized at +24h into the forecast)
    Variabels along the trajectory will be traced for T and q
    '''
    subprocess.call(shlex.split('./%s %s'%(filename,forecastdate)))

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


if __name__ == "__main__":
   
    ########################Arguments for the run##################

    Y=int(sys.argv[1]) #same year 2023
    M=int(sys.argv[2]) #month of forecast initialisation
    d=int(sys.argv[3]) # day of forecast initialisation
    h=00 #forecast initialized at 00 UTC default time

    ####################### Define start dates ###################
    Date_forecast=str(Y) + str(M).zfill(2) + str(d).zfill(2) + '_' + str(h).zfill(2)
    print(Date_forecast)

    ############################# Calculate trajectories with LAGRANTO ###################
    
    # LINK P files to the run directory for trajectory computations
    linkPfiles()
    print('P files linked to the run directory for the %s forecast'%Date_forecast)    
    # 4d forward from a manually given location, initialized from 6 vertical levels 50 - 200 hPa AGL ,dP = 50hPa, with a circle of 180km radius, start. points eq 50km
    calcEFI('runEFI_Xh.sh',Date_forecast)
    print('calculation done for 4-day forward trajectories initialized at given location')

    # remove the linked P files from the current run directory
    remlinkPfiles()

    ########################### PLOT trajectories ####################
    
    #### current location of Oden ########
    loc=pd.read_csv('Odenloc.txt', sep=' ',header=None)
    loc.columns=['lon','lat']
    
