#!/usr/bin/python
#!/usr/bin/env python

#import modules
import sys,glob,os
import numpy as np
import math
import cartopy.crs as ccrs
import pandas as pd
import matplotlib.dates as mdates
from cartopy.geodesic import Geodesic
import shapely.geometry as sgeom
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import matplotlib.path as mpath
import xarray as xr
import matplotlib.colors as mcolors #new
import matplotlib.gridspec as gridspec
import time,datetime
import string
#import cftime
import pandas as pd
from pathlib import Path

#define python funtions
#load modules
def getonetraj_multipleP(nSteps,nTrajs,df,plevs=None):
    '''
    Returns a dictionary of the trajectories for one initialization time of all trajectories.
    if plev=None: Returns a dictionary with trajectories for each location of a a certain block at ref time
    Assign plevs = None if there is only one initial pressure level chosen
    Assign plevs as a list of pressure levels if you want to have a dictionary of Trajs[plev] or Trajs[trajindex]
    Trajs[trajindex] - each trajectory as one part of the dictionary, max(trajindex) = nTrajs - 1
    plevs can also be given as hPa AGL, such as [50,100,150,200,250,300]
    nSteps = amount of steps (1-hourly) in one trajectory, i.e. if 4days trajs, then we have nSteps=4*24 = 96, if -2 days, then -2 * 24 = 48
    nTrajs = amount of trajectories calculated for a block area (over all pressure levels if you have P-levels)
    df = dataframe with all trajectory information (read in from a text file - LAGRANTO output)
    '''
    if plevs is None:
        #we have only one pressure level - this one might not work...
        size=int(nSteps +1) #size for one trajectory
        list_of_Trajs = [df.loc[i:i+size-1,:] for i in range(0, len(df),size)] #list of trajectories
        if len(list_of_Trajs) == int(nTrajs):
            Trajs={}
            for j,traj in enumerate(list_of_Trajs):
                Trajs[j]=traj #create the dictionary of trajectories
        else:
            raise ValueError('Not matching lenght!')
    else:
        #we have more pressure levels as initial level
        #print(len(df)/(len(plevs)))
        size = int(nTrajs/(len(plevs))*(nSteps +1)) #size for trajectories in one pressure level
        #print(size)
        if size != (len(df)/(len(plevs))):
               raise ValueError('Not matching lenght!')
        else:
            list_of_dfs = [df.loc[i:i+size-1,:] for i in range(0, len(df),size)] #divide into P levels
            if len(list_of_dfs) != len(plevs):
                raise ValueError('Not matching amount of given Pressure levels!')
            else:
                Trajs={} #per plev
                for p,df_p in enumerate(list_of_dfs):
                    if len(df_p)/(nSteps+1) == 1:
                        print('we have only one trajectory??')
                        Trajs[plevs[p]]=df_p.reset_index(drop=True) #dataframe per pressure level, only one traj
                    else:
                        size_p=int(nSteps +1) #size for one trajectory
                        list_of_Trajs = [df_p.reset_index(drop=True).loc[i:i+size_p-1,:] for i in range(0, len(df_p),size_p)]
                        #print(len(list_of_Trajs)) #should be 11
                        Trajs_01={}
                        for p1,df_p1 in enumerate(list_of_Trajs):
                            Trajs_01[p1]=df_p1.reset_index(drop=True)# dictionary with index 0-10 with df for trajectories in one p
                    Trajs[plevs[p]]=Trajs_01
    return Trajs


def getList(dict):
    '''
    This function returns a list of keys for a pandas dictionary
    '''
    list = []
    for key in dict.keys():
        list.append(key)
    return list

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100,returncolors=False):
    '''
    This function creates a new colormap from the cmap chosen, taking the limits from minval and maxval, and the number of colors as n.
    The default is a colormap as LinearSegmentedColormap, which can be used as colormap for imshow, contourf...
    If you assign returncolors as True, the there will be an array of n amount of colors from the chosen colormap. This can be used in plots when each event has a color (n=50)
    '''
    new_cmap = mcolors.LinearSegmentedColormap.from_list('new_cmap', plt.get_cmap(cmap)(np.linspace(minval, maxval, n)))
    if returncolors:
        new_cmap=plt.get_cmap(cmap)(np.linspace(minval, maxval, n))
    return new_cmap

def VirtTH_perrow(row):
    #θ(1+0.61(q/1000)), where potential temperature θ is in Kelvin and
    #specific humidity (q) in kg/kg (/1000 used here because q in g/kg)
    return row.TH*(1+(0.61*(row.Q/1000)))


def compute_VirtualTH(Traj_DICT):
    '''
    this function computes the virtual potential temperature (THv) from specific humidity
    It is not taking into consideration cloud liquid and ice water content atm
    '''
    #get pressure levels:
    plevs=getList(Traj_DICT)
    #get indexes for the different trajectories saved in dataframes
    indexes_traj=getList(Traj_DICT[plevs[0]])
    for p in plevs:
        for i in indexes_traj:
            Traj_DICT[p][i]['THv']=Traj_DICT[p][i].apply(lambda row: VirtTH_perrow(row), axis=1)
    return Traj_DICT

def plottracks(TrajDF,p,cmap=plt.get_cmap('Spectral_r'),var='p',Indx_ini=None,COL=None,markDays=False,
              colorbounds=np.linspace(500,1000,26),selshort=None,ma=None,ax=None,linesty=None,fillm=False,
              ms=11,hdif=24,lw=3,plotP=True,markEnd=False,plot24h_difcol=None,plotC=False):
    '''
    TrajDF = dictionary with trajectories, keys are pressure level and index (for startingpoints)
    p = give starting pressure level (AGL)
    If plotP = True, the trajectory is colored according to the pressure (or some other variable [var])
    along the trajectory. Colormap for the coloring given in cmap, colorbounds = levels for colorbar
    If plotC = True, the tracks are colored in one color
    If COL = None: the colors will be randomly chosen from the colormap "gist_rainbow". If not none, and a color is given, this will be the color of ALL trajectories
    Indx_ini : give a value between 0 - 38 (index for the starting point to select one specific one only)
    The selected marker style per trajectory can be manually selected by "ma",  eg 'o'
    The colors for each 24h step can also be chosen by plot24h_difcol (give as a list with length of days to mark)
    if markEnd = True, marks the end point with a black dot
    hdif=frequency for marks along the trajectory - default 24
    If markDays = True, marks every hdif (h) from the starting point with a transparent circle marker
    selshort = hours for a shorter trajectory than the full length
    lw= width of the trajectory
    ms= markersize
    fillm= True if markers shall be filled, default False
    linesty = linestyle, defaulte solid line
    '''
    if ax is None:
        ax=plt.gca()
    else:
        ax=ax
    if selshort:
        TR={}
        idxes=getList(TrajDF[p])
        for idx in idxes:
            T=TrajDF[p][idx].iloc[0:selshort+1]
            TR[idx]=T
        Trajs_plev=TR
    else:
        Trajs_plev=TrajDF[p]
        
    if linesty is None:
        lst='-'
    else:
        lst=linesty
        
    #id you want to select for a certain starting point
    if Indx_ini is None:
        #use all starting points
        Indxes=getList(Trajs_plev) # return the list of initial indexes
        Traj_sel=Trajs_plev
        
    if Indx_ini is not None:
        Trajs_sel={}
        Indxes = [Indx_ini]
        # check if the chosen index is in list
        if Indx_ini not in getList(Trajs_plev):
            raise ValueError('No matching indexes!')
        else:
            Traj_sel[Indx_ini]=Trajs_plev[Indx_ini]
    #marker type
    if ma is None:
        Ma = ['o' for i in range(len(Traj_sel))]
    else:
        Ma = [ma for i in range(len(Traj_sel))]
    
    #color of trajectories
    if COL is None:
        Colors=truncate_colormap('gist_rainbow',n=len(Indxes),returncolors=True) #discrete colors for trajectories
    else:
        if len(COL) > 1:
            if (isinstance(COL[0], list)) or (isinstance(COL[0], str)):
                Colors = COL
            else:
                Colors = [COL]
        else:
            Colors = [COL for i in range(len(Indxes))]
    # norm for colormap if plotP == True
    norm = mcolors.BoundaryNorm(boundaries=colorbounds, ncolors=256)
    for j,i in enumerate(Indxes):
        #go through all trajs initialized from same pressure level, different start points /times
        lons=Traj_sel[i]['lon'].values
        lats=Traj_sel[i]['lat'].values
        p=Traj_sel[i][var].values
        
        if plotC:
            #color tracks in one given color if COL given, else random colors, one color per track
            ax.plot(lons,lats,color=Colors[j], linewidth=lw, marker=None,linestyle=lst,
                 transform=ccrs.Geodetic(),zorder=6,alpha=0.7
                )
            if markEnd:
                
                # mark end point with black marker
                ax.plot([lons[-1],lons[-1]],[lats[-1],lats[-1]],
                     markerfacecolor='k',markersize=3, markeredgecolor='None',marker='o',transform=ccrs.Geodetic(),zorder=9,
                   )
           
        if plotP:
            points = np.array([lons, lats]).T.reshape(-1, 1, 2)
            #segments = np.concatenate([points[:-1], points[1:]], axis=1)
            segments = np.concatenate([points[:-2], points[1:-1], points[2:]], axis=1) #smoother
            #linecol=p[:-1] + p[1:]
            lc = LineCollection(segments, cmap=cmap, alpha=0.6,norm=norm,transform=ccrs.Geodetic(),
                                linestyles=lst,zorder=6,linewidths=lw)
            #color of the line
            lc.set_array(p)
            #lc.set_linewidth(2)
            ax.add_collection(lc)
            if markEnd:
                # mark end point with black marker
                ax.plot([lons[-1],lons[-1]],[lats[-1],lats[-1]],
                     markerfacecolor='k',markersize=3, markeredgecolor='None',marker='o',transform=ccrs.Geodetic(),zorder=9,
                   )
        if markDays:
            #mark every 24h step from starting point
            df=Traj_sel[i] #chose the current trajectory
            timedif=sorted(list(dict.fromkeys(df['time'].diff().values))) #difference in time of one timestep [hourly]
            Dur_onetraj=df['time'].iloc[-1] #hours
            
            #print(Dur_onetraj)
            if abs(Dur_onetraj)>=hdif:
                dT = [x for x in timedif if str(x) not in ['nan', str(Dur_onetraj)]][0] # delete nan and difference between trajectories
                nSte=int(abs(hdif/dT)) #get the intevall for daily data
                
                Daily=df.iloc[::nSte, :].iloc[1:].reset_index(drop=True)#chose evert nSte:th data (daily), not lag0
                for l in range(len(Daily)):
                    
                    if plotP:
                        if plot24h_difcol is None:
                            if fillm:
                                ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                         markeredgecolor=Colors[j],markerfacecolor=Colors[j],markeredgewidth=1., markersize=ms, marker=Ma[j],
                                        transform=ccrs.Geodetic(),zorder=9,
                                       )
                            else:
                                ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                         markeredgecolor=Colors[j],markerfacecolor="None",markeredgewidth=1., markersize=ms, marker=Ma[j],
                                        transform=ccrs.Geodetic(),zorder=9,
                                       )
                        else:
                            
                            if fillm:
                                ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                 markeredgecolor=plot24h_difcol[l],markerfacecolor=plot24h_difcol[l],markeredgewidth=1., markersize=ms, marker=Ma[j],
                                transform=ccrs.Geodetic(),zorder=9,
                               )
                                
                            else:
                                ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                     markeredgecolor=plot24h_difcol[l],markerfacecolor="None",markeredgewidth=1., markersize=ms, marker=Ma[j],
                                    transform=ccrs.Geodetic(),zorder=9,
                                   )
                        #ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']],
                        #         markerfacecolor=Colors_markers[j],markersize=3, markeredgecolor='None',marker='o',transform=ccrs.Geodetic(),zorder=9,
                        #       )

                    if plotC:
                        if plot24h_difcol is not None:
                            if fillm:
                                ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                 markeredgecolor=plot24h_difcol[l],markerfacecolor=plot24h_difcol[l],markeredgewidth=1.5, markersize=ms, marker=Ma[j],
                                transform=ccrs.Geodetic(),zorder=9,
                               )
                            else:
                                ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                     markeredgecolor=plot24h_difcol[l],markerfacecolor="None",markeredgewidth=1.5, markersize=ms, marker=Ma[j],
                                    transform=ccrs.Geodetic(),zorder=9,
                                   )
                        else:
                            if COL is not None:
                                if fillm:
                                    ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                     markeredgecolor=Colors[j],markerfacecolor=Colors[j],markeredgewidth=1.5, markersize=ms, marker=Ma[j],
                                    transform=ccrs.Geodetic(),zorder=9,
                                   )
                                else:
                                    ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                         markeredgecolor=Colors[j],markerfacecolor="None",markeredgewidth=1.5, markersize=ms, marker=Ma[j],
                                        transform=ccrs.Geodetic(),zorder=9,
                                       )
                            else:
                                #black markers when trajectories obtain different colors
                                if fillm:
                                    ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                            markeredgecolor='k',markerfacecolor='k',markeredgewidth=1.5, markersize=ms, marker=Ma[j],
                                            transform=ccrs.Geodetic(),zorder=9,
                                           )
                                else:
                                    ax.plot([Daily.iloc[l]['lon'],Daily.iloc[l]['lon']],[Daily.iloc[l]['lat'],Daily.iloc[l]['lat']], alpha=0.9,
                                            markeredgecolor='k',markerfacecolor="None",markeredgewidth=1.5, markersize=ms, marker=Ma[j],
                                            transform=ccrs.Geodetic(),zorder=9,
                                           )
    
    if plotP:
        return lc
