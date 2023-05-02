#!/usr/bin/python
#!/usr/bin/env python

#import modules
import subprocess
import shlex
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

# Required for running this plot script are:
# access to trajectory file(s) and the Odenloc.txt file with the current location of Oden
# Input traj file: traj_ARTofMELT_forecast_YYYYmmdd_hh_70Nlatband_0h_4dfw
# Arguments: Year month day (of forecast initialization)

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
    nSteps = amount of steps (1-hourly) in one trajectory, i.e. if 2days trajs, then we have nSteps=2*24 = 48
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
        if var == 'T':
            p=Traj_sel[i][var].values# - 273.15 #from Kelvin to Celcius
        else:
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


if __name__ == "__main__":
   
    ######################## Arguments for the run ##################
    # same arguments as for the trajectory calculations
    Y=int(sys.argv[1]) #same year 2023
    M=int(sys.argv[2]) #month of forecast initialisation
    d=int(sys.argv[3]) # day of forecast initialisation
    h=00 #forecast initialized at 00 UTC default time

    ####################### Define parameters ###################
    Date_forecast=str(Y) + str(M).zfill(2) + str(d).zfill(2) + '_' + str(h).zfill(2)
    Date=datetime.datetime.strptime(Date_forecast, '%Y%m%d_%H') #forecast date in different format for plots
        
    #open Oden Loc in pandas
    loc=pd.read_csv('Odenloc.txt', sep=' ',header=None)
    loc.columns=['lon','lat']
    
    #initialization pressure levels in LAGRANTO for startpoints
    hPa_AGL=np.arange(50,350,50).tolist()
    
    #panel numbering
    panels = list(string.ascii_lowercase)
    
    # global settings
    plt.rc('figure', titlesize=20)
    plt.rcParams['xtick.labelsize']=20
    plt.rcParams['ytick.labelsize']=20
    plt.rcParams['legend.fontsize']=20
    plt.rcParams['axes.labelsize']=20
    plt.rcParams['figure.facecolor']='white'

    #define the folder accordingly
    TrajFile='traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast
    
    # get the trajectories (in dataframe)
    #read the selected textfile with Traj info
    df_traj=pd.read_csv(TrajFile,encoding='latin-1', delimiter=r"\s+" ,skiprows=[0,1,3],header=[0])
    Dur_h=pd.read_csv(TrajFile,encoding='latin-1', delimiter=r"\s+",nrows=1,header=None).iloc[:,6][0]/60. #for hours
    # get the steps per trajectory & amount of trajectories
    w=sorted(list(dict.fromkeys(df_traj['time'].diff().values))) #difference in time of one timestep [hourly]
    cleanedListH = [x for x in w if str(x) not in ['nan', str(Dur_h), str(abs(Dur_h)), str(Dur_h*-1)]][0] # delete nan and difference between trajectories
    nSteps_h = abs(Dur_h /cleanedListH) #how many steps for one trajectory
    nTrajs=len(df_traj)/(nSteps_h+1) #amount of trajectories in the whole file (including all from different levels)

    #CHECK that the durations are correct
    if 96 !=Dur_h:
        raise ValueError('ERROR! %s does not equal file traj duration %s [h]'%(str(96),
                                                                           str(Dur_h)))
    #dictionary for each trajectory, keys = initial pressure levels
    Trajs_Plevs=getonetraj_multipleP(nSteps_h,nTrajs,df_traj,plevs=hPa_AGL)
    
    #compute virtual potential temperature
    Trajs_Plevs=compute_VirtualTH(Trajs_Plevs)
    
    # plot variables
    gd = Geodesic()
    props_loc=dict(facecolor='white', alpha=0.85, edgecolor='k')
    props=dict(facecolor='white', alpha=0.85, edgecolor='none')
    #new colormap
    colMap_T=truncate_colormap('RdYlBu_r',minval=0.01,maxval=0.99,n=256)
    projex=ccrs.NorthPolarStereo()
    
    ########################### PLOT trajectories ####################
   
    # get starting positions
    Xo= [Trajs_Plevs[50][i].iloc[0]['lon'] for i in getList(Trajs_Plevs[50])] #lon
    Yo= [Trajs_Plevs[50][i].iloc[0]['lat'] for i in getList(Trajs_Plevs[50])] #lat
    
    ####### PLOT TRACKS, mark every 24h, tracks colored by black and STATISTIC about actual initial PRESSURE & start locations ###########
    
    # get the initial pressure at t=0 for all 6 levels
    IniP_df=pd.DataFrame(columns=['p','pAGL'])
    for p in hPa_AGL:
        Ps=[Trajs_Plevs[p][i]['p'].values[0] for i in getList(Trajs_Plevs[p])]
        p_agl=[p for i in range(len(Ps))]
        df_oneP=pd.DataFrame(data={'p' : Ps , 'pAGL' : p_agl})
        IniP_df=IniP_df.append(df_oneP)
        
    ######## PLOT

    fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(25,15),
                             subplot_kw={'projection':ccrs.NorthPolarStereo()})
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.85,hspace=0.1)
    ax=ax.ravel()
    for ip,p in enumerate(hPa_AGL):
        # set figure axis
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax[ip].set_boundary(circle, transform=ax[ip].transAxes)
        ax[ip].set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND,zorder=4,alpha=0.5)
        ax[ip].coastlines(zorder=6,color='#a57e52')
        ax[ip].annotate(panels[ip]+')',(.02,.9),xycoords='axes fraction',
                        fontsize=22,fontweight='bold',color='k')
        gl = ax[ip].gridlines(linewidth=1, color='k', alpha=0.3,zorder=5)
        gl.ylocator = mticker.FixedLocator([-10,20,50,60,70,80,90])
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180+60,60))
        #gl.xlines = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.yformatter = LONGITUDE_FORMATTER
        gl.n_steps = 90
      
        #title
        #ax[ip].set_title('Trajectories initialized\nat %s hPa AGL, 5d'%str(p),fontsize=12)
        ax[ip].text(155, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)
        
        text_iniP='P$_{mean}$ = %s, P$_{median}$ = %s,\nP$_{max}$ = %s, P$_{min}$ = %s'%(np.round(IniP_df[IniP_df.pAGL==p]['p'].mean(),1),
                                                                                          np.round(IniP_df[IniP_df.pAGL==p]['p'].median(),1),
                                                                                          np.round(IniP_df[IniP_df.pAGL==p]['p'].max(),1),
                                                                                          np.round(IniP_df[IniP_df.pAGL==p]['p'].min(),1))
        ax[ip].text(-155, 55, text_iniP,transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=15,color='k',zorder=9)
        plot_C_dif_24=truncate_colormap('PuRd',n=len(Trajs_Plevs[50][0].iloc[::24, :].iloc[1:]),minval=0.3,returncolors=True)
        
        lc=plottracks(Trajs_Plevs,p,ax=ax[ip],ma='o',COL='k',plotC=True,plotP=False,lw=1.5,
                      markDays=True,markEnd=False,fillm=True,plot24h_difcol=plot_C_dif_24,hdif=24,ms=4)

        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=120,zorder=15, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='*')
        
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        ax[ip].scatter(Xo,Yo,c='k',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')

        lon_ref=0
        lat_ref=90
        r = 90 - 70 #radii from 90N
        def compute_radius(ortho, radius_degrees):
            phi1 = lat_ref + radius_degrees if lat_ref <= 0 else lat_ref - radius_degrees
            _, y1 = ortho.transform_point(lon_ref, phi1, ccrs.PlateCarree())
            return abs(y1)
        # Compute the required radius in projection native coordinates:
        r_ortho = compute_radius(ccrs.NorthPolarStereo(), r) #70 N
        projx1, projy1 = projex.transform_point(lon_ref, lat_ref, ccrs.Geodetic())
        ax[ip].add_patch(mpatches.Circle(xy=[projx1, projy1],radius=r_ortho,linewidth=1.5,linestyle='-',alpha=1, transform=ccrs.NorthPolarStereo(),
                                              zorder=8,facecolor=None,fill=False, edgecolor='k'))

        ax[ip].text(75, 70.5, '70'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        ax[ip].text(75, 80.5, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)
        
        ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

    #legend
    labels = [str(l)+' days' for l in np.arange(1,5,1).tolist()[0:len(Trajs_Plevs[50][0].iloc[::24, :].iloc[1:])]]
    markers = [plt.Line2D((0,1),(0,0), color=c, marker='o', linestyle='',markerfacecolor=c,ms=10,lw=3)
               for c in plot_C_dif_24]

    leg= ax[3].legend(markers,labels, ncol=1, loc='center left', bbox_to_anchor=(1.01,1.16),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
      
    fig.suptitle('Forecast trajectories at %s'%Date.strftime('%d %B, %Y'),fontsize=20)
    savename='PLOT_oneC_Pstart_traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast

    #fig.savefig(savefigs + savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=150)
    #fig.savefig(savefigs + savename + '.png',bbox_inches = 'tight',dpi=350)
    #plt.show()

    ############ PLOT TRAJECTORIES BY PRESSURE #################
    fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(25,15),
                             subplot_kw={'projection':ccrs.NorthPolarStereo()})
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.85,hspace=0.1)
    ax=ax.ravel()
    for ip,p in enumerate(hPa_AGL):
        # set figure axis
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax[ip].set_boundary(circle, transform=ax[ip].transAxes)
        ax[ip].set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND,zorder=4,alpha=0.5)
        ax[ip].coastlines(zorder=6,color='#a57e52')
        ax[ip].annotate(panels[ip]+')',(.02,.9),xycoords='axes fraction',
                        fontsize=22,fontweight='bold',color='k')
        gl = ax[ip].gridlines(linewidth=1, color='k', alpha=0.3,zorder=5)
        gl.ylocator = mticker.FixedLocator([-10,20,50,60,70,80,90])
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180+60,60))
        #gl.xlines = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.yformatter = LONGITUDE_FORMATTER
        gl.n_steps = 90
        
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=120,zorder=15, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='*')
        
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        ax[ip].scatter(Xo,Yo,c='k',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')

        lon_ref=0
        lat_ref=90
        r = 90 - 70 #radii from 90N
        def compute_radius(ortho, radius_degrees):
            phi1 = lat_ref + radius_degrees if lat_ref <= 0 else lat_ref - radius_degrees
            _, y1 = ortho.transform_point(lon_ref, phi1, ccrs.PlateCarree())
            return abs(y1)
        # Compute the required radius in projection native coordinates:
        r_ortho = compute_radius(ccrs.NorthPolarStereo(), r) #70 N
        projx1, projy1 = projex.transform_point(lon_ref, lat_ref, ccrs.Geodetic())
        ax[ip].add_patch(mpatches.Circle(xy=[projx1, projy1],radius=r_ortho,linewidth=1.5,linestyle='-',alpha=1, transform=ccrs.NorthPolarStereo(),
                                              zorder=8,facecolor=None,fill=False, edgecolor='k'))

        ax[ip].text(75, 70.5, '70'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        ax[ip].text(75, 80.5, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)
        ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        #title
        ax[ip].text(155, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

        lc=plottracks(Trajs_Plevs,p,var='p',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,
                      markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)


    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.linspace(500,1000,6)) # pressure ticks
    cbar.ax.set_title('P (hPa)', fontweight='bold',fontsize=20) #pressure label
    cbar.ax.invert_yaxis() #invert pressure axis
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories at %s'%Date.strftime('%d %B, %Y'),fontsize=20)
    savename='PLOT_pressureC_traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=150)
    #plt.show()


    ########### PLOT TRAJECTORIES BY SPECIFIC HUMIDITY ##############

    fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(25,15),
                             subplot_kw={'projection':ccrs.NorthPolarStereo()})
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.85,hspace=0.1)
    ax=ax.ravel()
    for ip,p in enumerate(hPa_AGL):
        # set figure axis
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax[ip].set_boundary(circle, transform=ax[ip].transAxes)
        ax[ip].set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND,zorder=4,alpha=0.5)
        ax[ip].coastlines(zorder=6,color='#a57e52')
        ax[ip].annotate(panels[ip]+')',(.02,.9),xycoords='axes fraction',
                        fontsize=22,fontweight='bold',color='k')
        gl = ax[ip].gridlines(linewidth=1, color='k', alpha=0.3,zorder=5)
        gl.ylocator = mticker.FixedLocator([-10,20,50,60,70,80,90])
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180+60,60))
        #gl.xlines = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.yformatter = LONGITUDE_FORMATTER
        gl.n_steps = 90
        
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=120,zorder=15, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='*')
        
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        ax[ip].scatter(Xo,Yo,c='k',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')

        lon_ref=0
        lat_ref=90
        r = 90 - 70 #radii from 90N
        def compute_radius(ortho, radius_degrees):
            phi1 = lat_ref + radius_degrees if lat_ref <= 0 else lat_ref - radius_degrees
            _, y1 = ortho.transform_point(lon_ref, phi1, ccrs.PlateCarree())
            return abs(y1)
        # Compute the required radius in projection native coordinates:
        r_ortho = compute_radius(ccrs.NorthPolarStereo(), r) #70 N
        projx1, projy1 = projex.transform_point(lon_ref, lat_ref, ccrs.Geodetic())
        ax[ip].add_patch(mpatches.Circle(xy=[projx1, projy1],radius=r_ortho,linewidth=1.5,linestyle='-',alpha=1, transform=ccrs.NorthPolarStereo(),
                                              zorder=8,facecolor=None,fill=False, edgecolor='k'))

        ax[ip].text(75, 70.5, '70'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        ax[ip].text(75, 80.5, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)
        ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        #title
        ax[ip].text(155, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

                      
        lc=plottracks(Trajs_Plevs,p,var='Q',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(0,5.1,0.1),
                      markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)

    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(0,6,1)) # for specific humidity Q
    cbar.ax.set_title('Q (gkg$^{-1}$)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories at %s'%Date.strftime('%d %B, %Y'),fontsize=20)
    savename='PLOT_spechumC_traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=150)
    #plt.show()

    ######### PLOT TRAJECTORIES BY TEMPERATURE ###################
    
    fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(25,15),
                             subplot_kw={'projection':ccrs.NorthPolarStereo()})
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.85,hspace=0.1)
    ax=ax.ravel()
    for ip,p in enumerate(hPa_AGL):
        # set figure axis
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax[ip].set_boundary(circle, transform=ax[ip].transAxes)
        ax[ip].set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND,zorder=4,alpha=0.5)
        ax[ip].coastlines(zorder=6,color='#a57e52')
        ax[ip].annotate(panels[ip]+')',(.02,.9),xycoords='axes fraction',
                        fontsize=22,fontweight='bold',color='k')
        gl = ax[ip].gridlines(linewidth=1, color='k', alpha=0.3,zorder=5)
        gl.ylocator = mticker.FixedLocator([-10,20,50,60,70,80,90])
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180+60,60))
        #gl.xlines = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.yformatter = LONGITUDE_FORMATTER
        gl.n_steps = 90
        
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=120,zorder=15, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='*')
        
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        ax[ip].scatter(Xo,Yo,c='k',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')

        lon_ref=0
        lat_ref=90
        r = 90 - 70 #radii from 90N
        def compute_radius(ortho, radius_degrees):
            phi1 = lat_ref + radius_degrees if lat_ref <= 0 else lat_ref - radius_degrees
            _, y1 = ortho.transform_point(lon_ref, phi1, ccrs.PlateCarree())
            return abs(y1)
        # Compute the required radius in projection native coordinates:
        r_ortho = compute_radius(ccrs.NorthPolarStereo(), r) #70 N
        projx1, projy1 = projex.transform_point(lon_ref, lat_ref, ccrs.Geodetic())
        ax[ip].add_patch(mpatches.Circle(xy=[projx1, projy1],radius=r_ortho,linewidth=1.5,linestyle='-',alpha=1, transform=ccrs.NorthPolarStereo(),
                                              zorder=8,facecolor=None,fill=False, edgecolor='k'))

        ax[ip].text(75, 70.5, '70'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        ax[ip].text(75, 80.5, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)
        ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        #title
        ax[ip].text(155, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

             
        lc=plottracks(Trajs_Plevs,p,var='T',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(-20,20.2,0.2),
                      cmap=colMap_T,markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)

    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(-20,25,5)) # for Temperature in Celcius (T)
    cbar.ax.set_title('T ($^{\circ}$C)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories at %s'%Date.strftime('%d %B, %Y'),fontsize=20)
    savename='PLOT_tempC_traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=150)
    #plt.show()
    
    
    ####################### PLOT TRAJECTORIES BY THETA ###################

    fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(25,15),
                             subplot_kw={'projection':ccrs.NorthPolarStereo()})
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.85,hspace=0.1)
    ax=ax.ravel()
    for ip,p in enumerate(hPa_AGL):
        # set figure axis
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax[ip].set_boundary(circle, transform=ax[ip].transAxes)
        ax[ip].set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND,zorder=4,alpha=0.5)
        ax[ip].coastlines(zorder=6,color='#a57e52')
        ax[ip].annotate(panels[ip]+')',(.02,.9),xycoords='axes fraction',
                        fontsize=22,fontweight='bold',color='k')
        gl = ax[ip].gridlines(linewidth=1, color='k', alpha=0.3,zorder=5)
        gl.ylocator = mticker.FixedLocator([-10,20,50,60,70,80,90])
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180+60,60))
        #gl.xlines = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.yformatter = LONGITUDE_FORMATTER
        gl.n_steps = 90
        
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=120,zorder=15, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='*')
        
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        ax[ip].scatter(Xo,Yo,c='k',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')

        lon_ref=0
        lat_ref=90
        r = 90 - 70 #radii from 90N
        def compute_radius(ortho, radius_degrees):
            phi1 = lat_ref + radius_degrees if lat_ref <= 0 else lat_ref - radius_degrees
            _, y1 = ortho.transform_point(lon_ref, phi1, ccrs.PlateCarree())
            return abs(y1)
        # Compute the required radius in projection native coordinates:
        r_ortho = compute_radius(ccrs.NorthPolarStereo(), r) #70 N
        projx1, projy1 = projex.transform_point(lon_ref, lat_ref, ccrs.Geodetic())
        ax[ip].add_patch(mpatches.Circle(xy=[projx1, projy1],radius=r_ortho,linewidth=1.5,linestyle='-',alpha=1, transform=ccrs.NorthPolarStereo(),
                                              zorder=8,facecolor=None,fill=False, edgecolor='k'))

        ax[ip].text(75, 70.5, '70'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        ax[ip].text(75, 80.5, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)
        ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        #title
        ax[ip].text(155, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

        lc=plottracks(Trajs_Plevs,p,var='TH',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(260,300.5,0.5),
                    cmap=colMap_T,markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)
                      
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(260,305,5)) # for specific humidity Q
    cbar.ax.set_title(r'$\theta$ (K)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories at %s'%Date.strftime('%d %B, %Y'),fontsize=20)
    savename='PLOT_thetaC_traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast
    #fig.savefig(savefigs + savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=150)
    #plt.show()

    ####################### PLOT TRAJECTORIES BY THETA_v ###################

    fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(25,15),
                             subplot_kw={'projection':ccrs.NorthPolarStereo()})
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.85,hspace=0.1)
    ax=ax.ravel()
    for ip,p in enumerate(hPa_AGL):
        # set figure axis
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax[ip].set_boundary(circle, transform=ax[ip].transAxes)
        ax[ip].set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND,zorder=4,alpha=0.5)
        ax[ip].coastlines(zorder=6,color='#a57e52')
        ax[ip].annotate(panels[ip]+')',(.02,.9),xycoords='axes fraction',
                        fontsize=22,fontweight='bold',color='k')
        gl = ax[ip].gridlines(linewidth=1, color='k', alpha=0.3,zorder=5)
        gl.ylocator = mticker.FixedLocator([-10,20,50,60,70,80,90])
        gl.xlocator = mticker.FixedLocator(np.arange(-180,180+60,60))
        #gl.xlines = False
        gl.yformatter = LATITUDE_FORMATTER
        gl.yformatter = LONGITUDE_FORMATTER
        gl.n_steps = 90
        
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=120,zorder=15, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='*')
        
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        ax[ip].scatter(Xo,Yo,c='k',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')

        lon_ref=0
        lat_ref=90
        r = 90 - 70 #radii from 90N
        def compute_radius(ortho, radius_degrees):
            phi1 = lat_ref + radius_degrees if lat_ref <= 0 else lat_ref - radius_degrees
            _, y1 = ortho.transform_point(lon_ref, phi1, ccrs.PlateCarree())
            return abs(y1)
        # Compute the required radius in projection native coordinates:
        r_ortho = compute_radius(ccrs.NorthPolarStereo(), r) #70 N
        projx1, projy1 = projex.transform_point(lon_ref, lat_ref, ccrs.Geodetic())
        ax[ip].add_patch(mpatches.Circle(xy=[projx1, projy1],radius=r_ortho,linewidth=1.5,linestyle='-',alpha=1, transform=ccrs.NorthPolarStereo(),
                                              zorder=8,facecolor=None,fill=False, edgecolor='k'))

        ax[ip].text(75, 70.5, '70'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        ax[ip].text(75, 80.5, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)
        ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=10)

        #title
        ax[ip].text(155, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

        lc=plottracks(Trajs_Plevs,p,var='THv',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(260,300.5,0.5),
                    cmap=colMap_T,markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)
                      
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(260,305,5)) # for specific humidity Q
    cbar.ax.set_title(r'$\theta_v$ (K)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories at %s'%Date.strftime('%d %B, %Y'),fontsize=20)
    savename='PLOT_thetavC_traj_ARTofMELT_forecast_%s_70Nlatband_0h_4dfw'%Date_forecast
    #fig.savefig(savefigs + savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=150)
    #plt.show()
