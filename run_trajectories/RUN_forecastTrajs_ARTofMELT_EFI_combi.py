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
from plot_functions_forecastTrajs_ARTofMELT import *

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
    Variabels along the trajectory will be traced for T q and TH
    '''
    subprocess.call(shlex.split('./%s %s'%(filename,forecastdate)))

def linkPfiles():
    '''
    This function links the forecast Pfiles from a directory defined in the bash script to the current running directory
    '''
    subprocess.call('./linkPfiles.sh')


def remlinkPfiles():
    '''
    This function removes the link to the forecast Pfiles from the current rundirectory  
    '''
    subprocess.call('./removelinkedPfiles.sh')


#python scripts in the .py script "plot_functions_forecastTrajs_ARTofMELT.py"


if __name__ == "__main__":
   
    ########################Arguments for the run##################

    Y=int(sys.argv[1]) #same year 2023
    M=int(sys.argv[2]) #month of forecast initialisation
    d=int(sys.argv[3]) # day of forecast initialisation
    h=00 #forecast initialized at 00 UTC default time

    ####################### Define dates ###################
    Date_forecast=str(Y) + str(M).zfill(2) + str(d).zfill(2) + '_' + str(h).zfill(2)
    print(Date_forecast)
    Date=datetime.datetime.strptime(Date_forecast, '%Y%m%d_%H') #forecast date in different format for plots
    
    ############################# Calculate trajectories with LAGRANTO ###################
    
    #### if trajectories are already computed, comment these following functions! ###
    
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
    
    # define plot variables
    gd = Geodesic()
    props_loc=dict(facecolor='white', alpha=0.85, edgecolor='k')
    props=dict(facecolor='white', alpha=0.85, edgecolor='none')
    #new colormap
    colMap_T=truncate_colormap('RdYlBu_r',minval=0.01,maxval=0.99,n=256)
    
    #get the file for location info (user-defined)
    #open given location-textfile in pandas
    loc_EFI=pd.read_csv('EFIloc.txt', sep=' ',header=None) #this file is an appended file (one row per run)
    loc_EFI.columns=['Date','hours','lon','lat','lon_int','lat_int'] #lat/lon "_int" refers to integral number of given decimal locations (used for the trajectory file names)
    LOC_sel=loc_EFI.iloc[-1] #get the most recent run (works better for the combined script)
    if LOC_sel['Date'] != Date_forecast:
        raise ValueError('No matching Date for %s'%Date_forecast)
        
    #extract information needed for the last run
    hours=LOC_sel.hours #hours into the forecast run to start trajs.
    lon_sel=LOC_sel.lon #initial location: longitude
    lat_sel=LOC_sel.lat #initial location: latitude
    lon_int=LOC_sel.lon_int #initial location in integer of longitude
    lat_int=LOC_sel.lat_int #initial location in integer of latitude
    
    ######################### READ in trajectories ###########################
    
    #define the folder accordingly (use this if the traj-file is in same folder as rundirectory)
    TrajFile='traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))
    
    #read the selected textfile with Traj info
    df_traj=pd.read_csv(TrajFile,encoding='latin-1', delimiter=r"\s+" ,skiprows=[0,1,3],header=[0])
    Dur_h=pd.read_csv(TrajFile,encoding='latin-1', delimiter=r"\s+",nrows=1,header=None).iloc[:,6][0]/60. #for hours
    # get the steps per trajectory & amount of trajectories
    w=sorted(list(dict.fromkeys(df_traj['time'].diff().values))) #difference in time of one timestep [hourly]
    cleanedListH = [x for x in w if str(x) not in ['nan', str(Dur_h), str(abs(Dur_h)), str(Dur_h*-1)]][0] # delete nan and difference between trajectories
    nSteps_h = abs(Dur_h /cleanedListH) #how many steps for one trajectory
    nTrajs=len(df_traj)/(nSteps_h+1) #amount of trajectories in the whole file (including all from different levels)
    
    #dictionary for each trajectory, keys = initial pressure levels
    Trajs_Plevs=getonetraj_multipleP(nSteps_h,nTrajs,df_traj,plevs=hPa_AGL)

    #compute virtual potential temperature
    Trajs_Plevs=compute_VirtualTH(Trajs_Plevs)
    
    ########################### PLOT trajectories ####################
    # get starting positions
    Xo= [Trajs_Plevs[50][i].iloc[0]['lon'] for i in getList(Trajs_Plevs[50])] #lon
    Yo= [Trajs_Plevs[50][i].iloc[0]['lat'] for i in getList(Trajs_Plevs[50])] #lat
    
    
    fig, ax = plt.subplots(ncols=1, nrows=1,figsize=(10,10),facecolor='white',
                         constrained_layout=True,subplot_kw={'projection':ccrs.NorthPolarStereo()})
    ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
    ax.coastlines(zorder=2,linewidth=0.9,color='k',alpha=1.)
    gl = ax.gridlines(linewidth=1, color='k', alpha=0.3,zorder=2.)
    gl.ylocator = mticker.FixedLocator(np.arange(0,95,5))
    gl.xlocator = mticker.FixedLocator(np.arange(-180,180+30,30))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.n_steps = 90

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    ax.add_feature(cfeature.LAND)
    ax.scatter(Xo,Yo,c='r',s=10,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='o')
    ax.scatter(lon_sel,lat_sel,c='k',s=50,zorder=10, edgecolor='g',transform=ccrs.PlateCarree(),marker='o')

    ax.scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
    #ax.text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

    ax.text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)
    ax.text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)
    
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax.legend(markers_2,labels_2, ncol=1, loc='upper right',
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 15)
                   
    cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
    ax.add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='k', alpha=0.8,lw=2,facecolor='none')

    savename='PLOT_traj_ARTofMELT_forecast_%s_EFI_startingpoints_%s_%s_%sh'%(Date_forecast,str(lon_int),str(lat_int),str(hours))
    #fig.savefig(savefigpath + savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)


    ####### PLOT TRACKS, mark every 24h, tracks colored by black and STATISTIC about actual initial PRESSURE ###########
    
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
        ax[ip].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
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
        
        cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
        ax[ip].add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='r', alpha=0.8,lw=2,facecolor='none')
        ax[ip].scatter(lon_sel,lat_sel,c='k',s=10,zorder=12, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='o')
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        #title
        ax[ip].text(140, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)
        
        text_iniP='P$_{mean}$ = %s, P$_{median}$ = %s,\nP$_{max}$ = %s, P$_{min}$ = %s'%(np.round(IniP_df[IniP_df.pAGL==p]['p'].mean(),1),
                                                                                          np.round(IniP_df[IniP_df.pAGL==p]['p'].median(),1),
                                                                                          np.round(IniP_df[IniP_df.pAGL==p]['p'].max(),1),
                                                                                          np.round(IniP_df[IniP_df.pAGL==p]['p'].min(),1))
        ax[ip].text(-155, 45, text_iniP,transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=15,color='k',zorder=9)
        plot_C_dif_24=truncate_colormap('PuRd',n=len(Trajs_Plevs[50][0].iloc[::24, :].iloc[1:]),minval=0.3,returncolors=True)
        
        lc=plottracks(Trajs_Plevs,p,ax=ax[ip],ma='o',COL='k',plotC=True,plotP=False,lw=1.5,
                      markDays=True,markEnd=False,fillm=True,plot24h_difcol=plot_C_dif_24,hdif=24,ms=4)

        if ip == 0:
            ax[ip].text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

            ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)


    #legend
    labels = [str(l)+' days' for l in np.arange(1,5,1).tolist()[0:len(Trajs_Plevs[50][0].iloc[::24, :].iloc[1:])]]
    markers = [plt.Line2D((0,1),(0,0), color=c, marker='o', linestyle='',markerfacecolor=c,ms=10,lw=3)
               for c in plot_C_dif_24]

    leg= ax[3].legend(markers,labels, ncol=1, loc='center left', bbox_to_anchor=(1.03,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
                   
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax[4].legend(markers_2,labels_2, ncol=1, loc='center left', bbox_to_anchor=(.92,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)

    fig.suptitle('Forecast trajectories +%sh from %s'%(str(hours),Date.strftime('%d %B, %Y')),fontsize=20)
    savename='PLOT_oneC_Pstart_traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)
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
        ax[ip].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
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
    
        cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
        ax[ip].add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='k', alpha=0.8,lw=2,facecolor='none')
        ax[ip].scatter(lon_sel,lat_sel,c='k',s=10,zorder=12, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='o')
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        #title
        ax[ip].text(140, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)
   
        lc=plottracks(Trajs_Plevs,p,var='p',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,
                      markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)

        if ip == 0:
            ax[ip].text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

            ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.linspace(500,1000,6)) # pressure ticks
    cbar.ax.set_title('P (hPa)', fontweight='bold',fontsize=20) #pressure label
    cbar.ax.invert_yaxis() #invert pressure axis
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)
    
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax[4].legend(markers_2,labels_2, ncol=1, loc='center left', bbox_to_anchor=(.92,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
                   
    fig.suptitle('Forecast trajectories +%sh from %s'%(str(hours),Date.strftime('%d %B, %Y')),fontsize=20)
    savename='PLOT_pressureC_traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)
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
        ax[ip].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
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
        
        cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
        ax[ip].add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='k', alpha=0.8,lw=2,facecolor='none')
        ax[ip].scatter(lon_sel,lat_sel,c='k',s=10,zorder=12, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='o')
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        #title
        ax[ip].text(140, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

        lc=plottracks(Trajs_Plevs,p,var='Q',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(0,5.1,0.1),
                      markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)
        
        if ip == 0:
            ax[ip].text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

            ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)
            
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax[4].legend(markers_2,labels_2, ncol=1, loc='center left', bbox_to_anchor=(.92,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
                   
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(0,6,1)) # for specific humidity Q
    cbar.ax.set_title('Q (gkg$^{-1}$)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories +%sh from %s'%(str(hours),Date.strftime('%d %B, %Y')),fontsize=20)
    savename='PLOT_spechumC_traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)
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
        ax[ip].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
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
        
        cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
        ax[ip].add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='k', alpha=0.8,lw=2,facecolor='none')
        ax[ip].scatter(lon_sel,lat_sel,c='k',s=10,zorder=12, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='o')
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        #title
        ax[ip].text(140, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)
    

        lc=plottracks(Trajs_Plevs,p,var='T',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(-20,20.2,0.2),
                      cmap=colMap_T,markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)
        
        if ip == 0:
            ax[ip].text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

            ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)
        
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax[4].legend(markers_2,labels_2, ncol=1, loc='center left', bbox_to_anchor=(.92,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
                   
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(-20,25,5)) # for Temperature in Celcius (T)
    cbar.ax.set_title('T ($^{\circ}$C)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)

    fig.suptitle('Forecast trajectories +%sh from %s'%(str(hours),Date.strftime('%d %B, %Y')),fontsize=20)
    savename='PLOT_tempC_traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))
    #fig.savefig(savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)
    #plt.show()
    
    
    ####################### PLOT TRAJECTORIES BY THETA ###################
    
    # make a plot of trajectories
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
        ax[ip].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
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
        
        cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
        ax[ip].add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='k', alpha=0.8,lw=2,facecolor='none')
        ax[ip].scatter(lon_sel,lat_sel,c='k',s=10,zorder=12, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='o')
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        #title
        #ax[ip].set_title('Trajectories initialized\nat %s hPa AGL, 5d'%str(p),fontsize=12)
        ax[ip].text(140, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

        lc=plottracks(Trajs_Plevs,p,var='TH',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(260,300.5,0.5),
                            cmap=colMap_T,markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)

        if ip == 0:
            ax[ip].text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

            ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)
        
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax[4].legend(markers_2,labels_2, ncol=1, loc='center left', bbox_to_anchor=(.92,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
                   
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(260,305,5)) # for specific humidity Q
    cbar.ax.set_title(r'$\theta$ (K)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)
      
    fig.suptitle('Forecast trajectories +%sh from %s'%(str(hours),Date.strftime('%d %B, %Y')),fontsize=20)
    savename='PLOT_thetaC_traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))

    #fig.savefig(savefigs + savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)

    #plt.show()
    
    ####################### PLOT TRAJECTORIES BY THETA_v ###################
    
    # make a plot of trajectories
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
        ax[ip].set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
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
        
        cp = gd.circle(lon=lon_sel, lat=lat_sel, radius=180000.)
        ax[ip].add_geometries([sgeom.Polygon(cp)], crs=ccrs.PlateCarree(),zorder=10, edgecolor='k', alpha=0.8,lw=2,facecolor='none')
        ax[ip].scatter(lon_sel,lat_sel,c='k',s=10,zorder=12, edgecolor='g',transform=ccrs.PlateCarree(),
                       marker='o')
        ax[ip].scatter(loc['lon'],loc['lat'],c='k',s=100,zorder=10, edgecolor='k',transform=ccrs.PlateCarree(),marker='*')
        #ax[ip].text(loc['lon'],loc['lat']+2, 'Oden',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

        #title
        #ax[ip].set_title('Trajectories initialized\nat %s hPa AGL, 5d'%str(p),fontsize=12)
        ax[ip].text(140, 60, str(p)+'hPa AGL',transform=ccrs.Geodetic(),bbox=props_loc,
                    fontsize=25,color='k',zorder=9)

        lc=plottracks(Trajs_Plevs,p,var='THv',ax=ax[ip],ma='o',plotC=False,plotP=True,lw=1.5,colorbounds=np.arange(260,300.5,0.5),
                            cmap=colMap_T,markDays=False,markEnd=True,fillm=False,hdif=24,ms=4)

        if ip == 0:
            ax[ip].text(75, 81, '80'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)

            ax[ip].text(75, 61, '60'+r'$^{\circ}$N',transform=ccrs.Geodetic(),bbox=props,fontsize=12,color='k',zorder=9)
        
    labels_2=['Oden','chosen start-\nlocation']
    markers_2=[plt.Line2D((0,1),(0,0), color='k', marker='*', linestyle='',markerfacecolor='k',ms=12,lw=3)]
    markers_2.append(plt.Line2D((0,1),(0,0), markerfacecolor='k', marker='o', linestyle='',
                                markeredgewidth=3,markeredgecolor='g',ms=7,lw=3))
    leg2=ax[4].legend(markers_2,labels_2, ncol=1, loc='center left', bbox_to_anchor=(.92,1.1),
                   columnspacing=1.0,  labelspacing=0.23, handletextpad=0.1, handlelength=1,fancybox=True,
                   shadow=True,fontsize = 20)
                   
    # make colorbar if return LC - colored trajectories with the third variable
    cbar_ax = fig.add_axes([0.89, 0.1, 0.015, 0.9-0.1])
    cbar=fig.colorbar(lc, cax=cbar_ax, extend='both')
    cbar.set_ticks(np.arange(260,305,5)) # for specific humidity Q
    cbar.ax.set_title(r'$\theta_v$ (K)', fontweight='bold',fontsize=20) #for specific humidity Q
    cbar.solids.set(alpha=1)
    plt.setp(cbar_ax.get_yticklabels(), fontsize=25)
      
    fig.suptitle('Forecast trajectories +%sh from %s'%(str(hours),Date.strftime('%d %B, %Y')),fontsize=20)
    savename='PLOT_thetavC_traj_ARTofMELT_forecast_%s_EFI_%s_%s_%sh_4dfw'%(Date_forecast,str(lon_int),str(lat_int),str(hours))

    #fig.savefig(savefigs + savename + '.pdf',bbox_inches = 'tight',format='pdf',dpi=150)
    fig.savefig(savename + '.png',bbox_inches = 'tight',dpi=100)

    #plt.show()


    print('Plotting done for forward trajectories done on %s forecast initialized at +%sh from (%s lon, %s lat)'%(Date_forecast,str(hours),str(lon_sel),str(lat_sel)))
