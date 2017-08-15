#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Script to visualize World Ocean Database data around Canary Islands.

data acquired from: https://www.nodc.noaa.gov/cgi-bin/OC5/SELECT/builder.pl



Main aim: visualize historic trends in seawater  temperature  around Canary Islands.
Period: ~1950 - 2017

WOD:
    
    Boyer, T.P., J. I. Antonov, O. K. Baranova, C. Coleman, H. E. Garcia, 
    A. Grodsky, D. R. Johnson, R. A. Locarnini, A. V. Mishonov, T.D. O'Brien, 
    C.R. Paver, J.R. Reagan, D. Seidov, I. V. Smolyar, and M. M. Zweng, 
    2013: World Ocean Database 2013, NOAA Atlas NESDIS 72, S. Levitus, Ed., 
    A. Mishonov, Technical Ed.; Silver Spring, MD, 209 pp., 
    http://doi.org/10.7289/V5NZ85MT

Created on Sun Aug 13 16:14:12 2017

:author: Ruyman Azzollini
:contact: ruyman.azzollini_at_gmail.com

"""

# IMPORT STUFF
import os
import numpy as np
from glob import glob
from  wodpy import  wod
from pdb import set_trace as stop
import geoplotlib
from astropy import table  as tab
from vissim.support.files import cPickleRead, cPickleDumpDictionary
from datetime import date, datetime
from scipy import interpolate

import matplotlib
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['image.interpolation'] = 'none'
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
from mpl_toolkits.basemap import Basemap


# END IMPORT


best_exts = ['XBT','PFL','CTD','MBT']


fkY = 2000 # dummy leap year to allow input X-02-29 (leap day)
seasons_dict = [('winter', (date(fkY,  1,  1),  date(fkY,  3, 20))),
           ('spring', (date(fkY,  3, 21),  date(fkY,  6, 20))),
           ('summer', (date(fkY,  6, 21),  date(fkY,  9, 22))),
           ('autumn', (date(fkY,  9, 23),  date(fkY, 12, 20))),
           ('winter', (date(fkY, 12, 21),  date(fkY, 12, 31)))]

def get_season(now):
    """Dummy function to return "season" for a date"""
    if isinstance(now, datetime):
        now = now.date()
    now = now.replace(year=fkY)
    return next(season for season, (start, end) in seasons_dict
                if start <= now <= end)


def get_profiles(ffile,N=-1):
    """Extracts profiles from a WOD file."""
    
    fid = open(ffile)
    pfs = []
    counter=0
    while True:
        pf = wod.WodProfile(fid)
        if pf.is_last_profile_in_file(fid): break
        pfs.append(pf)
        counter += 1
        if N>=0:
            if counter == N: break
    return pfs
    

#def explore_wodfile(ffile):
#    pfs = get_profiles(ffile,N=1)
#    stop()
#    print pfs[0].columns
#    stop()
    
def parse_profiles(pfs):
    """Lat, Long, datetime, and temperature, pH, etc. at 
    shallowest depth (surface) for all profiles in pfs."""
    
    Npfs = len(pfs)
    
    lat = []
    lon = []
    datetime = []
    temp = []
    phosphate = []
    pH = []
    
    for ipf in range(Npfs):
        nlevels = pfs[ipf].n_levels()
        if nlevels>0:
            lat.append(pfs[ipf].latitude())
            lon.append(pfs[ipf].longitude())
            datetime.append(pfs[ipf].datetime())
            temp.append(pfs[ipf].t()[0])
            phosphate.append(pfs[ipf].phosphate()[0])
            pH.append(pfs[ipf].pH()[0])
    
    
    t = tab.Table()
    t['lat'] = tab.Column(lat)
    t['lon'] = tab.Column(lon)
    t['datetime'] = tab.Column(datetime)
    t['temp'] = tab.Column(temp)
    t['phosphate'] = tab.Column(phosphate)
    t['pH'] = tab.Column(pH)
    
    return t
    


def get_surf_temps(ffile):
    """Returns table with surface properties extracted from a WOD file."""
    pfs = get_profiles(ffile,N=-1)
    table = parse_profiles(pfs)
    return table

def get_all_data(exts,path):
    """Retrieves all data from WOD files of given extensions in path."""
    
    data = dict()    
    for ext in exts:
        print 'Loading extension %s...' % ext
        ffile = glob(os.path.join(path,'ocldb*.%s' % ext))[0]
        table = get_surf_temps(ffile)
        data[ext] = table

    return data


def big_loader(outfile,path):
    """Data Loader, main script."""
    
    exts = ['CTD','DRB', 'MBT', 'PFL', 'SUR', 'XBT']
    data = get_all_data(exts,path)
    cPickleDumpDictionary(data,outfile)
    return data


def merge_datasets(data,best_exts=best_exts):
    """Mertes tables from different data-sets"""
    ext0 = best_exts[0]
    ret = data[ext0]
    for ext in best_exts[1:]:
        ret = tab.vstack([ret,data[ext]])
    
    ret.sort(['datetime'])
    
    vseason = map(get_season,ret['datetime'])
    vyear = np.array([item.year for item in ret['datetime']])
    
    ret['season'] = tab.Column(vseason)
    ret['year'] = tab.Column(vyear)
    
    return ret
    
def plot_temp_vs_datenseason(mdata):
    """Plots temperature vs. year, for each season."""
    
    ysdata = mdata.group_by(['season','year'])
    
    seasons = ['winter','spring','summer','autumn']
    colors = ['b','g','r','m']
    
    pdata = dict()
    for season in seasons:
        pdata[season] = dict()
        pdata[season]['year'] = []
        pdata[season]['temp'] = []
    
    for key, group in zip(ysdata.groups.keys,ysdata.groups):
        season = key[0]
        iy = key[1]
        it = np.nanmean(group['temp'])
        pdata[season]['year'].append(iy)
        pdata[season]['temp'].append(it)
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    
    for ic,season in enumerate(seasons):
        y = np.array(pdata[season]['year'])
        t = np.array(pdata[season]['temp'])
        
        ixnonan = np.where((~np.isnan(y)) & (~np.isnan(t)))
        
        z = np.polyfit(y[ixnonan],t[ixnonan],2)
        p = np.poly1d(z)
        
        ax.plot(y,t,label=season,color=colors[ic])
        ax.plot(y,p(y),ls='--',color=colors[ic])
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels,loc='best')
    ax.set_xlabel('Year')
    ax.set_ylabel('Temp [C]')
    ax.set_xlim([1950,2020])
    ax.set_title('Surface Ocean Temperature - Canary Islands',fontsize=16)
    
    plt.show()


def get_tpol_by_sectors(mdata):
    """Retrieves a dictionary with temperatures vs. year, for each season, and segregated by 
    sectors of 1 x 1 degree in Canary Islands region."""
    
    ysdata = mdata.group_by(['season','year'])
    
    seasons = ['winter','spring','summer','autumn']
    
    nlon = 6
    nlat = 4
    
    vlon = np.linspace(-19,-13,nlon+1)
    vlat = np.linspace(26,30,nlat+1)
    
    poldata = dict()
    
    sect_coos = dict()
    
    for season in seasons:
        poldata[season] = dict()
        
        for ilon in range(nlon):
            for ilat in range(nlat):
                lon0 = vlon[ilon]
                lon1 = vlon[ilon+1]
                lat0 = vlat[ilat]
                lat1 = vlat[ilat+1]
                
                sector = 'S%i%i' % (ilon,ilat)
                
                sect_coos[sector] = (lon0+lon1)/2.,(lat0+lat1)/2.
                #sect_coos[sector] = lon0,lat0
                         
                poldata[season][sector] = dict()
                
                poldata[season][sector]['year'] = []
                poldata[season][sector]['temp'] = []
                poldata[season][sector]['ptemp'] = None
    
    poldata['coos'] = sect_coos
    
    for key, group in zip(ysdata.groups.keys,ysdata.groups):
        
        for ilon in range(nlon):
            for ilat in range(nlat):
                lon0 = vlon[ilon]
                lon1 = vlon[ilon+1]
                lat0 = vlat[ilat]
                lat1 = vlat[ilat+1]
                
                sector = 'S%i%i' % (ilon,ilat)
        
                ixsel = np.where((group['lat'] >= lat0) & (group['lat'] <= lat1) &\
                         (group['lon'] >= lon0) & (group['lon'] <= lon1))
        
                if len(ixsel[0]) > 0:
                    season = key[0]
                    iy = key[1]
                    it = np.mean(group['temp'][ixsel])
                    poldata[season][sector]['year'].append(iy)
                    poldata[season][sector]['temp'].append(it)
                else:
                    pass
                    #print 'empty: %i-%s, %s' % (iy,season,sector)
    
    
    sectors = poldata['coos'].keys()
    
    for season in seasons:
        for sector in sectors:
            
            y = np.array(poldata[season][sector]['year'])
            t = np.array(poldata[season][sector]['temp'])
            
            ixnonan = np.where((~np.isnan(y)) & (~np.isnan(t)))
        
            z = np.polyfit(y[ixnonan],t[ixnonan],1)
            p = np.poly1d(z)

            poldata[season][sector]['ptemp'] = p
    
    return poldata
    

def plot_temp_vs_t_bysector(mdata):
    """Plots temperatures vs. year, for each season, and segregated by 
    sectors of 1 x 1 degree in Canary Islands region."""
    
    seasons = ['winter','spring','summer','autumn']
    colors = ['b','g','r','m']
    
    poldata = get_tpol_by_sectors(mdata)
    nlat = 4
    nlon = 6
    
    #sectors = [key for key in poldata.keys() if key[0] == 'S']
    #sectors = np.sort(sectors)
    
    fig = plt.figure(figsize=(9,10))
    axs = []
    counter = 1
    
    for ilat in range(nlat-1,-1,-1):
        for ilon in range(nlon):
            
            sector = 'S%i%i' % (ilon,ilat)
    
            axs.append(fig.add_subplot(nlat,nlon,counter))
    
            for ic,season in enumerate(seasons):
                y = poldata[season][sector]['year']
                
                if len(y)>0:
                    p = poldata[season][sector]['ptemp']
                    #axs[-1].plot(y,t,label=season,color=colors[ic])
                    axs[-1].plot(y,p(y),ls='--', color=colors[ic], label=season)
                    
                
                if ic ==0: 
                    #axs[-1].set_title(sector,fontsize=10)
                    axs[-1].set_xlim([1950,2020])
                    axs[-1].set_ylim([17,24])
                    coo = poldata['coos'][sector]
                    cootxt = '%.1fW %.1fN' % (np.abs(coo[0]),coo[1])
                    axs[-1].text(0.1, 0.8,cootxt, horizontalalignment='left',
                       verticalalignment='center',
                       transform=axs[-1].transAxes,fontsize=6)
                if ilat != 0:
                    axs[-1].get_xaxis().set_ticks([])
                else:
                    axs[-1].set_xticks([1960,1990])
                    
                if ilon !=0:
                    axs[-1].get_yaxis().set_ticks([])
                else:
                    axs[-1].set_yticks([18,20,22,24])
                    
                axs[-1].tick_params(axis='both', which='major', labelsize=10)
                axs[-1].tick_params(axis='both', which='minor', labelsize=8)
            
            counter += 1
    
   
    
    handles,labels = axs[-1].get_legend_handles_labels()
    plt.figlegend(handles,labels,loc='right')

    
    plt.tight_layout()
    
    plt.subplots_adjust(top=0.90)
    plt.subplots_adjust(right=0.80)
    plt.subplots_adjust(wspace=0,hspace=0)
    
    
    plt.suptitle('Surface Ocean Temperature - Canary Islands',fontsize=12)
    

    plt.show()
    plt.close()



def plot_temp_vs_date(data):
    """Plot temperature vs. date. All data available."""
    
    fig = plt.figure(figsize=(8,7))
    ax = fig.add_subplot(111)
    
    
    for ext in data.keys():
        date = data[ext]['datetime'].data.copy()
        temp = data[ext]['temp'].data.copy()
        ax.plot(date,temp,label=ext)
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels,loc='best')
    ax.set_xlabel('Date',fontsize=15)
    ax.set_ylabel('Temp [C]',fontsize=15)
    ax.set_title('Av. Temp. Ocean - Canary Islands')
    plt.tight_layout()
    
    plt.show()
    

def show_data_on_map(data):
    """Shows probes locations on map of Can. Is."""
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.05,0.8,0.9])
    # create polar stereographic Basemap instance.
    m = Basemap(projection='stere',lon_0=-15.4,lat_0=28.1,lat_ts=28.1,\
            llcrnrlat=26.5,urcrnrlat=29.5,\
            llcrnrlon=-18.5,urcrnrlon=-13.5,\
            rsphere=6371200.,resolution='l',area_thresh=10000,epsg=4326)
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    m.drawlsmask(resolution='f')
    m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', 
                 service='ESRI_Imagery_World_2D', xpixels=400, ypixels=None, dpi=96, verbose=False,)
    # draw parallels.
    parallels = np.arange(0.,90,1.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(180.,360.,2.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    
    for ext in best_exts:
        sdata = data[ext]
        lat = sdata['lat'].data.copy()
        lon = sdata['lon'].data.copy()
        m.plot(lon,lat,'o',label = ext,alpha=0.5)
        #ax.plot(lat,lon,'o',label=ext)
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels,loc='lower right')
    ax.set_title('Probes Geographical  Distribution',fontsize=14)
    plt.show()


def show_avtemps_on_map(mdata,Y,season,diffY=0,suptitle='',figname=''):
    """Shows average temperatures for selected years, or difference
    in av. temperature between 2 years, on map of Can. Is."""
    
    seasons = ['winter','spring','summer','autumn']
    poldata = get_tpol_by_sectors(mdata)
    #nlat = 2
    #nlon = 4
    sect_coos_dict = poldata['coos']
    sectors = np.sort(sect_coos_dict.keys())
    
    vlat = []
    vlon = []
    vtemp = []
    
    
    dlon = 0.5
    dlat = 0.5
    
    for sector in sectors:
        ilat = sect_coos_dict[sector][1]
        vlat += [ilat-dlat,ilat+dlat,ilat+dlat,ilat-dlat]
        ilon = sect_coos_dict[sector][0]
        vlon += [ilon-dlon,ilon-dlon,ilon+dlon,ilon+dlon]
        itemp = poldata[season][sector]['ptemp'](Y)
        if diffY !=0:
            itemp -= poldata[season][sector]['ptemp'](diffY)
        
        vtemp += [itemp]*4
        
    #f = interpolate.interp2d(vlon,vlat,vtemp,kind='linear')
    #f = interpolate.bisplrep(vlon,vlat,vtemp,s=1)

    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # create polar stereographic Basemap instance.
    m = Basemap(projection='stere',lon_0=-15.4,lat_0=28.1,lat_ts=28.1,\
            llcrnrlat=26.5,urcrnrlat=29.5,\
            llcrnrlon=-18.5,urcrnrlon=-13.5,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
    #4326
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    m.drawlsmask(resolution='f')
    #m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', 
    #             service='ESRI_Imagery_World_2D', xpixels=400, ypixels=None, dpi=96, verbose=False,)
    # draw parallels.
    parallels = np.arange(0.,90,1.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(180.,360.,2.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    
    lons, lats = m.makegrid(100, 100) # get lat/lons of ny by nx evenly space grid.
    
    points = np.stack((vlon,vlat),axis=1)
    
    data = interpolate.griddata(points,vtemp,(lons,lats),method='linear')
    
    #data = f(lons[0,:],lats[:,0]).reshape(lons.shape)
    x, y = m(lons, lats) # compute map proj coordinates.
    # draw filled contours.
    if diffY ==0:
        clevs = np.linspace(16,27,30)
    else:
        clevs = np.linspace(-2,+4,10)
        
    cs = m.contourf(x,y,data,clevs,cmap=cm.rainbow,alpha=0.5)
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad='5%')
    
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    
    cbar.set_label('T [C]')
    
    ax.set_title(suptitle)
    
    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()




if __name__ == '__main__':
    
    doLoad = False
    datapck = 'can_water_temp_wider.pick'
    path = 'data_wider'
    if doLoad:
        data = big_loader(datapck,path)
    else:
        data = cPickleRead(datapck)
    
    show_data_on_map(data)
    plot_temp_vs_date(data)
    

    
    mdata = merge_datasets(data,best_exts)
    
    plot_temp_vs_datenseason(mdata)
    
    plot_temp_vs_t_bysector(mdata)
    
    seasons = ['winter','spring','summer','autumn']
    for season in seasons:
        for Y in [1955,2015]:
            show_avtemps_on_map(mdata,Y,season,diffY=0,\
                                suptitle='%s - %i' % (season,Y),
                                                     figname='avTemp_%s_%s.png' %\
                                                     (season,Y))
        show_avtemps_on_map(mdata,2015,season,diffY=1955,
                            suptitle=r'$%s\ -\ \Delta T[2015-1955]$' % season,
                                                       figname='avDeltaTemp_%s_%i_%i.png' % \
                                                       (season,2015,1955))
    
    
    
    