import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap, cm, shiftgrid, addcyclic
from dateutils import splitdate,dayofyear,dateshift
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import collections
import pygrib
import sys, os, glob
from scipy.ndimage.filters import minimum_filter, maximum_filter

def extrema(mat,mode='wrap',window=10):
    """find the indices of local extrema (min and max)
    in the input array."""
    mn = minimum_filter(mat, size=window, mode=mode)
    mx = maximum_filter(mat, size=window, mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    # (mat == mn) true if pixel is equal to the local in
    # Return the indices of the maxima, minima
    return np.nonzero(mat == mn), np.nonzero(mat == mx)

def plotlows(mslpmean):
    # the window parameter controls the number of highs and lows detected.
    # (higher value, fewer highs and lows)
    local_min, local_max = extrema(mslpmean, mode='wrap', window=40)
    xlows = x[local_min]; xhighs = x[local_max]
    ylows = y[local_min]; yhighs = y[local_max]
    lowvals = mslpmean[local_min]; highvals = mslpmean[local_max]
    # plot lows as blue L's, with min pressure value underneath.
    xyplotted = []
    yoffset = 0.05*(map.ymax-map.ymin)
    dmin = yoffset
    for xx,yy,p in zip(xlows, ylows, lowvals):
        if xx < map.xmax and xx > map.xmin and yy < map.ymax and yy > map.ymin:
            dist = [np.sqrt((xx-x0)**2+(yy-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                plt.text(xx,yy,'L',fontsize=14,fontweight='bold',
                        ha='center',va='center',color='b')
                plt.text(xx,yy-yoffset,repr(int(p)),fontsize=7,
                        ha='center',va='top',color='b',
                   bbox=dict(boxstyle="square",ec='None',fc=(0.9,0.9,0.9,0.5)))
                xyplotted.append((xx,yy))

def makehash():
    return collections.defaultdict(makehash)

date = sys.argv[1]
fhr = sys.argv[2]
gribPath = sys.argv[3]
fimensPath = sys.argv[4]
tcvitalsPath = sys.argv[5]
nmembers = sys.argv[6]
nmembers=int(nmembers)
nanals = nmembers 
print 'in plotens_FIM.py: date: ' , date
print 'in plotens_FIM.py: gribPath: ' , gribPath
print 'in plotens_FIM.py: fimensPath: ' , fimensPath
print 'in plotens_FIM.py: tcvitalsPath: ' , tcvitalsPath
print 'in plotens_FIM.py: nmembers: ' , nmembers
print 'in plotens_FIM.py: nanals: ' , nanals
ifhr = int(fhr); ifhrinc = 6
print 'in plotens_FIM.py: fhr: ' , fhr
print 'in plotens_FIM.py: ifhr: ' , ifhr
yyyy,mm,dd,hh = splitdate(date)
julday = dayofyear(yyyy,mm,dd)+1
hr = date[8:10]

domains = []
tclats  = makehash()
tclons  = makehash()
directories  = makehash()

tcvitals = os.path.join(tcvitalsPath, "tcvitals.%s.txt" % date)
print 'tcvitals:  ', tcvitals

try:
    for line in open(tcvitals):
        print 'line: ',line
        storm = line.split()[1]
        stormName = line.split()[2]
        storm = storm.lstrip()
        storm = storm.rstrip()
        stormName = stormName.lstrip()
        tclat = line.split()[5]
        tclon = line.split()[6]
        if tclat[-1] == 'N':
           tclat = 0.1*float(tclat[:-1])
        else:
            tclat = -0.1*float(tclat[:-1])
        if tclon[-1] == 'E':
           tclon = 0.1*float(tclon[:-1])
        else:
           tclon = -0.1*float(tclon[:-1])
        directories[storm]='%s_%s' % (stormName,storm)
        print 'storm: ',storm,' directories[',storm,']: ', directories[storm]
        tclats[storm]=tclat
        tclons[storm]=tclon
        domains.append(storm)
        ensDir = os.path.join(fimensPath,directories[storm])
        print 'ensDir:  ', ensDir
        try:
           os.stat(ensDir)
        except:
           print 'directory ', ensDir, ' does not exist'
           try:
             print '....making directory', ensDir 
             os.makedirs(ensDir)
           except:
             print 'ERROR - COULD NOT MAKE: ', ensDir
except:
    print 'IN EXCEPT'
    pass

for domain in domains:
    print 'DOMAIN: ',domain

domains.insert(0,'West_Pacific')
domains.insert(0,'East_Pacific')
domains.insert(0,'Atlantic')
directories['West_Pacific']='West_Pacific'
directories['East_Pacific']='East_Pacific'
directories['Atlantic']='Atlantic'
tclats['East_Pacific']=0.0
tclats['West_Pacific']=0.0
tclats['Atlantic']=0.0
tclons['West_Pacific']=0.0
tclons['East_Pacific']=0.0
tclons['Atlantic']=0.0

# read forecast tracks.
# set default to FIM tracker file
datapath = gribPath
print 'datapath: ',datapath
trackerfile = os.path.join(datapath,"tracker_01/%i/track.%s00.FE01" % (ifhr,date))
print 'trackerfile',trackerfile

fcstracks=makehash()
print 'trackerfile: ',trackerfile
try:
    for line in open(trackerfile,'r'):
        linesplit = line.split(',')
      # WP, 04, 2011052500, 03, GE00, 000, 124N, 1281E,
        print 'line: ',line
        basin = linesplit[0]
        storm = linesplit[1]
        print 'basin: ',basin,' storm: ',storm
        if basin == 'WP':
            storm = storm+'W'
        elif basin == 'NA':
            storm = storm+'A'
        elif basin == 'EP':
            storm = storm+'E'
        elif basin == 'AL':
            storm = storm+'L'
        else:
            continue
        storm = storm.lstrip()
        fhour = int(linesplit[5])
        lat = linesplit[6]
        if lat[-1] == 'N':
            lat = 0.1*float(lat[:-1])
        else:
            lat = -0.1*float(lat[:-1])
        lon = linesplit[7]
        if lon[-1] == 'E':
            lon = 0.1*float(lon[:-1])
        else:
            lon = -0.1*float(lon[:-1])
        minp = int(linesplit[9])
        fcstracks[storm][fhour] = (lat,lon,minp)
        print 'storm: **',storm
        print 'fcstracks[',storm,']: ',fcstracks[storm][fhour]

except:
    pass

lat_0, lon_0 = None, None

for domain in domains:
    print '******************************'
    print 'DOMAIN',domain
    strmid = domain
    if domain == 'Atlantic':
        lat_0 = 22.; lon_0=-55.
        width=9000.e3; height=4500.e3
    elif domain == 'West_Pacific':
        lat_0 = 22.; lon_0=135.
        width=9000.e3; height=4500.e3
    elif domain == 'East_Pacific':
        lat_0 = 22.; lon_0=360.-125
        width=9000.e3; height=4500.e3
    else:
        width=6000.e3; height=3000.e3
        ifhrx = ifhr; foundit = False
        lat_0 = None; lon_0 = None
        while not foundit:
            try:
                lat_0,lon_0,minp = fcstracks[strmid][ifhrx]
                foundit = True
                print domain, 'found !'
            except:
                lat_0=tclats[strmid]
                lon_0=tclons[strmid]
                foundit = False
                print domain, 'not found - skipping.....!'
                ifhrx = ifhrx - ifhrinc
                if ifhrx <  6: break

    if lat_0 is not None:
        print strmid, lat_0, lon_0
        
        
        map = Basemap(projection='lcc',lon_0=lon_0,lat_0=lat_0,width=width,height=height)
        # map2 =\
        # Basemap(projection='stere',lon_0=lon_0,lat_0=lat_0,width=1.5*width,height=1.5*height,resolution='l')
        map2 = map
        xx2, yy2 = map2(lon_0, lat_0)
        xx, yy = map(lon_0, lat_0)
        lats = None; CS=None
        for nanal in range(1,nanals+1):
            print '*** NANAL: ',nanal

            datapath = gribPath
            filename = os.path.join(datapath,"post_%02i/fim/NAT/grib1/%s%03i%s000%03i"%(nanal,date[2:4],julday,hr,int(fhr)))
            print 'datapath: ',datapath
            print 'filename: ',filename

            try:
                grbs = pygrib.open(filename)
            except:
                print "can't open ", filename
                continue
            try:
                grbu=grbs.select(shortName='10u')[0]
                print '10u read ..'
                grbv=grbs.select(shortName='10v')[0]
                print '10v read ..'
                grbpwat=grbs.select(shortName='pwat')[0]
                print 'pwat read ..'
                grbmslp=grbs.select(typeOfLevel='meanSea')[0]
                print 'msl read ..'
            except:
                print 'grib variable is missing ..'
            umean = grbu.values
            vmean = grbv.values
            # swap data
            umean = umean[::-1,:]        
            vmean = vmean[::-1,:]              
            # convert mps to knots.
            umean = 1.943844*umean
            vmean = 1.943844*vmean
            try:
              mslpmean = 0.01*grbmslp.values
              mslpmean = mslpmean[::-1,:]                #swap data
              isMSL = True
            except:
              print 'mean sea level pressure missing!'
              isMSL = False
            try:
              pwat = grbpwat.values
              pwat = pwat[::-1,:]                        #swap data
              isPWAT = True
            except:
              print 'precipitable water missing!'
              isPWAT = False
            if fhr > '000':
                try:
                    isPrecip = True
                    # use for older grib output where 6h acc used timerange = 0
                    #grbapcp = grbs.select(parameterName='TP Total precipitation kg m**-2',timeRangeIndicator=0)[0]

                    # assumes 6h accum is second APCP record in file;  [0] would get the total precip field (0-fhr acc)
                    grbapcp = grbs.select(parameterName='TP Total precipitation kg m**-2',timeRangeIndicator=4)[1]
                    precip = grbapcp.values
                    precip = precip[::-1,:]             #swap data
                    print 'APCP read ..'
                except:
                    isPrecip = False
                    print 'APCP missing ..'
            if nanal == 1:
                lats, lons = grbu.latlons()
                lats1d = -lats[:,0]
                lons1d = (360./lons.shape[1])*np.arange(lons.shape[1])
                lons, lats = np.meshgrid(lons1d,lats1d)
            # create 3 figs
            if nanal == 1:
                fig1 = plt.figure(num=1,figsize=(17,14))
                fig1.subplots_adjust(bottom=0.1)
                fig1.subplots_adjust(left=0.125)
                fig1.subplots_adjust(top=0.9)
                fig1.subplots_adjust(right=0.9)
                fig1.subplots_adjust(hspace=0.15)
                fig1.subplots_adjust(wspace=0.15)
                fig2 = plt.figure(num=2,figsize=(17,14))
                fig2.subplots_adjust(bottom=0.05)
                fig2.subplots_adjust(left=0.125)
                fig2.subplots_adjust(top=0.9)
                fig2.subplots_adjust(right=0.9)
                fig2.subplots_adjust(hspace=0.15)
                fig2.subplots_adjust(wspace=0.15)
                fig3 = plt.figure(num=3,figsize=(17,14))
                fig3.subplots_adjust(bottom=0.1)
                fig3.subplots_adjust(left=0.125)
                fig3.subplots_adjust(top=0.9)
                fig3.subplots_adjust(right=0.9)
                fig3.subplots_adjust(hspace=0.15)
                fig3.subplots_adjust(wspace=0.15)
                fig4 = plt.figure(num=4,figsize=(17,14))
                fig4.subplots_adjust(bottom=0.1)
                fig4.subplots_adjust(left=0.125)
                fig4.subplots_adjust(top=0.9)
                fig4.subplots_adjust(right=0.9)
                fig4.subplots_adjust(hspace=0.15)
                fig4.subplots_adjust(wspace=0.15)
            # 10m wind plot
            plt.figure(1)
            ax1=plt.subplot(4,3,nanal)           ## 4 rows x 3 cols
            nx = int((map.xmax-map.xmin)/50000.)+1; ny = int((map.ymax-map.ymin)/50000.)+1
            umean, lons1 = addcyclic(umean, lons1d)
            umean,lons2 = shiftgrid(180.,umean,lons1,start=False)
            vmean, lons1 = addcyclic(vmean, lons1d)
            vmean,lons2 = shiftgrid(180.,vmean,lons1,start=False)
            umean, vmean, x, y = map.transform_vector(umean,vmean,lons2,lats1d,nx,ny,returnxy=True)
            spd = np.sqrt(umean**2 + vmean**2)
            print spd.shape,spd.min(), spd.max()
            map.drawcoastlines(ax=ax1)
            map.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=8,ax=ax1)
            map.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=8,ax=ax1)
            CS0=map.contourf(x,y,spd,np.arange(8,81,8),cmap=cm.GMT_haxby_r,extend='max',ax=ax1)
            nskip = 5
            brb =\
            map.barbs(x[::nskip,::nskip],y[::nskip,::nskip],umean[::nskip,::nskip],vmean[::nskip,::nskip],length=3,barbcolor='k',flagcolor='r',linewidth=0.5,ax=ax1)
            plt.title('member %02i' % (nanal),fontsize=10)

            # mslp plot
            if isMSL:
              plt.figure(2)
              ax2 = plt.subplot(4,3,nanal)           ## 4 rows x 3 cols
              map.drawcoastlines(ax=ax2)
              map.fillcontinents(ax=ax2,color='0.8',lake_color='0.8')
              map.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=8,ax=ax2)
              map.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=8,ax=ax2)
              xg,yg = map(lons,lats)
              xg2,yg2 = map2(lons,lats)
              CS1=map.contour(xg,yg,mslpmean,levels=np.arange(950,1030,2),colors='k',linewdith=0.5,ax=ax2)
              CS2=map.contour(xg,yg,mslpmean,levels=[1006],colors='b',linewidth=0.5,ax=ax2)
              print mslpmean.min(), mslpmean.max(), mslpmean.shape
              #plotlows(mslpmean)
              plt.title('member %02i' % (nanal),fontsize=10)

            # pwat plot
            if isPWAT:
              plt.figure(4)
              ax4 = plt.subplot(4,3,nanal)           ## 4 rows x 3 cols
              map2.drawcoastlines(ax=ax4)
              map2.fillcontinents(ax=ax4,color='0.8',lake_color='0.8')
              map2.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=8,ax=ax4)
              map2.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=8,ax=ax4)
              if not isMSL:
                  xg2,yg2 = map2(lons,lats)
              CS4=map2.contourf(xg2,yg2,pwat,np.arange(15,66,5),cmap=plt.cm.jet,extend='both',ax=ax4)
              plt.title('member %02i' % (nanal),fontsize=10)

            #precip plot
            if fhr > '000' and isPrecip:
               plt.figure(3)
               ax3 = plt.subplot(4,3,nanal)           ## 4 rows x 3 cols
               map.drawcoastlines(ax=ax3)
               #map.fillcontinents()
               map.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=8,ax=ax3)
               map.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=8,ax=ax3)
               if not isMSL:
                  xg,yg = map(lons,lats)
               ### Susan's levels
               ### clevs=[0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
               ### Phil's levels
               clevs=[0,0.1,1,2,4,6,8,10,12.5,15,20,30,40,50,70,100,150,200,250,300,400,500]
               CS3=map.contourf(xg,yg,precip,clevs,cmap=cm.s3pcpn,extend='max',ax=ax3)
               plt.title('member %02i' % (nanal),fontsize=10)

            grbs.close()

        ### end for nanal

        # settings for plots
        fig1.text(0.5,0.96,'FIM-40km 10-m wind ens %s-hr fcst for %s from %s' %(fhr,strmid,date),
                horizontalalignment='center',fontsize=16)
        fig1.text(0.5,0.075,'NOAA/ESRL Global Systems Division',
                horizontalalignment='center',fontsize=12)
        fig2.text(0.5,0.96,'FIM-40km MSLP ens %s-hr fcst for %s from %s' %(fhr,strmid,date),
                horizontalalignment='center',fontsize=16)
        fig2.text(0.5,0.075,'NOAA/ESRL Global Systems Division',
                horizontalalignment='center',fontsize=12)
        fig4.text(0.5,0.96,'FIM-40km PWAT ens %s-hr fcst for %s from %s' %(fhr,strmid,date),
                horizontalalignment='center',fontsize=16)
        fig4.text(0.5,0.075,'NOAA/ESRL Global Systems Division',
                horizontalalignment='center',fontsize=12)
        plt.figure(1); plt.axes(ax1)
        cax1 = fig1.add_axes([0.1, 0.035, 0.8, 0.02])
        cb = plt.colorbar(CS0,cax=cax1, orientation='horizontal')
        cb.ax.set_xlabel('wind speed (knots)')

        # save plots to disk
        dirName=directories[domain]
        print '***** dirName: ',dirName
        plt.savefig('%s/%s/10mwind_f%03i.png'%\
                     (fimensPath,dirName,int(fhr)))
        plt.close(1)

        if isMSL:
            plt.figure(2); plt.axes(ax2)
                        # (date,strmid,int(fhr)))
            plt.savefig('%s/%s/mslp_f%03i.png'%\
                          (fimensPath,dirName,int(fhr)))
            plt.close(2)

        if isPWAT:
            plt.figure(4); plt.axes(ax4)
            cax4 = fig4.add_axes([0.1, 0.035, 0.8, 0.02])
            cb = plt.colorbar(CS4,cax=cax4, orientation='horizontal')
            cb.ax.set_xlabel('Precipitable water (mm)')
            plt.savefig('%s/%s/pwat_f%03i.png'%\
                         (fimensPath,dirName,int(fhr)))
            plt.close(4)

        if fhr > '000' and isPrecip:
            plt.figure(3); plt.axes(ax3)
            fig3.text(0.5,0.96,'FIM-40km Precip (mm) ens %s-%s hr fcst for %s from %s' %\
            (repr(int(fhr)-6),fhr,strmid,date),horizontalalignment='center',fontsize=16)
            fig3.text(0.5,0.075,'NOAA/ESRL Global Systems Division',
                horizontalalignment='center',fontsize=12)
            #left, bottom, width, height] 
            cax = plt.axes([0.1, 0.035, 0.8, 0.02])
            cb = plt.colorbar(CS3,cax=cax, orientation='horizontal')
            cb.ax.set_xlabel('Precipitation (mm)')
            # plt.savefig('precip%s_%s_gfsenkf_ens_f%03i.png' %\
                    # (date,strmid,int(fhr)))
            plt.savefig('%s/%s/precip_f%03i.png' %\
                    (fimensPath,dirName,int(fhr)))
            plt.close(3)

    #plt.show()
    #raise SystemExit
