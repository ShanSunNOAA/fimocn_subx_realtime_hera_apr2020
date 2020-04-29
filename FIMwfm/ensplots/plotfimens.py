import matplotlib
matplotlib.use('Agg')
from dateutils import splitdate,dayofyear, dateshift
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid, cm
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import pygrib
import sys, os, glob, collections
from scipy.ndimage.filters import minimum_filter, maximum_filter

def plot_tracks(map,fcstracks,obtracks,ifhr,ifhrinc):
    for storm in obtracks.keys():
        xfpos = []; yfpos = []
        xopos = []; yopos = []
        for ihr in range(0,ifhr+1,ifhrinc):
            try:
                ftrackdat = fcstracks[storm][ihr]
                xf,yf = map(ftrackdat[1],ftrackdat[0])
                xfpos.append(xf); yfpos.append(yf)
            except:
                pass
            try:
                otrackdat = obtracks[storm][ihr]
                xo,yo = map(otrackdat[1],otrackdat[0])
                xopos.append(xo); yopos.append(yo)
            except:
                pass
        if xfpos: 
            map.plot(xfpos,yfpos,'bo-',ms=3)
        if xopos: 
            map.plot(xopos,yopos,'ro-',ms=3)

def extrema(mat,mode='wrap',window=10):
    """find the indices of local extrema (min and max)
    in the input array."""
    mn = minimum_filter(mat, size=window, mode=mode)
    mx = maximum_filter(mat, size=window, mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    # (mat == mn) true if pixel is equal to the local in
    # Return the indices of the maxima, minima
    return np.nonzero(mat == mn), np.nonzero(mat == mx)

def plotlows(mslpgrid):
    # the window parameter controls the number of highs and lows detected.
    # (higher value, fewer highs and lows)
    local_min, local_max = extrema(mslpgrid, mode='wrap', window=40)
    xlows = x[local_min]; xhighs = x[local_max]
    ylows = y[local_min]; yhighs = y[local_max]
    lowvals = mslpgrid[local_min]; highvals = mslpgrid[local_max]
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

date = sys.argv[1]
ifhr1 = int(sys.argv[2])
ifhr2 = int(sys.argv[3])
ifhrinc = int(sys.argv[4])
if len(sys.argv) > 5:
    run = sys.argv[5]
else:
    run = 'FIMY'
expname = os.getenv('EXPT')
pngdir = os.getenv('PNGDIR')
if expname is None:
    expname = 'gfsenkf_t574'
    pngdir = '/lfs1/projects/gfsenkf/hurrplots/%s' % date
hr = date[8:10]
yyyy,mm,dd,hh = splitdate(date)
julday = dayofyear(yyyy,mm,dd)+1
datapath = '/lfs1/projects/gfsenkf/tcvitals'
globstring = datapath+'/reftrk*%s*' % date
print globstring
reftrks = glob.glob(globstring)
print reftrks

reftrks.insert(0,'WPAC')
reftrks.insert(0,'EPAC')
reftrks.insert(0,'ATL')

# just use the three basins (no domains centered on storms)
#reftrks = ['ATL','EPAC','WPAC']

# read forecast tracks.
def makehash():
    return collections.defaultdict(makehash)
datapath = "/lfs1/projects/fim/whitaker/%s/%s" % (expname,date)
if run == 'FIMY':
    trackerfile = os.path.join(datapath,"control/track.%s.GE00" % date)
else:
    trackerfile = os.path.join(datapath,"control/track.%s.AVNO" % date)
print trackerfile
fcstracks=makehash()
try:
    for line in open(trackerfile,'r'):
        linesplit = line.split(',')
# WP, 04, 2011052500, 03, GE00, 000, 124N, 1281E,
        basin = linesplit[0]
        storm = linesplit[1]
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
except:
    pass
#print 'fcst tracks'
#print fcstracks

# read observed tracks
obtracks = makehash()
for fhour in range(0,168,6):
    datev = dateshift(date,fhour)
    try:
        for line in open('/lfs1/projects/gfsenkf/tcvitals/tcvitals.%s.txt' % datev):
            linesplit = line.split()
            storm = linesplit[1].lstrip()
            lat = linesplit[5]
            lon = linesplit[6]
            if lat[-1] == 'N':
                lat = 0.1*float(lat[:-1])
            else:
                lat = -0.1*float(lat[:-1])
            if lon[-1] == 'E':
                lon = 0.1*float(lon[:-1])
            else:
                lon = -0.1*float(lon[:-1])
            minp = int(linesplit[9])
            obtracks[storm][fhour] = (lat,lon,minp)
    except:
        continue
#print 'observed tracks'
#print obtracks

fig = None
for reftrk in reftrks:
    if reftrk == 'ATL':
        strmid = reftrk
        lat_0 = 20.; lon_0=-65.
        width=7500.e3; height=4500.e3
    elif reftrk == 'WPAC':
        strmid = reftrk
        lat_0 = 20.; lon_0=135.
        width=7500.e3; height=4500.e3
    elif reftrk == 'EPAC':
        strmid = reftrk
        lat_0 = 20.; lon_0=-125
        width=7500.e3; height=4500.e3
    else:
        width=5000.e3; height=3000.e3
        for nline,line in enumerate(open(reftrk)):
            if not nline:
               strmid = line.split()[3][:3]
            else:
               fhrin = line.split()[0]
               lat_0 = float(line.split()[1])
               lon_0 = float(line.split()[2])
               if int(fhr) == int(fhrin):
                   break
               if int(fhrin) > 0 and int(fhrin) > int(fhr):
                   lon_0 = 0.5*(lon_0 + lon_0_prev)
                   lat_0 = 0.5*(lat_0 + lat_0_prev)
                   break
               lon_0_prev = lon_0
               lat_0_prev = lat_0
    print strmid, lat_0, lon_0
    datapath = "/lfs1/projects/fim/whitaker/%s/%s" % (expname,date)
    map =\
    Basemap(projection='stere',lon_0=lon_0,lat_0=lat_0,width=width,height=height,resolution='l')
    xx, yy = map(lon_0, lat_0)
    for ifhr in range(ifhr1,ifhr2+1,ifhrinc):
        if ifhr < 100:
            fhr = '%02i' % ifhr
        else:
            fhr = repr(ifhr)
    
        if run == 'FIMY':
            filename = os.path.join(datapath,"control/pgb_gfscntl_%s_fhr%s" %
                    (date,fhr))
            filename = os.path.join(datapath,"FIMY/fim_C/%s%03i%s000%03i" %\
                    (date[2:4],julday,hr,int(fhr)))
            #filename =\	    
        #os.path.join(datapath,"control/GE00_d%s_%s.grb" % (date[0:4],date))
        else:
            filename=\
            '/public/data/grids/gfs/0p5deg/grib2/%s%03i%s000%03i' %\
        (date[2:4],julday,hr,int(fhr))
        print filename
        try:
            grbs = pygrib.open(filename)
        except: # file missing, skip to next forecast time.
            continue
        grbmslp=grbs.select(shortName=['msl','prmsl'])[0]
        grbu=grbs.select(shortName='10u')[0]
        grbv=grbs.select(shortName='10v')[0]
        if ifhr:
            grbapcp = grbs.select(shortName='tp')[0]
            precip = grbapcp.values
        ugrid = grbu.values
        vgrid = grbv.values
        mslpgrid = 0.01*grbmslp.values
        lats, lons = grbmslp.latlons()
        lats1d = lats[:,0]
        lons1d = lons[0,:]
        nx = int((map.xmax-map.xmin)/50000.)+1; ny = int((map.ymax-map.ymin)/50000.)+1
        ugrid, lons1 = addcyclic(ugrid, lons1d)
        ugrid,lons2 = shiftgrid(180.,ugrid,lons1,start=False)
        vgrid, lons1 = addcyclic(vgrid, lons1d)
        vgrid,lons2 = shiftgrid(180.,vgrid,lons1,start=False)
        ugrid, vgrid, xx, yy = map.transform_vector(ugrid,vgrid,lons2,lats1d,nx,ny,returnxy=True)
        # convert mps to knots.
        ugrid = 1.943844*ugrid
        vgrid = 1.943844*vgrid
        spd = np.sqrt(ugrid**2 + vgrid**2)
        x,y = map(lons,lats)

        # mslp plot
        fig = plt.figure(figsize=(10,8))
        map.drawcoastlines()
        map.fillcontinents()
        map.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=12)
        map.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=12)
        CS=map.contour(x,y,mslpgrid,levels=np.arange(950,1050,2),colors='k',linewdith=0.5)
        CS=map.contour(x,y,mslpgrid,levels=[1006],colors='b',linewidth=0.5)
        plotlows(mslpgrid)
        plot_tracks(map,fcstracks,obtracks,ifhr,ifhrinc)
        if run == 'gfsenkf':
            plt.title('T574 GFS/EnKF MSLP cntl %s-hr fcst for %s from %s' %\
                    (fhr,strmid,date))
            plt.savefig(os.path.join(pngdir,'mslp%s_%s_gfsenkf_cntrl_f%03i.png'%\
                    (date,strmid,ifhr)))
        else:
            plt.title('T574 GFS/Opnl MSLP cntl %s-hr fcst for %s from %s' %\
                    (fhr,strmid,date))
            plt.savefig(os.path.join(pngdir,'mslp%s_%s_gfsopnl_cntrl_f%03i.png'%\
                    (date,strmid,ifhr)))
        plt.close()

        # 10m wind plot
        fig = plt.figure(figsize=(10,8))
        map.drawcoastlines()
        map.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=12)
        map.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=12)
        CS=map.contourf(xx,yy,spd,np.arange(8,81,8),cmap=cm.GMT_haxby_r,extend='max')
        nskip = 4
        brb = map.barbs(xx[::nskip,::nskip],yy[::nskip,::nskip],ugrid[::nskip,::nskip],vgrid[::nskip,::nskip],length=5,barbcolor='k',flagcolor='r',linewidth=0.5)
        plot_tracks(map,fcstracks,obtracks,ifhr,ifhrinc)
        cb = plt.colorbar(CS, orientation='horizontal')
        cb.ax.set_xlabel('wind speed (knots)')
        if run == 'gfsenkf':
            plt.title('T574 GFS/EnKF 10-m wind cntl %s-hr fcst for %s from %s' %\
                    (fhr,strmid,date))
            plt.savefig(os.path.join(pngdir,'10mwind%s_%s_gfsenkf_cntrl_f%03i.png'%\
                    (date,strmid,ifhr)))
        else:
            plt.title('T574 GFS/Opnl 10-m wind cntl %s-hr fcst for %s from %s' %\
                    (fhr,strmid,date))
            plt.savefig(os.path.join(pngdir,'10mwind%s_%s_gfsopnl_cntrl_f%03i.png'%\
                    (date,strmid,ifhr)))
        plt.close()

        # precip plot
        if ifhr:
            fig = plt.figure(figsize=(10,8))
            map.drawcoastlines()
            map.drawmeridians(np.arange(0,360,15),labels=[0,0,0,1],yoffset=80000,fontsize=12)
            map.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=12)
            clevs=[0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
            CS=map.contourf(x,y,precip,clevs,cmap=cm.s3pcpn,extend='max')
            cb = plt.colorbar(CS, orientation='horizontal',ticks=clevs,\
                    format='%g')
            plot_tracks(map,fcstracks,obtracks,ifhr,ifhrinc)
            if run == 'gfsenkf':
                plt.title('T574 GFS/EnKF Precip cntl %s-%s hr fcst for %s from %s' %\
                        (repr(grbapcp.startStep),repr(grbapcp.endStep),strmid,date))
                plt.savefig(os.path.join(pngdir,'precip%s_%s_gfsenkf_cntrl_f%03i.png'%\
                        (date,strmid,ifhr)))
            else:
                plt.title('T574 GFS/Opnl Precip cntl %s-%s hr fcst for %s from %s' %\
                        (repr(grbapcp.startStep),repr(grbapcp.endStep),strmid,date))
                plt.savefig(os.path.join(pngdir,'precip%s_%s_gfsopnl_cntrl_f%03i.png'%\
                        (date,strmid,ifhr)))
            plt.close()
