import matplotlib
matplotlib.use('Agg')
import pygrib, sys, os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

date = sys.argv[1] # read in anal date YYYYMMDDHH on command line

fhr = sys.argv[2]
datapath = sys.argv[3]
datapathGefs = sys.argv[4]
yyjdyhh  = sys.argv[5]
fimensPath = sys.argv[6]
tcvitalsPath = sys.argv[7]
nmembers = sys.argv[8]
t2_str = sys.argv[9]
nmembers = int(nmembers)
nanals = nmembers*2
#JKH  for now since can't read GRIB2 files
nanals = nmembers

t2 = int(t2_str)
print 'in plotprob.py: date: ' , date
print 'in plotprob.py: fhr: ' , fhr
print 'in plotprob.py: datapath: ' , datapath
print 'in plotprob.py: datapathGefs: ' , datapathGefs
print 'in plotprob.py: yyjdyhh: ' , yyjdyhh
print 'in plotprob.py: fimensPath: ' , fimensPath
print 'in plotprob.py: tcvitalsPath: ' , tcvitalsPath
print 'in plotprob.py: nmembers: ' , nmembers
print 'in plotprob.py: nanals: ' , nanals
print 'in plotprob.py: t2: ' , t2

# two prob thresholds
thresh  = 34 # tropical storm (kts)
#thresh2 = 64 # hurricane

# path to ens grib output
expname = 'FIMENS_realtime2014'

# forecast hours
fhrs = range(0,t2,6)

# initialize some vars
prob = None; mapdict = {}
mapdict['Atlantic']=None
mapdict['West_Pacific']=None
mapdict['East_Pacific']=None

ifhr = 0
while (ifhr <= t2):
    # for ifhr in fhrs: # loop over fcst hours
    print 'ifhr: ',ifhr
    fhr = '%03i' % ifhr
    # read in ensemble.
    for nanal in range(1,nanals+1):
        # open grib file for nanal'th ens member, extract 10m winds.

        if nanal<=nmembers:
            filename = os.path.join("%s/post_%02i/fim/NAT/grib1/%s000%03i"%(datapath,nanal,yyjdyhh,ifhr))
        else:
            expname = 'hybda_realtime2014'
            dir = datapathGefs
            filename = os.path.join(dir,"%s000%03i_mem%02i"%(yyjdyhh,ifhr,nanal))
            
        print 'nanal: ',nanal
        print 'filename: ',filename
        
        grbs = pygrib.open(filename)
        grbu = grbs.select(shortName='10u')[0]
        grbv = grbs.select(shortName='10v')[0]
        grbu_val = grbu.values
        grbv_val = grbv.values
        grbu_val = grbu_val[::-1,:]         # swap data
        grbv_val = grbv_val[::-1,:]              
        u10m = 1.943844*grbu_val            # mps to knots
        v10m = 1.943844*grbv_val
        # wind speed (converted to 1d array)
        spd = np.sqrt(u10m**2+v10m**2).ravel()
        # calculate prob
        if prob is None: # only define these once.
            lats, lons = grbu.latlons()
            lats1d = -lats[:,0]
            lons1d = (360./lons.shape[1])*np.arange(lons.shape[1])
            lons, lats = np.meshgrid(lons1d,lats1d)
	    nlats, nlons = lats.shape
            flag = np.zeros((nanals,nlats*nlons), np.bool)
            prob = np.zeros(nlats*nlons, np.float)
            flag2 = np.zeros((nanals,nlats*nlons), np.bool)
            #prob2 = np.zeros(nlats*nlons, np.float)
        pts = np.argwhere(np.logical_and(spd>thresh,flag[nanal-1]==False)) 
        #pts2 = np.argwhere(np.logical_and(spd>thresh2,flag2[nanal-1]==False)) 
        prob[pts]=prob[pts]+1.0 # ts prob
        # each ens member only gets counted once per grid point.
        flag[nanal-1,pts]=True 
        #prob2[pts2]=prob2[pts2]+1.0 # hurr prob
        #flag2[nanal-1,pts2]=True
        print nanal,spd.min(),spd.max(),len(pts)#,len(pts2)
        grbs.close() # close grib file.

    # make a copy of prob array, normalize and reshape to 2D
    tcprob = 100.*prob.copy()/nanals
    tcprob.shape = (nlats,nlons)
    print ifhr,tcprob.min(),tcprob.max()

    # make plots
    # if ifhr%6==0:
    if ifhr==t2:
        for basin in ['Atlantic','East_Pacific','West_Pacific']:

            # setup Basemap instances for each region
            # on first pass thru, cache in dict.
            if mapdict[basin]==None:
                if basin == 'Atlantic':
                    lat_0 = 22.5; lon_0=-60.
                    width=10000.e3; height=5000.e3
                elif basin == 'East_Pacific':
                    lat_0 = 22.5; lon_0=360.-130
                    width=10000.e3; height=5000.e3
                elif basin == 'West_Pacific':
                    lat_0 = 22.5; lon_0=135.
                    width=10000.e3; height=5000.e3
                map =\
                Basemap(projection='lcc',lon_0=lon_0,lat_0=lat_0,width=width,height=height,resolution='l')
                x, y = map(lons,lats)
                mapdict[basin]=map,x,y
            else:
                map,x,y = mapdict[basin]
    
            # ts prob plot.
            # create figure (size in inches).
            fig = plt.figure(figsize=(10,8))
            # set contour levs
            clevs = np.arange(5,101,5)
            # set colormap
            cmap = plt.cm.spectral_r
            # set lat and lon grid lines
            lonlines = np.arange(0,360,15)
            latlines = np.arange(-90,90,10)
            # draw geography
            map.drawcoastlines()
            map.drawcountries()
            map.drawstates()
            # draw lat/lon grid lines
            map.drawmeridians(lonlines,labels=[0,0,0,1],fontsize=12)
            map.drawparallels(latlines,labels=[1,0,0,0],fontsize=12)
            # draw filled contours
            CS=map.contourf(x,y,tcprob,clevs,cmap=cmap)
            # draw colorbar
            cb = plt.colorbar(CS, orientation='horizontal',format='%g')
            cb.ax.set_xlabel('probability of wind speed > %s kts (percent)' % thresh)
            # draw title
            plt.title('FIM-40km (1-10) and GEFS (11-20) Ens \nTrop Storm Force Wind Probabilities %s hours from %s' %\
                    (ifhr,date),fontsize=12)
            # save to png

            outputFile= '%s/%s/tcprob_f%03i.png' % (fimensPath,basin,int(fhr))
            print 'output file: ',outputFile
            plt.savefig('%s/%s/tcprob_f%03i.png'%\
                    (fimensPath,basin,int(fhr)))
            # close figure
            plt.close()
    ifhr = ifhr + 6
print 'plotprob.py finished!'
