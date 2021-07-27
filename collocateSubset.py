###############################################
# Script to determine bounds of input GEDI
# h5py files and generate commands for
# collocateWaves for spatial subsets of the
# input data.
###############################################

import h5py
import numpy as np
from pyproj import Proj, transform
import os
import glob

###############################################

dir="/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/subset/2009/"
inputs='*.h5'
path=os.path.join(dir,inputs)
files=glob.glob(path)

beamlist=['BEAM0101','BEAM1000','BEAM0010','BEAM0000','BEAM0011','BEAM0110','BEAM0001','BEAM1011']
inProj=Proj(init="epsg:4326")
outProj=Proj(init="epsg:32616")

for file in files:
    outFile=file.split('/')
    outFile=outFile[-1].split('_')
    outFile='shell/laselvaCollo/'+outFile[3]+'.sh'
    alsList='/exports/csce/datastore/geos/groups/3d_env/data/ALS/lists/ground_bounds/alsBounds.nasa_laselva_2009.txt'
    waveforms=0
    f=h5py.File(file,'r')
    print('Reading file',file)
    lats=np.empty(0)
    lons=np.empty(0)
    for beam in beamlist:
        try:
            lats=np.append(lats,f[beam]['geolocation']['latitude_lastbin'])
            lons=np.append(lons,f[beam]['geolocation']['longitude_lastbin'])
        except:
            print('No data in',beam)
    lons,lats=transform(inProj,outProj,lons,lats)
    minLat=np.min(lats)
    maxLat=np.max(lats)
    minLon=np.min(lons)
    maxLon=np.max(lons)

    if lons.shape[0] < 2000:
        steps=10
    else:
        steps=20

    stepLon=(maxLon-minLon)/steps

    with open(outFile, 'a') as outf:
        outf.write('#!/usr/bin/env bash\n\necho "Processing file ' + file + '"\n\n')

    for i in range(steps):
        waves=file.split('/')
        wave=waves[-1]
        wave='sim'+str(i)+wave[6:]
        output=wave[:-2]+'correl.txt'
        with open(outFile,'a') as outf:
            outf.write('''collocateWaves -listAls {0} -gedi {1} -solveCofG -aEPSG 32616 -geoError 20 5 -writeWaves /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/sim/2009/{2} -output /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/sim/2009/{3} -fixFsig -minSense 0.9 -readPulse ~/Dissertation/code/pulses/meanPulse.BEAM1000.filt -bounds {4} {5} {6} {7}\n\n'''.format(alsList, file, wave, output, minLon+(stepLon*i), minLat, minLon+(stepLon*(i+1)), maxLat))
            #outf.write('''collocateWaves -listAls ~/Dissertation/data/{0} -gedi {1} -solveCofG -aEPSG 32617 -geoError 40 10 -writeWaves /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/sim/{2} -output /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/sim/{3} -fixFsig -minSense 0.9 -maxScanAng 30 -decimate 100 -readPulse ~/Dissertation/code/pulses/meanPulse.BEAM1000.filt -bounds {4} {5} {6} {7}\n\n'''.format(alsList, file, wave, output, minLon+(stepLon*i),
        #print('Longitude range:',minLon+(stepLon*i),minLon+(stepLon*(i+1)))

    with open(outFile, 'a') as outf:
        outf.write('echo "All done"')
