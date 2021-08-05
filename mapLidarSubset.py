###############################################
# Script to generate commands for mapLidar
# to generate DSM and DTM tiles from ALS data.
###############################################

import h5py
import numpy as np
from pyproj import Proj, transform
import os
import glob

###############################################

latBounds=(45.775,45.892)
lonBounds=(-121.775,-122.103)

inProj=Proj(init="epsg:4326")
outProj=Proj(init="epsg:32610")
lons,lats=transform(inProj,outProj,lonBounds,latBounds)

minLat=np.min(lats)
maxLat=np.max(lats)
minLon=np.min(lons)
maxLon=np.max(lons)

stepLon=(maxLon-minLon)/5
stepLat=(maxLat-minLat)/5

with open('shell/mapWREF.sh', 'a') as outf:
    outf.write('#!/usr/bin/env bash\n\n')

i=0
for x in np.arange(minLon,maxLon,stepLon):
    for y in np.arange(minLat,maxLat,stepLat):
        with open('shell/mapWREF.sh', 'a') as outf:
            outf.write('mapLidar -output ../data/wref/dems/DSM{0} -inList ../data/groundBounds.neon_wref.txt -res 10 -height -float -epsg 32610 -bounds {1} {2} {3} {4}\n'.format(i,x,y,x+stepLon,y+stepLat))
        i+=1
