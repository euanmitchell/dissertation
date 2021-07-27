###################################
# Read through a "groundBounds" ALS
# file and find the ultimate bounds
# of the dataset.
###################################

import numpy as np
from pyproj import Proj, transform

inFile='/exports/csce/datastore/geos/groups/3d_env/data/ALS/lists/ground_bounds/alsBounds.nasa_laselva_2009.txt'

minXArray=np.empty(0,dtype=float)
minYArray=np.empty(0,dtype=float)
maxXArray=np.empty(0,dtype=float)
maxYArray=np.empty(0,dtype=float)

minX=10000000.0
minY=10000000.0
maxX=-10000000.0
maxY=-10000000.0

counter=0

with open (inFile,'r') as f:
    for line in f:
        lowX=line.split()[1]
        minXArray=np.append(minXArray,float(lowX))
        lowY=line.split()[2]
        minYArray=np.append(minYArray,float(lowY))
        hiX=line.split()[4]
        maxXArray=np.append(maxXArray,float(hiX))
        hiY=line.split()[5]
        maxYArray=np.append(maxYArray,float(hiY))
        if float(lowX) < minX:
            minX=float(lowX)
        if float(lowY) < minY:
            minY=float(lowY)
        if float(hiX) > maxX:
            maxX=float(hiX)
        if float(hiY) > maxY:
            maxY=float(hiY)
        counter+=1
    print('Lines read',counter)
    print('minX minY maxX maxY',minX,minY,maxX,maxY)
    print('array values',np.amin(minXArray),np.amin(minYArray),np.amax(maxXArray),np.amax(maxYArray))

latBounds=(minY,maxY)
lonBounds=(minX,maxX)

inProj=Proj(init="epsg:32616")
outProj=Proj(init="epsg:4326")
lons,lats=transform(inProj,outProj,lonBounds,latBounds)
print(lons,lats)
