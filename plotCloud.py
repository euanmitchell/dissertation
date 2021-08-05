######################################
# Script to take the .pts file output
# from lasPoints and make vertical
# section plots of ground and canopy
######################################

import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj, transform

inProj=Proj(init="epsg:4326")
outProj=Proj(init="epsg:32616")
lon=-84.018017
lat=10.428154
lon,lat=transform(inProj,outProj,lon,lat)
print('lon',lon)
print('lat',lat)

inGround='../data/gridClouds/825580.1152770.ground.pts'
inCanopy='../data/gridClouds/825580.1152770.can.pts'
inRaw=''

canopyData=np.loadtxt(inCanopy, usecols=(0,1,2,4),comments='#') #0=x, 1=y, 2=z, 4=scanAng
groundData=np.loadtxt(inGround, usecols=(0,1,2,4),comments='#')
print(groundData.shape[0])

plt.plot(canopyData[:,0],canopyData[:,2],'x')
plt.plot(groundData[:,0],groundData[:,2],'x')
plt.savefig('../data/gridClouds/825580.1152770.lon.png')
plt.close()
plt.clf()
#plt.show()

plt.plot(canopyData[:,1],canopyData[:,2],'x')
plt.plot(groundData[:,1],groundData[:,2],'x')
plt.savefig('../data/gridClouds/825580.1152770.lat.png')
plt.close()
plt.clf()
#plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''inGround='clouds/19970600200347295.ground.pts'
inCanopy='clouds/19970600200347295.can.pts'
inRaw=''

canopyData=np.loadtxt(inCanopy, usecols=(0,1,2,4),comments='#') #0=x, 1=y, 2=z, 4=scanAng
groundData=np.loadtxt(inGround, usecols=(0,1,2,4),comments='#')
print(groundData.shape[0])

plt.plot(canopyData[:,0],canopyData[:,2],'x')
plt.plot(groundData[:,0],groundData[:,2],'x')
plt.savefig('clouds/19970600200347295lon.png')
plt.close()
plt.clf()

plt.plot(canopyData[:,1],canopyData[:,2],'x')
plt.plot(groundData[:,1],groundData[:,2],'x')
plt.savefig('clouds/19970600200347295lat.png')
plt.close()
plt.clf()'''
