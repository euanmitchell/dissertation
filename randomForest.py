########################################
# Script to run Random Forest to
# predict canopy height from Sentinel-2.
########################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
from osgeo import gdal
from pyproj import Proj, transform
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
import scipy.stats as sps
import osr
import argparse
import os
import glob

#######################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify input ALS and GEDI data files and programme control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--inDir', dest='inDir', type=str, default=' ',
        help=('The path to the directory containing the metric data.\nDefault is not set'))
    p.add_argument('--inFile', dest='inFile', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/grid20/cov99wid0.5/metricAll.txt',
        help=('The path to a single simulated data file.\nDefault is not set'))
    p.add_argument('--inTiff', dest='inTiff', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/laselva20mComposite.tif',
        help=('The path to a single tiff image.\nDefault is '))
    p.add_argument('--tiff', dest='tiff', action='store_true', default=False,
        help=('Read the specified tiff and make relevant plots.'))
    p.add_argument('--simPlots', dest='simPlots', action='store_true', default=False,
        help=('Make histogram and box plots from the simulated data.'))
    cmdargs=p.parse_args()
    return cmdargs

#######################################

class gridMetrics():
    '''
    A class to read output of gediRat simulations and imagery
    and make comparison plots.
    '''

    def __init__(self):
        '''Class initialiser'''

    def readTiff(self,filename):
        '''Read in a satellite image'''

        # Read and process NDVI image
        ds=gdal.Open(filename)

        proj=osr.SpatialReference(wkt=ds.GetProjection())
        self.epsg=int(proj.GetAttrValue('Authority',1))
        print('Raster EPSG:',self.epsg)
        self.nX=ds.RasterXSize
        self.nY=ds.RasterYSize
        print('Raster size:',self.nX,self.nY)

        # Get image origin and resolution info
        transform=ds.GetGeoTransform()
        self.xOrigin=transform[0]
        self.yOrigin=transform[3]
        self.xPixel=transform[1]
        self.yPixel=transform[5]
        print('Raster origin:',self.xOrigin,self.yOrigin)
        print('Raster resolution:',self.xPixel,self.yPixel)

        # Read each raster band into 3D numpy array
        self.tiffData=ds.ReadAsArray(0,0,self.nX,self.nY)
        print('Array shape:',self.tiffData.shape)
        print('Array data:',self.tiffData[:,0,0])
        #print(self.tiffData[0,0])
        #print(self.tiffData[0,10979])

    def readMetric(self,inFile):

        print('Reading ALS metric file',inFile)

        self.data=np.loadtxt(inFile, usecols=(1,3,4,5,10,32,95,106,107,112))
        # Each column is a different metric as follows
        #0 1 = 'true ground' (from ALS)
        #1 3 = 'ground slope' (from ALS)
        #2 4 = 'ALS cover' (from ALS)
        #3 5 = 'gHeight' (from waveform)
        #4 10 = 'cover' (from waveform)
        #5 32 = 'rhGauss 95' (waveform RH95)
        #6 95 = 'rhReal 95' (ALS RH95)
        #7 106 = 'lon'
        #8 107 = 'lat'
        #9 112 = 'blairSense' (from ?)

        print(self.data.shape)
        bad=np.where(self.data[:,0]==-1000000)
        #print(bad)

        self.use=np.delete(self.data,bad,axis=0)
        print(self.use.shape)

        self.residual = (self.use[:,3] - self.use[:,0])

        plt.plot(self.use[:,6],self.use[:,5],'o',markersize=1)
        plt.plot([0,60], [0,60],ls='--',color='r')
        plt.xlim([0,60])
        plt.ylim([0,60])
        plt.title('ALS vs Waveform Height')
        plt.xlabel('ALS RH95 (m)')
        plt.ylabel('Waveform RH95 (m)')
        #plt.savefig('../data/laselvaFigures/randomForest1.png')
        #plt.close()
        #plt.clf()
        plt.show()

    def tiffExtract(self):
        '''Extract data from raster at footprint locations'''

        self.tiffExtract=np.zeros((self.residual.shape[0],14))
        print(self.tiffExtract.shape)

        for i in range(self.use.shape[0]):  #self.use.shape[0]
            tempLon=self.use[i,7]
            #tempLon=self.lons[i]
            #tempLon=826030
            tempLat=self.use[i,8]
            #tempLat=self.lats[i]
            #tempLat=1149590

            xDist=tempLon-self.xOrigin
            yDist=tempLat-self.yOrigin
            #print(xDist,yDist)

            xInd=int((tempLon-self.xOrigin)//self.xPixel)
            yInd=int((tempLat-self.yOrigin)//self.yPixel)
            #print(xInd,yInd)

            self.tiffExtract[i,:]=self.tiffData[:,yInd,xInd]
            #print('Tiff value',self.tiffData[:,yInd,xInd])

        x=self.tiffExtract[:,:]
        print('Shape of x',x.shape)
        y=self.use[:,6]
        print('Shape of y',y.shape)

        x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.1,random_state=0)

        print('Running Random Forest algorithm ...')
        regressor=RandomForestRegressor(n_estimators=10,random_state=0)
        regressor.fit(x_train,y_train)
        y_pred=regressor.predict(x_test)
        print('Random Forest coefficient',regressor.score(x_train,y_train))

        plt.plot(y_test,y_pred,'o',markersize=1)
        plt.plot([0,60], [0,60],ls='--',color='r')
        plt.xlim([0,60])
        plt.ylim([0,60])
        plt.title('Random Forest Height Prediction')
        plt.xlabel('ALS RH95 (m)')
        plt.ylabel('RF Predicted Height (m)')
        plt.savefig('../data/laselvaFigures/randomForest1.png')
        plt.close()
        plt.clf()
        #plt.show()

        # Independently fit a linear regression to the RF output
        slope,intercept,r,p,se = sps.linregress(y_test,y_pred)
        print('Slope: {}, intercept: {}, r: {}, p: {}, std_err: {}'.format(slope,intercept,r,p,se))

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    comp=gridMetrics()
    comp.readMetric(args.inFile)
    comp.readTiff(args.inTiff)
    comp.tiffExtract()
