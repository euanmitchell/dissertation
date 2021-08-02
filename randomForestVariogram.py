########################################
# Script to run Random Forest to
# predict canopy height from Sentinel-2
# for real La Selva GEDI data.
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
from sklearn.metrics import mean_squared_error
#import skgstat as skg
from math import sqrt
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
    p.add_argument('--inFile', dest='inFile', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/grid20/cov99wid0.5/metricAll.txt',
        help=('The path to a single simulated data file.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/grid20/cov99wid0.5/metricAll.txt'))
    p.add_argument('--inTiff', dest='inTiff', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/laselva20mComposite.tif',
        help=('The path to a single tiff image.\nDefault is '))
    p.add_argument('--tiff', dest='tiff', action='store_true', default=False,
        help=('Read the specified tiff and make relevant plots.'))
    p.add_argument('--simPlots', dest='simPlots', action='store_true', default=False,
        help=('Make histogram and box plots from the simulated data.'))
    p.add_argument('--beams', dest='beams', type=str, default='all',
        help=('Which beams to use: "all", "power", or "coverage".\nDefault is all'))
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

    def readMetric(self,inFile):
        '''Read the metric data file'''

        print('Reading ALS metric file',inFile)

        self.metricData=np.loadtxt(inFile, usecols=(1,3,4,5,10,32,95,106,107))
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

        print('metricData shape',self.metricData.shape)
        bad=np.where(self.metricData[:,0]==-1000000)
        #print(bad)

        self.use=np.delete(self.metricData,bad,axis=0)
        print('self.use shape',self.use.shape)

    def readTiff(self,filename):
        '''Read in a satellite image'''

        # Read and process NDVI image
        ds=gdal.Open(filename)

        proj=osr.SpatialReference(wkt=ds.GetProjection())
        self.epsg=int(proj.GetAttrValue('Authority',1))
        self.nX=ds.RasterXSize
        self.nY=ds.RasterYSize

        # Get image origin and resolution info
        transform=ds.GetGeoTransform()
        self.xOrigin=transform[0]
        self.yOrigin=transform[3]
        self.xPixel=transform[1]
        self.yPixel=transform[5]

        # Read each raster band into 3D numpy array
        self.tiffData=ds.ReadAsArray(0,0,self.nX,self.nY)

    def tiffExtract(self):
        '''Extract data from raster at footprint locations'''

        self.tiffExtract=np.zeros((self.use.shape[0],14))
        print('tiffExtract shape',self.tiffExtract.shape)

        for i in range(self.use.shape[0]):  #self.use.shape[0]
            tempLon=self.use[i,7]
            tempLat=self.use[i,8]

            xDist=tempLon-self.xOrigin
            yDist=tempLat-self.yOrigin

            xInd=int((tempLon-self.xOrigin)//self.xPixel)
            yInd=int((tempLat-self.yOrigin)//self.yPixel)

            self.tiffExtract[i,:]=self.tiffData[:,yInd,xInd]

    def generateForest(self):
        '''Run Random Forest algorithm with the imported datasets'''

        plt.rcParams['figure.figsize']=(8,6)
        plt.rcParams['xtick.labelsize']=16
        plt.rcParams['ytick.labelsize']=16
        plt.rcParams['axes.labelsize']=18
        plt.rcParams['axes.labelpad']=8.0
        plt.rcParams['axes.titlesize']=20

        # Assign training data as a particular row of cells from the grid
        self.train=np.where((self.use[:,8]==1152610) | (self.use[:,8]==1152630) | (self.use[:,8]==1152650))
        print('train length',len(self.train[0]))

        # Independent variables are Sentinel reflectances and indices
        x_train=self.tiffExtract[self.train[0],:]
        print('Shape of x',x_train.shape)
        # Dependent variable is ALS RH95
        y_train=self.use[self.train[0],6]
        print('Shape of y',y_train.shape)

        distance=np.zeros(49)
        scores=np.zeros(49)

        for i in range(2,51):
            dist=i*140
            min=1152630-dist
            max=1152630+dist
            distance[i-2]=dist

            self.test=np.where((self.use[:,8]==min) | (self.use[:,8]==max))
            print('test length',len(self.test[0]))

            if len(self.test[0]) > 0:
                x_test=self.tiffExtract[self.test[0],:]
                y_test=self.use[self.test[0],6]

                print('Running Random Forest algorithm iteration {} ...'.format(i))
                #x=self.use[:,7:9]
                #y=self.use[:,6]
                #x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.1,random_state=0)
                regressor=RandomForestRegressor(n_estimators=200,random_state=0)
                regressor.fit(x_train,y_train)
                self.y_pred=regressor.predict(x_test)
                print('y_pred shape',self.y_pred.shape)
                print('Random Forest coefficient',regressor.score(x_test,y_test))
                print('Random Forest RMSE',sqrt(mean_squared_error(y_test,self.y_pred)))
                print('Random Forest Bias',np.sum(self.y_pred-y_test)/self.y_pred.shape[0])
                #print('RF Feature Importance',regressor.feature_importances_)
                scores[i-2]=regressor.score(x_test,y_test)

        plt.plot(distance,scores,'o')
        plt.xlabel('Distance (m)')
        plt.ylabel('RF Correlation Coefficient')
        plt.ylim([-0.8,0.8])
        plt.show()

        # Independently fit a linear regression to the RF output
        '''slope,intercept,r,p,se = sps.linregress(y_test,self.y_pred)
        print('Slope: {}, intercept: {}, r: {}, p: {}, std_err: {}'.format(slope,intercept,r,p,se))

        #print('Height minimum:',np.amin(self.y_pred))
        #print('Height maximum:',np.amax(self.y_pred))
        #print('Height mean:',np.mean(self.y_pred))
        #print('Height std. dev.:',np.std(self.y_pred))

        # The ground residual as defined in alsComp.py
        #self.residual=(self.gediExtract[self.test[0],0] - self.use[self.test[0],0])

        # The difference between the RF prediction of canopy height and the GEDI L2A height value
        #self.RFresidual=(self.y_pred-y_test)'''

    def makeVariogram(self):
        '''Compute and plot semi-variogram of RH95 values'''

        '''variogram=skg.Variogram(coordinates=self.use[::10,7:9],values=self.use[::10,6])
        print(variogram)
        bins,values=variogram.get_empirical()
        plt.plot(bins,values,'o')
        plt.show()'''

        plt.rcParams['figure.figsize']=(8,6)
        plt.rcParams['xtick.labelsize']=16
        plt.rcParams['ytick.labelsize']=16
        plt.rcParams['axes.labelsize']=18
        plt.rcParams['axes.labelpad']=8.0
        plt.rcParams['axes.titlesize']=20

        maxDist=10000
        res=20
        nBins=int(maxDist/res+1)
        meanVar=np.zeros(nBins,dtype=float)
        contN=np.zeros(nBins,dtype=int)

        for i in range(0,self.use.shape[0],10):
            dists=np.sqrt((self.use[i,7]-self.use[::10,7])**2+(self.use[i,8]-self.use[::10,8])**2)
            diffs=self.use[::10,6]-self.use[i,6]

            bins=np.array(dists//res,dtype=int)

            for j in range(0,bins.shape[0]):
                if (j==i) | (bins[j]<0) | (bins[j]>=contN.shape[0]):
                    continue
                meanVar[bins[j]]=meanVar[bins[j]]+diffs[j]**2
                contN[bins[j]]=contN[bins[j]]+1

        meanVar[contN>0]=np.sqrt(meanVar[contN>0]/contN[contN>0])

        dists=np.arange(0,maxDist,res)
        plt.plot(dists[1:],meanVar[1:-1])
        plt.xlabel('Range (m)')
        plt.ylabel('RH95 Semi-Variance (m)')
        plt.show()

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    comp=gridMetrics()
    comp.readMetric(args.inFile)
    comp.readTiff(args.inTiff)
    comp.tiffExtract()
    comp.generateForest()
    #comp.makeVariogram()
