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
    p.add_argument('--inDir', dest='inDir', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/',
        help=('The path to the directory containing the GEDI L2A data.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/'))
    p.add_argument('--inFile', dest='inFile', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/metricAll.txt',
        help=('The path to a single simulated data file.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/metricAll.txt'))
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

        def parse(id):
            return str(id).split('.')[2]

        self.waveID=np.loadtxt(inFile, usecols=(0,),converters={0: parse},dtype=np.uint64)
        #print('metric waveID type',self.waveID.dtype)
        #print('First wave ID',self.waveID[0])

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

        self.useWaveID=np.delete(self.waveID,bad,axis=0)
        print('self.useWaveID shape',self.useWaveID.shape)

    def readGedi(self,inDir):
        '''Read the GEDI L2A file'''

        print('Reading GEDI data from', inDir)

        if args.beams == 'all':
            beamlist=['BEAM0101','BEAM1000','BEAM0010','BEAM0000','BEAM0011','BEAM0110','BEAM0001','BEAM1011']
        elif args.beams == 'power':
            beamlist=['BEAM0101','BEAM1000','BEAM0110','BEAM1011']
        elif args.beams == 'coverage':
            beamlist=['BEAM0010','BEAM0000','BEAM0011','BEAM0001']

        inputs='*.h5'
        path=os.path.join(inDir,inputs)
        fileList=glob.glob(path)

        self.shotNumb=np.empty(0,dtype=np.uint64)
        lowestMode=np.empty(0,dtype=float)
        rh95=np.empty(0,dtype=float)
        qualFlag=np.empty(0,dtype=float)
        ground1=np.empty(0,dtype=float)
        ground2=np.empty(0,dtype=float)
        ground3=np.empty(0,dtype=float)
        ground4=np.empty(0,dtype=float)
        ground5=np.empty(0,dtype=float)
        ground6=np.empty(0,dtype=float)

        for file in fileList:
            f=h5py.File(file,'r')
            for beam in beamlist:
                try:
                    self.shotNumb=np.append(self.shotNumb,f[beam]['shot_number'])
                    lowestMode=np.append(lowestMode,f[beam]['elev_lowestmode'])
                    rh95=np.append(rh95,f[beam]['rh'][:,95])
                    qualFlag=np.append(qualFlag,f[beam]['quality_flag'])
                    ground1=np.append(ground1,f[beam]['geolocation']['elev_lowestmode_a1'])
                    ground2=np.append(ground2,f[beam]['geolocation']['elev_lowestmode_a2'])
                    ground3=np.append(ground3,f[beam]['geolocation']['elev_lowestmode_a3'])
                    ground4=np.append(ground4,f[beam]['geolocation']['elev_lowestmode_a4'])
                    ground5=np.append(ground5,f[beam]['geolocation']['elev_lowestmode_a5'])
                    ground6=np.append(ground6,f[beam]['geolocation']['elev_lowestmode_a6'])
                except:
                    continue
                    #print('Empty beam',beam)

        self.gediData=np.stack((lowestMode,rh95,qualFlag,ground1,ground2,ground3,ground4,ground5,ground6),axis=1)
        print('gediData shape',self.gediData.shape)
        #print('first gedi data',self.gediData[0])
        #print('gedi shotnumber type',self.shotNumb.dtype)
        #print('First gedi shotnumber',self.shotNumb[0])

    def readTiff(self,filename):
        '''Read in a satellite image'''

        # Read and process NDVI image
        ds=gdal.Open(filename)

        proj=osr.SpatialReference(wkt=ds.GetProjection())
        self.epsg=int(proj.GetAttrValue('Authority',1))
        #print('Raster EPSG:',self.epsg)
        self.nX=ds.RasterXSize
        self.nY=ds.RasterYSize
        #print('Raster size:',self.nX,self.nY)

        # Get image origin and resolution info
        transform=ds.GetGeoTransform()
        self.xOrigin=transform[0]
        self.yOrigin=transform[3]
        self.xPixel=transform[1]
        self.yPixel=transform[5]
        #print('Raster origin:',self.xOrigin,self.yOrigin)
        #print('Raster resolution:',self.xPixel,self.yPixel)

        # Read each raster band into 3D numpy array
        self.tiffData=ds.ReadAsArray(0,0,self.nX,self.nY)
        #print('Array shape:',self.tiffData.shape)
        #print('Array data:',self.tiffData[:,0,0])
        #print(self.tiffData[0,0])
        #print(self.tiffData[0,10979])

    def tiffExtract(self):
        '''Extract data from raster at footprint locations'''

        self.tiffExtract=np.zeros((self.use.shape[0],14))
        self.gediExtract=np.zeros((self.use.shape[0],9))
        print('tiffExtract shape',self.tiffExtract.shape)
        print('gediExtract shape',self.gediExtract.shape)

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

            tempID=self.useWaveID[i]
            #print('tempID',tempID)
            tempInd=np.where(self.shotNumb==tempID)
            #print('tempInd',tempInd)
            self.gediExtract[i,:]=self.gediData[tempInd,:]
            #print('Extracted gedi shot number',self.shotNumb[tempInd])
            #print('Extracted gedi data',self.gediData[tempInd,:])

        # Investigate those instances where L2A RH95 is v. large
        '''outliers=np.where((self.gediExtract[:,1] > 80))
        print('Numb of outliers',len(outliers[0]))
        print('Outlier ground elev.',self.gediExtract[outliers,0])
        print('Outlier ground residual',self.residual[outliers])
        print('Outlier metric ALS RH95',self.use[outliers,6])
        print('Outlier metric waveform RH95',self.use[outliers,5])'''


    def generateForest(self):
        '''Run Random Forest algorithm with the imported datasets'''

        # Apply a quality_flag = 1 mask to the training data
        mask=np.where((self.gediExtract[:,2]==1.0))
        print('mask length',len(mask[0]))

        # Independent variables are Sentinel reflectances and indices
        x=self.tiffExtract
        x_train=self.tiffExtract[mask[0],:]
        print('Shape of x',x_train.shape)
        # Dependent variable is GEDI RH95
        y_train=self.gediExtract[mask[0],1]
        print('Shape of y',y_train.shape)

        #x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.1,random_state=0)

        print('Running Random Forest algorithm ...')
        regressor=RandomForestRegressor(n_estimators=10)
        regressor.fit(x_train,y_train)
        self.y_pred=regressor.predict(x)
        print('y_pred shape',self.y_pred.shape)
        print('Random Forest coefficient',regressor.score(x_train,y_train))
        print('RF Feature Importance',regressor.feature_importances_)

    def utiliseForest(self):
        '''Use Random Forest output to correct real data'''

        # The ground residual as defined in alsComp.py
        self.residual=(self.gediExtract[:,0] - self.use[:,0])

        # The difference between the RF prediction of canopy height and the GEDI L2A height value
        self.RFresidual=(self.y_pred-self.gediExtract[:,1])

        # Identify waveforms where RF height prediction is more than 5m greater than L2A RH95
        badGround=np.where((self.RFresidual>5.0))
        print('Number of badGround waveforms',len(badGround[0]))

        # Histogram of ground residual values for these waveforms
        plt.hist(self.residual[badGround],bins=np.arange(-20,40+1,1),edgecolor='white',linewidth=2)
        plt.title('Ground Residual Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Ground Residual (m)')
        plt.savefig('../data/laselvaRF/badGroundHist.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

        # Get the preferred ground elevation for these waveforms and the minimum of the other six elevations
        badElev=self.gediExtract[badGround,0]
        print('badElev shape',badElev.shape)

        badMinElev=self.gediExtract[:,3:].min(axis=1)
        print('badMinElev shape',badMinElev.shape)
        badMinElev=badMinElev[badGround]
        print('badMinElev shape',badMinElev.shape)

        # How many unique elevations are there among the six algorithms for each waveform?
        '''elevs=self.gediExtract[:,3:]
        counts=np.zeros(elevs.shape[0],dtype=int)
        for i in range(elevs.shape[0]):
            unique=np.unique(elevs[i,:])
            counts[i]=(len(unique))
        uniqueElevs,uniqueCounts=np.unique(counts,return_counts=True)
        frequencies=np.asarray((uniqueElevs,uniqueCounts))
        print(frequencies)
        plt.hist(counts,bins=np.arange(0.5,5.5+1,1),edgecolor='white',linewidth=2)
        plt.title('Unique Ground Elevations')
        plt.ylabel('Frequency')
        plt.xlabel('Number of Unique Ground Elevations per Waveform')
        plt.savefig('../data/laselvaRF/grdElevDist.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()'''

    def compareForest(self):
        '''Make comparison plots of RF algorithm output'''

        # Histograms of height residual values
        plt.hist(self.RFresidual,bins=np.arange(-30.0,30.0+1,1),edgecolor='white',linewidth=2)
        plt.title('RF Height Residual Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('RF Height Residual (m)')
        plt.savefig('../data/laselvaRF/heightResHist.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Box plots of height residual as a function of canopy cover
        x1=self.RFresidual[np.where((self.use[:,2] < 0.1))]
        x2=self.RFresidual[np.where((self.use[:,2] > 0.1) & (self.use[:,2] < 0.2))]
        x3=self.RFresidual[np.where((self.use[:,2] > 0.2) & (self.use[:,2] < 0.3))]
        x4=self.RFresidual[np.where((self.use[:,2] > 0.3) & (self.use[:,2] < 0.4))]
        x5=self.RFresidual[np.where((self.use[:,2] > 0.4) & (self.use[:,2] < 0.5))]
        x6=self.RFresidual[np.where((self.use[:,2] > 0.5) & (self.use[:,2] < 0.6))]
        x7=self.RFresidual[np.where((self.use[:,2] > 0.6) & (self.use[:,2] < 0.7))]
        x8=self.RFresidual[np.where((self.use[:,2] > 0.7) & (self.use[:,2] < 0.8))]
        x9=self.RFresidual[np.where((self.use[:,2] > 0.8) & (self.use[:,2] < 0.9))]
        x10=self.RFresidual[np.where((self.use[:,2] > 0.9) & (self.use[:,2] < 1.0))]
        data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

        plt.boxplot(data,sym='')
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('Canopy Cover vs. RF Height Residual')
        plt.ylabel('RF Height Residual (m)')
        #plt.ylim([-30,40])
        plt.xlabel('Canopy Cover')
        plt.xticks(np.arange(1,11,step=1),['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
        plt.savefig('../data/laselvaRF/heightResCovBox.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Box plots of height residual as a function of GEDI L2A height
        x1=self.RFresidual[np.where((self.gediExtract[:,1] < 10))]
        x2=self.RFresidual[np.where((self.gediExtract[:,1] > 10) & (self.gediExtract[:,1] < 20))]
        x3=self.RFresidual[np.where((self.gediExtract[:,1] > 20) & (self.gediExtract[:,1] < 30))]
        x4=self.RFresidual[np.where((self.gediExtract[:,1] > 30) & (self.gediExtract[:,1] < 40))]
        x5=self.RFresidual[np.where((self.gediExtract[:,1] > 40) & (self.gediExtract[:,1] < 50))]
        x6=self.RFresidual[np.where((self.gediExtract[:,1] > 50) & (self.gediExtract[:,1] < 60))]
        data=[x1, x2, x3, x4, x5, x6]

        plt.boxplot(data,sym='')
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('Canopy Height vs. RF Height Residual')
        plt.ylabel('RF Height Residual (m)')
        #plt.ylim([-30,40])
        plt.xlabel('GEDI L2A RH95 (m)')
        plt.xticks(np.arange(1,7,step=1),['5', '15', '25', '35', '45', '55'])
        plt.savefig('../data/laselvaRF/heightResBox.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Random Forest predicted canopy height versus GEDI L2A canopy height
        plt.plot(self.gediExtract[:,1],self.y_pred,'o',markersize=1)
        plt.plot([0,80], [0,80],ls='--',color='r')
        plt.xlim([0,80])
        plt.ylim([0,80])
        plt.title('Random Forest Height Prediction')
        plt.xlabel('GEDI L2A RH95 (m)')
        plt.ylabel('RF Predicted Height (m)')
        plt.savefig('../data/laselvaRF/realVsPred.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Random Forest predicted canopy height versus ALS canopy height
        plt.plot(self.use[:,6],self.y_pred,'o',markersize=1)
        plt.plot([0,80], [0,80],ls='--',color='r')
        plt.xlim([0,80])
        plt.ylim([0,80])
        plt.title('Random Forest Height Prediction')
        plt.xlabel('ALS RH95 (m)')
        plt.ylabel('RF Predicted Height (m)')
        plt.savefig('../data/laselvaRF/ALSvsPred.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

        # Independently fit a linear regression to the RF output
        '''slope,intercept,r,p,se = sps.linregress(y_test,self.y_pred)
        print('Slope: {}, intercept: {}, r: {}, p: {}, std_err: {}'.format(slope,intercept,r,p,se))'''

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    comp=gridMetrics()
    comp.readMetric(args.inFile)
    comp.readGedi(args.inDir)
    comp.readTiff(args.inTiff)
    comp.tiffExtract()
    comp.generateForest()
    comp.utiliseForest()
    #comp.compareForest()
