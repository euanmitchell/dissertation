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
    p.add_argument('--inDir', dest='inDir', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/',
        help=('The path to the directory containing the GEDI L2A data.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/'))
    p.add_argument('--inDir2', dest='inDir2', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/gabon/l2a/',
        help=('The path to the directory containing the second GEDI L2A data.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/gabon/l2a/'))
    p.add_argument('--inFile', dest='inFile', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/metricAll.txt',
        help=('The path to a single simulated data file.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/metricAll.txt'))
    p.add_argument('--inTiff', dest='inTiff', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/laselva20mComposite.tif',
        help=('The path to a single tiff image.\nDefault is '))
    p.add_argument('--inTiff2', dest='inTiff2', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/gabon/gabon20mComposite.tif',
        help=('The path to the second tiff image.\nDefault is '))
    p.add_argument('--compare', dest='compare', action='store_true', default=False,
        help=('Run compareForest() method to make comparison plots.'))
    p.add_argument('--newGround', dest='newGround', action='store_true', default=False,
        help=('Run utiliseForest() method to generate new ground elevation estimates.'))
    p.add_argument('--beams', dest='beams', type=str, default='all',
        help=('Which beams to use: "all", "power", or "coverage".\nDefault is all'))
    cmdargs=p.parse_args()
    return cmdargs

#######################################

class RandomForest():
    '''
    A class to read output of gediMetric, GEDI L2A files, and imagery rasters
    and implement the Random Forest algorithm.
    '''

    def __init__(self):
        '''Class initialiser'''

    def readMetric(self,inFile):
        '''Read the metric data file'''

        print('Reading ALS metric file',inFile)

        def parse(id):
            return str(id).split('.')[2]

        self.waveID=np.loadtxt(inFile, usecols=(0,),converters={0: parse},dtype=np.uint64)

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

        self.use=np.delete(self.metricData,bad,axis=0)
        print('self.use shape',self.use.shape)

        self.useWaveID=np.delete(self.waveID,bad,axis=0)
        print('self.useWaveID shape',self.useWaveID.shape)

        # Define a diagonal to split the data into NW and SE datasets
        '''x=np.array((823600,830150))
        y=((x*1.389)+6500)
        plt.plot(self.use[:,7],self.use[:,8],'o',markersize=2)
        #plt.plot((823600,830150),(1150650,1159750))
        plt.plot(x,y)
        plt.show()
        nw=np.where((self.use[:,8])>((self.use[:,7]*1.389)+6500))
        se=np.where((self.use[:,8])<((self.use[:,7]*1.389)+6500))
        print('Len NW:',len(nw[0]))
        print('Len SE:',len(se[0]))'''

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

        # Read data from the second GEDI L2A directory
        '''print('Reading GEDI data from', inDir2)
        inputs='*.h5'
        path=os.path.join(inDir2,inputs)
        fileList=glob.glob(path)

        lowestMode=np.empty(0,dtype=float)
        rh95=np.empty(0,dtype=float)
        qualFlag=np.empty(0,dtype=float)
        ground1=np.empty(0,dtype=float)
        ground2=np.empty(0,dtype=float)
        ground3=np.empty(0,dtype=float)
        ground4=np.empty(0,dtype=float)
        ground5=np.empty(0,dtype=float)
        ground6=np.empty(0,dtype=float)
        lon84=np.empty(0,dtype='float64')
        lat84=np.empty(0,dtype='float64')

        for file in fileList:
            f=h5py.File(file,'r')
            for beam in beamlist:
                try:
                    lowestMode=np.append(lowestMode,f[beam]['elev_lowestmode'])
                    rh95=np.append(rh95,f[beam]['rh'][:,95])
                    qualFlag=np.append(qualFlag,f[beam]['quality_flag'])
                    ground1=np.append(ground1,f[beam]['geolocation']['elev_lowestmode_a1'])
                    ground2=np.append(ground2,f[beam]['geolocation']['elev_lowestmode_a2'])
                    ground3=np.append(ground3,f[beam]['geolocation']['elev_lowestmode_a3'])
                    ground4=np.append(ground4,f[beam]['geolocation']['elev_lowestmode_a4'])
                    ground5=np.append(ground5,f[beam]['geolocation']['elev_lowestmode_a5'])
                    ground6=np.append(ground6,f[beam]['geolocation']['elev_lowestmode_a6'])
                    lon84=np.append(lon84,f[beam]['lon_lowestmode'])
                    lat84=np.append(lat84,f[beam]['lat_lowestmode'])
                except:
                    continue
                    #print('Empty beam',beam)

        # Coordinates need converting out of EPSG 4326!
        inProj=Proj(init="epsg:4326")
        outProj=Proj(init="epsg:32732")
        lon,lat=transform(inProj,outProj,lon84,lat84)

        self.gediDataGabon=np.stack((lowestMode,rh95,qualFlag,ground1,ground2,ground3,ground4,ground5,ground6,lon,lat),axis=1)
        print('gediDataGabon shape',self.gediDataGabon.shape)

        self.useGabon=np.where((self.gediDataGabon[:,2]==1))
        print('useGabon shape',len(self.useGabon[0]))

        self.gabonGedi=self.gediDataGabon[self.useGabon[0],:]
        print('gabonData shape',self.gabonGedi.shape)'''

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

        # Read the second tiff into array
        '''ds=gdal.Open(filename2)
        proj=osr.SpatialReference(wkt=ds.GetProjection())
        self.epsgG=int(proj.GetAttrValue('Authority',1))
        self.nXG=ds.RasterXSize
        self.nYG=ds.RasterYSize
        transform=ds.GetGeoTransform()
        self.xOriginG=transform[0]
        self.yOriginG=transform[3]
        self.xPixelG=transform[1]
        self.yPixelG=transform[5]
        self.tiffDataGabon=ds.ReadAsArray(0,0,self.nXG,self.nYG)'''

    def tiffExtract(self):
        '''Extract data from raster at footprint locations'''

        self.tiffExtract=np.zeros((self.use.shape[0],14))
        self.gediExtract=np.zeros((self.use.shape[0],9))
        print('tiffExtract shape',self.tiffExtract.shape)
        print('gediExtract shape',self.gediExtract.shape)

        for i in range(self.use.shape[0]):  #self.use.shape[0]
            tempLon=self.use[i,7]
            tempLat=self.use[i,8]

            xDist=tempLon-self.xOrigin
            yDist=tempLat-self.yOrigin

            xInd=int((tempLon-self.xOrigin)//self.xPixel)
            yInd=int((tempLat-self.yOrigin)//self.yPixel)

            self.tiffExtract[i,:]=self.tiffData[:,yInd,xInd]

            tempID=self.useWaveID[i]
            tempInd=np.where(self.shotNumb==tempID)
            self.gediExtract[i,:]=self.gediData[tempInd,:]

        # Extract data from the second tiff
        '''self.tiffExtractGabon=np.zeros((self.gabonGedi.shape[0],14))
        print('tiffExtract shape',self.tiffExtractGabon.shape)

        for i in range(self.tiffExtractGabon.shape[0]):
            tempLon=self.gabonGedi[i,9]
            tempLat=self.gabonGedi[i,10]

            xDist=tempLon-self.xOriginG
            yDist=tempLat-self.yOriginG

            xInd=int((tempLon-self.xOriginG)//self.xPixelG)
            yInd=int((tempLat-self.yOriginG)//self.yPixelG)

            self.tiffExtractGabon[i,:]=self.tiffDataGabon[:,yInd,xInd]

        print('Gabon tiff shape',self.tiffExtractGabon.shape)'''

    def generateForest(self):
        '''Run Random Forest algorithm with the imported datasets'''

        plt.rcParams['figure.figsize']=(16,12)
        plt.rcParams['xtick.labelsize']=16
        plt.rcParams['ytick.labelsize']=16
        plt.rcParams['axes.labelsize']=18
        plt.rcParams['axes.labelpad']=8.0
        plt.rcParams['axes.titlesize']=20

        # Apply a quality_flag = 1 mask to the training data and define test data as the other data
        self.train=np.where((self.gediExtract[:,2]==1.0))
        #self.test=np.where((self.gediExtract[:,2]==0.0))

        # Split into train and test datasets based on location instead
        #self.train=np.where(((self.use[:,8])<((self.use[:,7]*1.389)+6500)))
        #self.test=np.where((self.use[:,8])>((self.use[:,7]*1.389)+6500))

        #print('train length',len(self.train[0]))
        #print('test length',len(self.test[0]))

        # Make histograms of the ALS and Sentinel properties for the two datasets
        '''plt.hist(self.tiffExtract[self.train[0],2],bins=np.arange(0.0,1000+20,20))
        plt.title('Sentinel Red Distribution Train')
        plt.ylabel('Frequency')
        plt.xlabel('Intensity')
        #plt.savefig(args.outRoot+'coverHist.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()'''

        x=self.tiffExtract[self.train[0],:]
        y=self.gediExtract[self.train[0],1]
        x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.5,random_state=0)

        # Independent variables are Sentinel reflectances and indices
        #x_train=self.tiffExtract[self.train[0],:]
        #x_test=self.tiffExtract[self.test[0],:]
        print('Shape of x',x_train.shape)
        # Dependent variable is GEDI RH95
        #y_train=self.gediExtract[self.train[0],1]
        #y_test=self.gediExtract[self.test[0],1]
        print('Shape of y',y_train.shape)

        #for i in range(1,21):
        print('Running Random Forest algorithm ...')
        regressor=RandomForestRegressor(n_estimators=200,random_state=0)
        regressor.fit(x_train,y_train)
        self.y_pred=regressor.predict(x_test)
        print('y_pred shape',self.y_pred.shape)
        print('Random Forest coefficient',regressor.score(x_test,y_test))
        print('Random Forest MSE',mean_squared_error(y_test,self.y_pred))
        print('Random Forest RMSE',sqrt(mean_squared_error(y_test,self.y_pred)))
        print('Random Forest ALS RMSE',sqrt(mean_squared_error(self.use[self.test[0],6],self.y_pred)))
        print('Random Forest Bias',np.sum(self.y_pred-y_test)/self.y_pred.shape[0])
        print('Random Forest ALS Bias',np.sum(self.y_pred-self.use[self.test[0],6])/self.y_pred.shape[0])
        #print('RF Feature Importance',regressor.feature_importances_)

        # Make canopy height map of RF model output
        '''plt.scatter(self.use[self.test[0],7],self.use[self.test[0],8],s=8,c=self.y_pred,cmap='Greens')
        plt.xlabel('Easting (m)')
        plt.ylabel('Northing (m)')
        plt.legend()
        plt.show()'''

        # Independently fit a linear regression to the RF output
        slope,intercept,r,p,se = sps.linregress(y_test,self.y_pred)
        print('Slope: {}, intercept: {}, r: {}, p: {}, std_err: {}'.format(slope,intercept,r,p,se))

        # Random Forest predicted canopy height versus GEDI L2A canopy height
        plt.plot(y_test,self.y_pred,'o',markersize=1)
        plt.plot([0,80], [0,80],ls='--',color='r')
        plt.xlim([0,80])
        plt.ylim([0,80])
        #plt.title('Random Forest Height Prediction')
        plt.xlabel('GEDI L2A RH95 (m)')
        plt.ylabel('RF Predicted Height (m)')
        #plt.savefig('../data/laselvaRF/realVsPred.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

        # Plot histogram of real and predicted canopy height
        '''plt.hist([self.gabonGedi[:,1],self.y_pred],bins=np.arange(0.0,80+5,5),edgecolor='white',linewidth=2,label=['GEDI L2A','Predicted'])
        #plt.title('Canopy Height Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Canopy Height (m)')
        plt.legend()
        #plt.savefig('../data/gabonFigs/heightHist.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()'''

        # The ground residual as defined in alsComp.py
        self.residual=(self.gediExtract[self.test[0],0] - self.use[self.test[0],0])

        # The difference between the RF prediction of canopy height and the GEDI L2A height value
        self.RFresidual=(self.y_pred-y_test)

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
        #plt.title('Canopy Height vs. RF Height Residual')
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
        #plt.title('Random Forest Height Prediction')
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

    def utiliseForest(self):
        '''Use Random Forest output to correct real data'''

        # Identify waveforms where RF height prediction is more than 5m greater than L2A RH95
        '''badGround=np.where((self.RFresidual>5.0))
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
        elevs=self.gediExtract[:,3:]
        counts=np.zeros(elevs.shape[0],dtype=int)
        for i in range(elevs.shape[0]):
            unique=np.unique(elevs[i,:])
            counts[i]=(len(unique))
        counts=counts[badGround]
        uniqueElevs,uniqueCounts=np.unique(counts,return_counts=True)
        frequencies=np.asarray((uniqueElevs,uniqueCounts))
        print(frequencies)
        plt.hist(counts,bins=np.arange(0.5,5.5+1,1),edgecolor='white',linewidth=2)
        plt.title('Unique Bad Ground Elevations')
        plt.ylabel('Frequency')
        plt.xlabel('Number of Unique Ground Elevations per Bad Waveform')
        plt.savefig('../data/laselvaRF/grdElevDistBad.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()'''

        # Plot all the different ground elevation estimates from the L2A data
        '''badPref=self.gediExtract[badGround,0]
        badAlg1=self.gediExtract[badGround,3]
        badAlg2=self.gediExtract[badGround,4]
        badAlg3=self.gediExtract[badGround,5]
        badAlg4=self.gediExtract[badGround,6]
        badAlg5=self.gediExtract[badGround,7]
        badAlg6=self.gediExtract[badGround,8]
        plt.plot(np.arange(badElev.shape[1]),badAlg1[0],'o',markersize=3)
        plt.plot(np.arange(badElev.shape[1]),badAlg2[0],'o',markersize=3)
        plt.plot(np.arange(badElev.shape[1]),badAlg3[0],'o',markersize=3)
        plt.plot(np.arange(badElev.shape[1]),badAlg4[0],'o',markersize=3)
        plt.plot(np.arange(badElev.shape[1]),badAlg5[0],'o',markersize=3)
        plt.plot(np.arange(badElev.shape[1]),badAlg6[0],'o',markersize=3)
        plt.plot(np.arange(badElev.shape[1]),badPref[0],'o',color='red',markersize=3)
        plt.xlabel('Waveform Index')
        plt.ylabel('Ground Elevation (m)')
        plt.savefig('../data/laselvaRF/elevations.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()'''

        # Build a new array of ground elevations
        newGround=np.zeros(self.use[self.test[0]].shape[0])
        missedThreshold=0
        madeThreshold=0
        updated=0
        for i in range(self.use[self.test[0]].shape[0]):
            # If predicted height within threshold of gedi value, use the L2A preferred algorithm elevation
            if self.RFresidual[i]<5.0:
                newGround[i]=self.gediExtract[i,0]
                missedThreshold+=1

            # Original implementation - get lowest value (if 1 or 2) or the mean (if 3 or more)
            else:
                madeThreshold+=1
                # Get the unique ground values from the six L2A algorithms - this will be sorted!
                options=np.unique(self.gediExtract[i,3:])
                # If there are only 1 or 2 options, take the lowest
                if options.shape[0]<=2:
                    newGround[i]=options[0]
                # If there are 3 or more options take a mean after filtering anything higher than the preferred setting
                else:
                    if options[0]==self.gediExtract[i,0]:
                        newGround[i]=self.gediExtract[i,0]
                    else:
                        lowOptions=np.where((options<self.gediExtract[i,0]))
                        newGround[i]=np.mean(options[lowOptions])
                if newGround[i] != self.gediExtract[i,0]:
                    updated+=1

            '''# Updated implementation - calculating a 'predicted ground' elevation
            else:
                madeThreshold+=1
                # Get the unique ground values from the six L2A algorithms - this will be sorted!
                options=np.unique(self.gediExtract[i,3:])
                # Get the RFresidual for this waveform
                residual=self.RFresidual[i]
                # Subtract this from preferred ground elevation
                inferredGround=self.gediExtract[i,0]-residual
                # Calculate difference btw inferred ground and six options
                diffs=options-inferredGround
                # Get index of minimum diff
                minInd=np.argmin(np.absolute(diffs))
                updatedGround=options[minInd]
                newGround[i]=updatedGround
                if updatedGround != self.gediExtract[i,0]:
                    updated+=1'''

        print('Count that missed the threshold',missedThreshold)
        print('Count that made the threshold',madeThreshold)
        print('Count that were actually changed',updated)

        #print('New ground elevation array',newGround)

        newResidual=(newGround-self.use[self.test[0],0])

        mask=np.where((self.gediExtract[self.test[0],2]==1.0))

        manualRMSE=sqrt(np.sum(self.residual**2)/self.use.shape[0])
        manualBias=(np.sum(self.residual))/self.use.shape[0]
        manualQualBias=(np.sum(self.residual[mask]))/self.use[mask].shape[0]
        newGroundBias=(np.sum(newResidual))/self.use.shape[0]

        print('Original Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,0])))
        print('Original Ground Bias',manualBias)
        print('A1 Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,3])))
        print('A1 Ground Bias',(np.sum(self.gediExtract[:,3]-self.use[:,0]))/self.use.shape[0])
        print('A2 Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,4])))
        print('A2 Ground Bias',(np.sum(self.gediExtract[:,4]-self.use[:,0]))/self.use.shape[0])
        print('A3 Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,5])))
        print('A3 Ground Bias',(np.sum(self.gediExtract[:,5]-self.use[:,0]))/self.use.shape[0])
        print('A4 Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,6])))
        print('A4 Ground Bias',(np.sum(self.gediExtract[:,6]-self.use[:,0]))/self.use.shape[0])
        print('A5 Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,7])))
        print('A5 Ground Bias',(np.sum(self.gediExtract[:,7]-self.use[:,0]))/self.use.shape[0])
        print('A6 Ground RMSE',sqrt(mean_squared_error(self.use[:,0],self.gediExtract[:,8])))
        print('A6 Ground Bias',(np.sum(self.gediExtract[:,8]-self.use[:,0]))/self.use.shape[0])
        print('Original Ground Qual RMSE',sqrt(mean_squared_error(self.use[mask,0],self.gediExtract[mask,0])))
        print('Original Ground Qual Bias',manualQualBias)
        print('New Ground RMSE',sqrt(mean_squared_error(self.use[self.test[0],0],newGround)))
        print('New Ground Bias',newGroundBias)


        # Scatter plot of ALS cover versus residuals
        plt.plot(self.use[self.test[0],2],self.residual,'o',label='Original Data')
        plt.plot(self.use[self.test[0],2],newResidual,'o',markersize=3,label='Corrected Data')
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.legend()
        #plt.ylim([-30,40])
        plt.xlabel('ALS Cover')
        plt.ylabel('Ground Residual (m)')
        plt.show()

        # Box plots of new residual vs. ALS cover
        x1=newResidual[np.where((self.use[self.test[0],2] < 0.1))]
        x2=newResidual[np.where((self.use[self.test[0],2] > 0.1) & (self.use[self.test[0],2] < 0.2))]
        x3=newResidual[np.where((self.use[self.test[0],2] > 0.2) & (self.use[self.test[0],2] < 0.3))]
        x4=newResidual[np.where((self.use[self.test[0],2] > 0.3) & (self.use[self.test[0],2] < 0.4))]
        x5=newResidual[np.where((self.use[self.test[0],2] > 0.4) & (self.use[self.test[0],2] < 0.5))]
        x6=newResidual[np.where((self.use[self.test[0],2] > 0.5) & (self.use[self.test[0],2] < 0.6))]
        x7=newResidual[np.where((self.use[self.test[0],2] > 0.6) & (self.use[self.test[0],2] < 0.7))]
        x8=newResidual[np.where((self.use[self.test[0],2] > 0.7) & (self.use[self.test[0],2] < 0.8))]
        x9=newResidual[np.where((self.use[self.test[0],2] > 0.8) & (self.use[self.test[0],2] < 0.9))]
        x10=newResidual[np.where((self.use[self.test[0],2] > 0.9) & (self.use[self.test[0],2] < 1.0))]
        data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

        plt.boxplot(data,sym='')
        plt.axhline(y=0.0,color='r',linestyle='--')
        #plt.title('Canopy Cover vs. New Ground Residual')
        plt.ylabel('Ground Residual (m)')
        #plt.ylim([-30,40])
        plt.xlabel('Canopy Cover')
        plt.xticks(np.arange(1,11,step=1),['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
        #plt.savefig('../data/laselvaRF/heightResCovBox.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    comp=RandomForest()
    comp.readMetric(args.inFile)
    comp.readGedi(args.inDir)
    comp.readTiff(args.inTiff)
    comp.tiffExtract()
    comp.generateForest()
    if args.compare:
        comp.compareForest()
    if args.newGround:
        comp.utiliseForest()
