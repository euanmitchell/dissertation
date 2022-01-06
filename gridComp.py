#######################################
# Script to plot output from running
# gediMetric over the simulated
# waveforms from gediRat.
#######################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
from osgeo import gdal
from pyproj import Proj, transform
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
    p.add_argument('--inTiff', dest='inTiff', type=str, default='../data/laselva/laselvaNDVI_20m.tif',
        help=('The path to a single tiff image.\nDefault is '))
    p.add_argument('--tiff', dest='tiff', action='store_true', default=False,
        help=('Read the specified tiff and make relevant plots.'))
    p.add_argument('--simPlots', dest='simPlots', action='store_true', default=False,
        help=('Make histogram and box plots from the simulated data.'))
    p.add_argument('--outRoot', dest='outRoot', type=str, default='../data/figures/',
        help=('Output root for plots.\nDefault is "../data/figures/"'))
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

        # Read raster image into numpy array
        self.tiffData=ds.GetRasterBand(1).ReadAsArray(0,0,self.nX,self.nY)
        self.tiffDataZeros=np.where(self.tiffData==101, 0, self.tiffData)
        print('Array shape:',self.tiffData.shape)
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
        print(bad)

        self.use=np.delete(self.data,bad,axis=0)
        print(self.use.shape)

        self.residual = (self.use[:,3] - self.use[:,0])
        sorted = np.sort(self.residual)
        #print(sorted[:25])
        #print(sorted[-25:])

        #print(residual.shape)
        #print(np.amin(residual))
        #print(np.amax(residual))

        inProj=Proj(init="epsg:32616")
        outProj=Proj(init="epsg:4326")
        self.lons,self.lats=transform(inProj,outProj,self.use[:,7],self.use[:,8])
        print('Re-projected lons array shape',self.lons.shape)

        if args.simPlots:
            plt.hist(self.use[:,2],bins=np.arange(0.0,1.0+0.02,0.02),edgecolor='white',linewidth=2)
            plt.title('Canopy Cover Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('ALS Canopy Cover')
            plt.savefig(args.outRoot + 'cov98wid05.coverHist.png')
            plt.close()
            plt.clf()
            #plt.show()

            plt.hist(self.residual,bins=np.arange(-20.0,20.0+1,1),edgecolor='white',linewidth=2)
            plt.title('Ground Residual Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('Ground Residual (m)')
            plt.savefig(args.outRoot + 'cov98wid05.resHist.png')
            plt.close()
            plt.clf()
            #plt.show()

            x1=self.residual[np.where((self.use[:,2] < 0.1))]
            x2=self.residual[np.where((self.use[:,2] > 0.1) & (self.use[:,2] < 0.2))]
            x3=self.residual[np.where((self.use[:,2] > 0.2) & (self.use[:,2] < 0.3))]
            x4=self.residual[np.where((self.use[:,2] > 0.3) & (self.use[:,2] < 0.4))]
            x5=self.residual[np.where((self.use[:,2] > 0.4) & (self.use[:,2] < 0.5))]
            x6=self.residual[np.where((self.use[:,2] > 0.5) & (self.use[:,2] < 0.6))]
            x7=self.residual[np.where((self.use[:,2] > 0.6) & (self.use[:,2] < 0.7))]
            x8=self.residual[np.where((self.use[:,2] > 0.7) & (self.use[:,2] < 0.8))]
            x9=self.residual[np.where((self.use[:,2] > 0.8) & (self.use[:,2] < 0.9))]
            x10=self.residual[np.where((self.use[:,2] > 0.9) & (self.use[:,2] < 1.0))]
            values=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

            plt.boxplot(values,sym='')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title('Canopy Cover vs. Residual')
            plt.ylabel('Ground Residual (m)')
            plt.ylim([-15,15])
            plt.xlabel('ALS Canopy Cover')
            plt.xticks(np.arange(1,11,step=1),['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
            plt.savefig(args.outRoot + 'cov98wid05.box.png')
            plt.close()
            plt.clf()
            #plt.show()

            plt.plot(self.use[:,2],self.residual,'o')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title('ALS Cover vs. Ground Residual')
            plt.xlabel('ALS Cover')
            plt.ylabel('Ground Residual (m)')
            #plt.xlim([200,1300])
            plt.ylim([-30,30])
            plt.savefig(args.outRoot + 'cov98wid05.scat.png')
            plt.close()
            plt.clf()
            #plt.show()

            '''nbins=1000
            x=self.use[:,2]
            y=self.residual
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.title('ALS Cover vs Ground Residual')
            plt.xlabel('ALS ALS Cover')
            plt.ylabel('Ground Residual (m)')
            plt.colorbar()
            plt.savefig(args.outRoot + 'cov99wid05.dens.png')
            plt.close()
            plt.clf()
            #plt.show()'''

    def tiffPlot(self):
        '''Plot data from metric file against imagery'''

        self.tiffExtract=np.zeros(self.residual.shape[0])
        print(self.tiffExtract.shape[0])

        for i in range(self.use.shape[0]):  #self.use.shape[0]
            tempLon=self.use[i,7]
            #tempLon=self.lons[i]
            #tempLon=-84.011111
            tempLat=self.use[i,8]
            #tempLat=self.lats[i]
            #tempLat=10.430643
            #print('First footprint lon',tempLon)
            #print('First footprint lat',tempLat)

            xDist=tempLon-self.xOrigin
            yDist=tempLat-self.yOrigin
            #print(xDist,yDist)

            xInd=int((tempLon-self.xOrigin)//self.xPixel)
            yInd=int((tempLat-self.yOrigin)//self.yPixel)
            #print(xInd,yInd)

            #self.tiffExtract[i]=np.mean(self.tiffDataZeros[yInd-1:yInd+2,xInd-1:xInd+2])
            self.tiffExtract[i]=self.tiffDataZeros[yInd,xInd]
            #if ((self.use[i,2] > 0.01) & (self.use[i,2] < 0.15)):
            #print('Tiff value',self.tiffData[yInd,xInd])
            #print('Footprint Long.',tempLon)
            #print('Footprint Lat.',tempLat)
            #print('ALS Cover',self.use[i,2])
            #print('Footprint Cover',self.use[i,4])

            #print('First footprint NDVI values',self.tiffData[yInd-1:yInd+1,xInd-1:xInd+1])
            #print('First footprint cover value',self.use[i,2])

        plt.plot(self.use[:,6],self.tiffExtract,'o',markersize=1)
        plt.title('ALS Height vs. Footprint NDVI')
        plt.xlabel('ALS RH95 (m)')
        plt.ylabel('Footprint NDVI')
        #plt.savefig(args.outRoot + 'ndwi_11CoverMean5x5.png')
        #plt.close()
        #plt.clf()
        plt.show()

        '''x1=self.tiffExtract[np.where((self.use[:,2] < 0.1))]
        x2=self.tiffExtract[np.where((self.use[:,2] > 0.1) & (self.use[:,2] < 0.2))]
        x3=self.tiffExtract[np.where((self.use[:,2] > 0.2) & (self.use[:,2] < 0.3))]
        x4=self.tiffExtract[np.where((self.use[:,2] > 0.3) & (self.use[:,2] < 0.4))]
        x5=self.tiffExtract[np.where((self.use[:,2] > 0.4) & (self.use[:,2] < 0.5))]
        x6=self.tiffExtract[np.where((self.use[:,2] > 0.5) & (self.use[:,2] < 0.6))]
        x7=self.tiffExtract[np.where((self.use[:,2] > 0.6) & (self.use[:,2] < 0.7))]
        x8=self.tiffExtract[np.where((self.use[:,2] > 0.7) & (self.use[:,2] < 0.8))]
        x9=self.tiffExtract[np.where((self.use[:,2] > 0.8) & (self.use[:,2] < 0.9))]
        x10=self.tiffExtract[np.where((self.use[:,2] > 0.9) & (self.use[:,2] < 1.0))]
        values=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]
        plt.boxplot(values,sym='')
        #plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('ALS Cover vs. Footprint NDWI')
        plt.xlabel('ALS Cover')
        plt.ylabel('Footprint Band 11 NDWI')
        plt.ylim([-1,1])
        plt.xticks(np.arange(1,11,step=1),['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
        plt.savefig(args.outRoot + 'ndwi_11CoverBoxMean5x5.png')
        plt.close()
        plt.clf()
        #plt.show()

        plt.plot(self.use[:,6],self.tiffExtract,'o',markersize=1)
        plt.plot([-2,60], [-2,60],ls='--',color='r')
        plt.title('ALS Height vs. Potapov Height')
        plt.xlabel('ALS RH95 (m)')
        plt.ylabel('Potapov Forest Height (m)')
        plt.xlim([-2,60])
        plt.ylim([-2,60])
        plt.savefig(args.outRoot + 'PotapovComp1z.png')
        plt.close()
        plt.clf()
        #plt.show()

        x1=self.tiffExtract[np.where((self.use[:,6] < 5))]
        x2=self.tiffExtract[np.where((self.use[:,6] > 5) & (self.use[:,6] < 10))]
        x3=self.tiffExtract[np.where((self.use[:,6] > 10) & (self.use[:,6] < 15))]
        x4=self.tiffExtract[np.where((self.use[:,6] > 15) & (self.use[:,6] < 20))]
        x5=self.tiffExtract[np.where((self.use[:,6] > 20) & (self.use[:,6] < 25))]
        x6=self.tiffExtract[np.where((self.use[:,6] > 25) & (self.use[:,6] < 30))]
        x7=self.tiffExtract[np.where((self.use[:,6] > 30) & (self.use[:,6] < 35))]
        x8=self.tiffExtract[np.where((self.use[:,6] > 35) & (self.use[:,6] < 40))]
        x9=self.tiffExtract[np.where((self.use[:,6] > 40) & (self.use[:,6] < 45))]
        x10=self.tiffExtract[np.where((self.use[:,6] > 45) & (self.use[:,6] < 50))]
        x11=self.tiffExtract[np.where((self.use[:,6] > 50) & (self.use[:,6] < 55))]
        x12=self.tiffExtract[np.where((self.use[:,6] > 55) & (self.use[:,6] < 60))]
        values=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12]
        plt.boxplot(values,sym='')
        #plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('ALS Height vs. Footprint NDWI')
        plt.ylabel('Footprint Band 11 NDWI')
        plt.ylim([-1,1])
        plt.xlabel('ALS RH95 (m)')
        plt.xticks(np.arange(1,13,step=1),['2.5', '7.5', '12.5', '17.5', '22.5', '27.5', '32.5', '37.5', '42.5', '47.5', '52.5', '57.5'])
        plt.savefig(args.outRoot + 'ndwi_11HeightBoxMean5x5.png')
        plt.close()
        plt.clf()
        #plt.show()'''

        '''nbins=1000
        x=self.use[:,1]
        y=self.tiffExtract
        k=kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(),yi.flatten()]))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        plt.title('ALS Slope vs. Footprint Slope')
        plt.xlabel('ALS Slope')
        plt.ylabel('Footprint SRTM Slope')
        plt.colorbar()
        plt.savefig(args.outRoot + 'slopeCoverDens.png')
        plt.close()
        plt.clf()
        #plt.show()'''

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    comp=gridMetrics()
    comp.readMetric(args.inFile)
    if args.tiff:
        comp.readTiff(args.inTiff)
        comp.tiffPlot()
