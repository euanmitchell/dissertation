#######################################
# Script to plot Version 2.0 preferred
# ground estimate against collocated
# ALS ground estimate.
#######################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
from sklearn.metrics import mean_squared_error
from math import sqrt
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
    p.add_argument('--als', dest='als', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/',
        help=('The path to the directory containing the ALS metric files.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/'))
    p.add_argument('--alsFile', dest='alsFile', type=str, default='',
        help=('The path to a single ALS metric file.\nDefault is not set'))
    p.add_argument('--gedi', dest='gedi', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/',
        help=('The path to the directory containing the GEDI data.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/'))
    p.add_argument('--gediFile', dest='gediFile', type=str, default='',
        help=('The path to a single GEDI data file.\nDefault is not set'))
    p.add_argument('--beams', dest='beams', type=str, default='all',
        help=('Which beams to use: "all", "power", or "coverage".\nDefault is all'))
    p.add_argument('--outRoot', dest='outRoot', type=str, default='../data/figures/',
        help=('Output root for plots.\nDefault is "../data/figures/"'))
    p.add_argument('--site', dest='site', type=str, default='',
        help=('Site name for plot titles.\nDefault is not set'))
    p.add_argument('--plots', dest='plots', action='store_true', default=False,
        help=('Call the plotData() method to make comparision scatter plots'))
    p.add_argument('--hist', dest='hist', action='store_true', default=False,
        help=('Make histogram plots'))
    p.add_argument('--box', dest='box', action='store_true', default=False,
        help=('Make box plots'))
    p.add_argument('--dens', dest='dens', action='store_true', default=False,
        help=('Make density plots'))
    p.add_argument('--waveforms', dest='waveforms', action='store_true', default=False,
        help=('Generate shell script to copy waveforms meeting criteria specified in file'))
    cmdargs=p.parse_args()
    return cmdargs

#######################################

class compareGround(object):
    '''A class to handle reading and plotting ALS and GEDI ground estimates'''

    def __init__(self):
        '''Initialise the class'''

    def readALS(self,input):
        '''Read the output text file from gediMetric'''

        '''SHOULD BE ABLE TO CAT ALL THE METRIC TEXT FILES INTO ONE AND
           READ IN WITH numpy.loadfromtxt. THAT WILL GENERATE A SINGLE
           2D ARRAY INSTEAD OF ALL THESE SEPARATE ARRAYS.'''

        print('Reading ALS metric file',input)
        # Create empty array for shotnumbers
        self.alsShot=np.empty(0)
        # Create empty arrays for ALS data
        self.alsGround=np.empty(0,dtype=float)      # 'true ground'
        self.slope=np.empty(0,dtype=float)          # 'ground slope'
        self.alsCover=np.empty(0,dtype=float)       # 'ALS cover'
        self.alsRH95=np.empty(0,dtype=float)        # 'rhReal 95'
        self.overlap=np.empty(0,dtype=float)        # 'groundOverlap'

        with open(input, 'r') as f:
            for line in f:
                if line[0]=='#':
                    continue
                else:
                    data=line.split(' ')
                    ident=data[0].split('.')
                    shot=ident[2]
                    self.alsShot=np.append(self.alsShot,shot)
                    self.alsGround=np.append(self.alsGround,data[1])
                    self.slope=np.append(self.slope,data[3])
                    self.alsCover=np.append(self.alsCover,data[4])
                    self.alsRH95=np.append(self.alsRH95,data[95])
                    self.overlap=np.append(self.overlap,data[108])


    def readGEDI(self,input):
        '''Read the GEDI L2A file'''

        print('Reading GEDI L2A file',input)
        if args.beams == 'all':
            beamlist=['BEAM0101','BEAM1000','BEAM0010','BEAM0000','BEAM0011','BEAM0110','BEAM0001','BEAM1011']
        elif args.beams == 'power':
            beamlist=['BEAM0101','BEAM1000','BEAM0110','BEAM1011']
        elif args.beams == 'coverage':
            beamlist=['BEAM0010','BEAM0000','BEAM0011','BEAM0001']

        # Create empty array for shotnumbers
        self.gediShot=np.empty(0,dtype=str)
        # Create empty arrays for GEDI ground estimate and beam sensitivity
        self.gediGround=np.empty(0,dtype=float)
        self.gediSens=np.empty(0,dtype=float)
        self.gediSolar=np.empty(0,dtype=float)
        self.gediQual=np.empty(0,dtype=int)
        self.algorithm=np.empty(0,dtype=int)
        self.ground1=np.empty(0,dtype=float)
        self.ground2=np.empty(0,dtype=float)
        self.ground3=np.empty(0,dtype=float)
        self.ground4=np.empty(0,dtype=float)
        self.ground5=np.empty(0,dtype=float)
        self.ground6=np.empty(0,dtype=float)

        f=h5py.File(input,'r')
        for beam in beamlist:
            # Need to handle empty beams
            try:
                self.gediShot=np.append(self.gediShot,f[beam]['shot_number'])
                self.gediGround=np.append(self.gediGround,f[beam]['elev_lowestmode'])
                self.gediSens=np.append(self.gediSens,f[beam]['sensitivity'])
                self.gediSolar=np.append(self.gediSolar,f[beam]['solar_elevation'])
                self.gediQual=np.append(self.gediQual,f[beam]['quality_flag'])
                self.algorithm=np.append(self.algorithm,f[beam]['selected_algorithm'])
                self.ground1=np.append(self.ground1,f[beam]['geolocation']['elev_lowestmode_a1'])
                self.ground2=np.append(self.ground2,f[beam]['geolocation']['elev_lowestmode_a2'])
                self.ground3=np.append(self.ground3,f[beam]['geolocation']['elev_lowestmode_a3'])
                self.ground4=np.append(self.ground4,f[beam]['geolocation']['elev_lowestmode_a4'])
                self.ground5=np.append(self.ground5,f[beam]['geolocation']['elev_lowestmode_a5'])
                self.ground6=np.append(self.ground6,f[beam]['geolocation']['elev_lowestmode_a6'])
                print('Data in',beam)
            except:
                print('Empty beam',beam)

    def plotData(self):
        '''Identify and plot corresponding data'''

        self.alsShotClean=np.empty(0)
        self.alsGroundClean=np.empty(0,dtype=float)
        self.slopeClean=np.empty(0,dtype=float)
        self.alsCoverClean=np.empty(0,dtype=float)
        self.alsRH95Clean=np.empty(0,dtype=float)
        self.overlapClean=np.empty(0,dtype=float)
        self.gediShotClean=np.empty(0)
        self.gediGroundClean=np.empty(0,dtype=float)
        self.gediSensClean=np.empty(0,dtype=float)
        self.gediSolarClean=np.empty(0,dtype=float)
        self.gediQualClean=np.empty(0,dtype=int)
        self.algorithmClean=np.empty(0,dtype=int)
        self.ground1Clean=np.empty(0,dtype=float)
        self.ground2Clean=np.empty(0,dtype=float)
        self.ground3Clean=np.empty(0,dtype=float)
        self.ground4Clean=np.empty(0,dtype=float)
        self.ground5Clean=np.empty(0,dtype=float)
        self.ground6Clean=np.empty(0,dtype=float)

        # Match each record in the ALS metric file to the corresponding L2A GEDI record
        print('Sorting and matching arrays ....')
        for i in range(alsShotMaster.shape[0]):
            for j in range(gediShotMaster.shape[0]):
                if alsShotMaster[i] == gediShotMaster[j]:
                    self.alsShotClean=np.append(self.alsShotClean,alsShotMaster[i])
                    self.alsGroundClean=np.append(self.alsGroundClean,alsGroundMaster[i])
                    self.slopeClean=np.append(self.slopeClean,slopeMaster[i])
                    self.alsCoverClean=np.append(self.alsCoverClean,alsCoverMaster[i])
                    self.alsRH95Clean=np.append(self.alsRH95Clean,alsRH95Master[i])
                    self.overlapClean=np.append(self.overlapClean,overlapMaster[i])
                    self.gediShotClean=np.append(self.gediShotClean,gediShotMaster[j])
                    self.gediGroundClean=np.append(self.gediGroundClean,gediGroundMaster[j])
                    self.gediSensClean=np.append(self.gediSensClean,gediSensMaster[j])
                    self.gediSolarClean=np.append(self.gediSolarClean,gediSolarMaster[j])
                    self.gediQualClean=np.append(self.gediQualClean,gediQualMaster[j])
                    self.algorithmClean=np.append(self.algorithmClean,algorithmMaster[j])
                    self.ground1Clean=np.append(self.ground1Clean,ground1Master[j])
                    self.ground2Clean=np.append(self.ground2Clean,ground2Master[j])
                    self.ground3Clean=np.append(self.ground3Clean,ground3Master[j])
                    self.ground4Clean=np.append(self.ground4Clean,ground4Master[j])
                    self.ground5Clean=np.append(self.ground5Clean,ground5Master[j])
                    self.ground6Clean=np.append(self.ground6Clean,ground6Master[j])

        print('self.alsShotClean length',self.alsShotClean.shape[0])
        print('self.gediShotClean length',self.gediShotClean.shape[0])

        good=np.where((self.gediQualClean == 1))
        bad=np.where((self.gediQualClean == 0))
        print('Number of good waveforms',len(good[0]))
        print('Number of bad waveforms',len(bad[0]))

        #cov90=np.where((self.alsCoverClean.astype(float) > 0.9))
        #cov95=np.where((self.alsCoverClean.astype(float) > 0.95))
        #print('Cov90 waveforms',len(cov90[0]))
        #print('Cov95 waveforms',len(cov95[0]))

        residual=np.subtract(self.gediGroundClean,self.alsGroundClean.astype(float))
        residual1=np.subtract(self.ground1Clean,self.alsGroundClean.astype(float))
        residual2=np.subtract(self.ground2Clean,self.alsGroundClean.astype(float))
        residual3=np.subtract(self.ground3Clean,self.alsGroundClean.astype(float))
        residual4=np.subtract(self.ground4Clean,self.alsGroundClean.astype(float))
        residual5=np.subtract(self.ground5Clean,self.alsGroundClean.astype(float))
        residual6=np.subtract(self.ground6Clean,self.alsGroundClean.astype(float))

        # Specify a subset of data to plot here with useInd?
        #useInd=np.arange(self.gediGroundClean.shape[0])
        #useInd=np.where((self.alsGroundClean.astype(float) > 0) & (residual < 50))
        #useInd=np.where((self.alsGroundClean.astype(float) > 0) & (self.gediQualClean == 1))                  # Extra residual threshold for GRSM data
        #useInd=np.where((self.alsGroundClean.astype(float) > 0) & (residual < 50) & (self.gediQualClean == 1))       # Default is all data, otherwise
        useInd=np.where((self.alsGroundClean.astype(float) > 0) & (residual < 50) & (self.gediSolarClean.astype(float) < 0.0) & (self.gediQualClean == 1))
        #print('useInd length',useInd.shape[0])
        print('useInd length',len(useInd[0]))

        preferredRMSE=sqrt(np.sum(residual[useInd]**2)/residual[useInd].shape[0])
        preferredBias=(np.sum(residual[useInd]))/residual[useInd].shape[0]
        '''ground1RMSE=sqrt(np.sum(residual1[useInd]**2)/residual1[useInd].shape[0])
        ground1Bias=(np.sum(residual1[useInd]))/residual1[useInd].shape[0]
        ground2RMSE=sqrt(np.sum(residual2[useInd]**2)/residual2[useInd].shape[0])
        ground2Bias=(np.sum(residual2[useInd]))/residual2[useInd].shape[0]
        ground3RMSE=sqrt(np.sum(residual3[useInd]**2)/residual3[useInd].shape[0])
        ground3Bias=(np.sum(residual3[useInd]))/residual3[useInd].shape[0]
        ground4RMSE=sqrt(np.sum(residual4[useInd]**2)/residual4[useInd].shape[0])
        ground4Bias=(np.sum(residual4[useInd]))/residual4[useInd].shape[0]
        ground5RMSE=sqrt(np.sum(residual5[useInd]**2)/residual5[useInd].shape[0])
        ground5Bias=(np.sum(residual5[useInd]))/residual5[useInd].shape[0]
        ground6RMSE=sqrt(np.sum(residual6[useInd]**2)/residual6[useInd].shape[0])
        ground6Bias=(np.sum(residual6[useInd]))/residual6[useInd].shape[0]'''
        print('RMSE:',preferredRMSE)
        print('Bias:',preferredBias)
        '''print('Algorithm 1 RMSE:',ground1RMSE)
        print('Algorithm 1 Bias:',ground1Bias)
        print('Algorithm 2 RMSE:',ground2RMSE)
        print('Algorithm 2 Bias:',ground2Bias)
        print('Algorithm 3 RMSE:',ground3RMSE)
        print('Algorithm 3 Bias:',ground3Bias)
        print('Algorithm 4 RMSE:',ground4RMSE)
        print('Algorithm 4 Bias:',ground4Bias)
        print('Algorithm 5 RMSE:',ground5RMSE)
        print('Algorithm 5 Bias:',ground5Bias)
        print('Algorithm 6 RMSE:',ground6RMSE)
        print('Algorithm 6 Bias:',ground6Bias)

        counts=np.bincount(self.algorithmClean[useInd])
        print(counts)'''

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        '''plt.rcParams['figure.figsize']=(8,6)
        plt.rcParams['xtick.labelsize']=16
        plt.rcParams['ytick.labelsize']=16
        plt.rcParams['axes.labelsize']=18
        plt.rcParams['axes.labelpad']=8.0
        plt.rcParams['axes.titlesize']=20
        #print(plt.rcParams)

        plt.plot(self.alsGroundClean[useInd].astype(float),self.gediGroundClean[useInd],'o',markersize=1)
        plt.plot((25,225),(25,225))
        plt.title(args.site + ' ALS vs GEDI Ground Elevation')
        plt.xlabel('ALS Ground Elevation (m)')
        plt.ylabel('GEDI Ground Elevation (m)')
        plt.xlim([25,225])
        plt.ylim([25,225])
        plt.savefig(args.outRoot+'ALSvsGEDI.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

        if args.dens:
            nbins=1000
            x=self.alsGroundClean[useInd].astype(float)
            y=self.gediGroundClean[useInd]
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.title(args.site + ' ALS vs GEDI Ground Elevation')
            plt.xlabel('ALS Ground Elevation (m)')
            plt.ylabel('GEDI Ground Elevation (m)')
            plt.colorbar()
            plt.savefig(args.outRoot+'ALSvsGEDIDensity.png')
            plt.close()
            plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Slope Plots
        plt.plot(self.slopeClean[useInd].astype(float),residual[useInd],'o',markersize=1)
        #plt.plot((580,1100),(580,1100))
        plt.title(args.site + ' Slope vs Residual')
        plt.xlabel('Slope')
        plt.ylabel('Ground Residual (m)')
        #plt.xlim([0,30])
        #plt.ylim([-5,25])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'SlopeVsRes.png',dpi=300)
        plt.close()
        plt.clf()

        if args.box:
            x1=residual[np.where((self.slopeClean[useInd].astype(float) < 10))]
            x2=residual[np.where((self.slopeClean[useInd].astype(float) > 10) & (self.slopeClean[useInd].astype(float) < 20))]
            x3=residual[np.where((self.slopeClean[useInd].astype(float) > 20) & (self.slopeClean[useInd].astype(float) < 30))]
            x4=residual[np.where((self.slopeClean[useInd].astype(float) > 30) & (self.slopeClean[useInd].astype(float) < 40))]
            #x5=residual[np.where((self.slopeClean[useInd].astype(float) > 40) & (self.slopeClean[useInd].astype(float) < 50))]
            x5=residual[np.where((self.slopeClean[useInd].astype(float) > 40))]
            data=[x1, x2, x3, x4, x5]

            plt.boxplot(data,sym='')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title(args.site + ' Slope vs. Residual')
            plt.ylabel('Ground Residual (m)')
            plt.ylim([-30,40])
            plt.xlabel('Slope (degrees)')
            plt.xticks(np.arange(1,6,step=1),['5', '15', '25', '35', '45'])
            plt.savefig(args.outRoot+'slopeBox.png',dpi=300)
            plt.close()
            plt.clf()

        if args.hist:
            plt.hist(self.slopeClean[useInd].astype(float),bins=np.arange(0.0,70.0+5,5),color='grey',edgecolor='white',linewidth=2)
            plt.title(args.site + ' Slope Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('Slope (degrees)')
            plt.savefig(args.outRoot+'slopeHist.png',dpi=300)
            plt.close()
            plt.clf()

        if args.dens:
            nbins=1000
            x=self.slopeClean[useInd].astype(float)
            y=residual[useInd]
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.axhline(y=0.0,color='r')
            plt.title('Slope vs Residual')
            plt.xlabel('Slope')
            plt.ylabel('Ground Residual (m)')
            plt.colorbar()
            plt.savefig(args.outRoot+'SlopeVsResDensity.png')
            plt.close()
            plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Cover Plots
        plt.plot(self.alsCoverClean[useInd].astype(float),residual[useInd],'o',markersize=1)
        #plt.plot((580,1100),(580,1100))
        plt.title(args.site + ' Canopy Cover vs Residual')
        plt.xlabel('Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        #plt.xlim([0.86,1.0])
        #plt.ylim([-5,25])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'CoverVsRes.png',dpi=300)
        plt.close()
        plt.clf()

        if args.box:
            x1=residual[np.where((self.alsCoverClean[useInd].astype(float) < 0.1))]
            x2=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.1) & (self.alsCoverClean[useInd].astype(float) < 0.2))]
            x3=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.2) & (self.alsCoverClean[useInd].astype(float) < 0.3))]
            x4=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.3) & (self.alsCoverClean[useInd].astype(float) < 0.4))]
            x5=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.4) & (self.alsCoverClean[useInd].astype(float) < 0.5))]
            x6=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.5) & (self.alsCoverClean[useInd].astype(float) < 0.6))]
            x7=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.6) & (self.alsCoverClean[useInd].astype(float) < 0.7))]
            x8=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.7) & (self.alsCoverClean[useInd].astype(float) < 0.8))]
            x9=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.8) & (self.alsCoverClean[useInd].astype(float) < 0.9))]
            x10=residual[np.where((self.alsCoverClean[useInd].astype(float) > 0.9) & (self.alsCoverClean[useInd].astype(float) < 1.00))]
            data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

            plt.boxplot(data,sym='')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title(args.site + ' Canopy Cover vs. Residual')
            plt.ylabel('Ground Residual (m)')
            plt.ylim([-30,40])
            plt.xlabel('Canopy Cover (fraction)')
            plt.xticks(np.arange(1,11,step=1),['0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
            plt.savefig(args.outRoot+'coverBox.png',dpi=300)
            plt.close()
            plt.clf()

        if args.hist:
            plt.hist(self.alsCoverClean[useInd].astype(float),bins=np.arange(0.0,1.0+0.05,0.05),color='grey',edgecolor='white',linewidth=2)
            plt.title(args.site + ' Canopy Cover Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('Canopy Cover (fraction)')
            plt.savefig(args.outRoot+'coverHist.png',dpi=300)
            plt.close()
            plt.clf()

        if args.dens:
            nbins=1000
            x=self.alsCoverClean[useInd].astype(float)
            y=residual[useInd]
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.axhline(y=0.0,color='r')
            plt.title(args.site + ' Canopy Cover vs Residual')
            plt.xlabel('Canopy Cover')
            plt.ylabel('Ground Residual (m)')
            plt.colorbar()
            plt.savefig(args.outRoot+'CoverVsResDensity.png')
            plt.close()
            plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Height Plots
        plt.plot(self.alsRH95Clean[useInd].astype(float),residual[useInd],'o',markersize=1)
        #plt.plot((580,1100),(580,1100))
        plt.title(args.site + ' Height vs Residual')
        plt.xlabel('RH95 (m)')
        plt.ylabel('Ground Residual (m)')
        #plt.xlim([10,50])
        #plt.ylim([-5,25])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'HeightVsRes.png',dpi=300)
        plt.close()
        plt.clf()

        if args.box:
            x1=residual[np.where((self.alsRH95Clean[useInd].astype(float) < 10))]
            x2=residual[np.where((self.alsRH95Clean[useInd].astype(float) > 10) & (self.alsRH95Clean[useInd].astype(float) < 20))]
            x3=residual[np.where((self.alsRH95Clean[useInd].astype(float) > 20) & (self.alsRH95Clean[useInd].astype(float) < 30))]
            x4=residual[np.where((self.alsRH95Clean[useInd].astype(float) > 30) & (self.alsRH95Clean[useInd].astype(float) < 40))]
            #x5=residual[np.where((self.alsRH95Clean[useInd].astype(float) > 40) & (self.alsRH95Clean[useInd].astype(float) < 50))]
            #x6=residual[np.where((self.alsRH95Clean[useInd].astype(float) > 50) & (self.alsRH95Clean[useInd].astype(float) < 60))]
            x5=residual[np.where((self.alsRH95Clean[useInd].astype(float) > 40))]
            data=[x1, x2, x3, x4, x5]

            plt.boxplot(data,sym='')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title(args.site + ' Canopy Height vs. Residual')
            plt.ylabel('Ground Residual (m)')
            plt.ylim([-30,40])
            plt.xlabel('Canopy Height (m)')
            plt.xticks(np.arange(1,6,step=1),['5', '15', '25', '35', '45'])
            plt.savefig(args.outRoot+'heightBox.png',dpi=300)
            plt.close()
            plt.clf()

        if args.hist:
            plt.hist(self.alsRH95Clean[useInd].astype(float),bins=np.arange(0.0,60.0+5,5),color='grey',edgecolor='white',linewidth=2)
            plt.title(args.site + ' Canopy Height Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('Canopy Height (m)')
            plt.savefig(args.outRoot+'heightHist.png',dpi=300)
            plt.close()
            plt.clf()

        if args.dens:
            nbins=1000
            x=self.alsRH95Clean[useInd].astype(float)
            y=residual[useInd]
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.axhline(y=0.0,color='r')
            plt.title(args.site + ' Height vs Residual')
            plt.xlabel('RH95 (m)')
            plt.ylabel('Ground Residual (m)')
            plt.colorbar()
            plt.savefig(args.outRoot+'HeightVsResDensity.png')
            plt.close()
            plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Ground Overlap Plots
        plt.plot(self.overlapClean[useInd].astype(float),residual[useInd],'o',markersize=1)
        #plt.plot((580,1100),(580,1100))
        plt.title(args.site + ' Ground Overlap vs Residual')
        plt.xlabel('Ground Overlap')
        plt.ylabel('Ground Residual (m)')
        #plt.xlim([0.2,1.0])
        #plt.ylim([-5,25])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'OverlapVsRes.png',dpi=300)
        plt.close()
        plt.clf()

        if args.box:
            x1=residual[np.where((self.overlapClean[useInd].astype(float) < 0.15))]
            x2=residual[np.where((self.overlapClean[useInd].astype(float) > 0.15) & (self.overlapClean[useInd].astype(float) < 0.25))]
            x3=residual[np.where((self.overlapClean[useInd].astype(float) > 0.25) & (self.overlapClean[useInd].astype(float) < 0.35))]
            x4=residual[np.where((self.overlapClean[useInd].astype(float) > 0.35) & (self.overlapClean[useInd].astype(float) < 0.45))]
            x5=residual[np.where((self.overlapClean[useInd].astype(float) > 0.45) & (self.overlapClean[useInd].astype(float) < 0.55))]
            x6=residual[np.where((self.overlapClean[useInd].astype(float) > 0.55) & (self.overlapClean[useInd].astype(float) < 0.65))]
            x7=residual[np.where((self.overlapClean[useInd].astype(float) > 0.65) & (self.overlapClean[useInd].astype(float) < 0.75))]
            x8=residual[np.where((self.overlapClean[useInd].astype(float) > 0.75) & (self.overlapClean[useInd].astype(float) < 0.85))]
            x9=residual[np.where((self.overlapClean[useInd].astype(float) > 0.85) & (self.overlapClean[useInd].astype(float) < 0.95))]
            x10=residual[np.where((self.overlapClean[useInd].astype(float) > 0.95) & (self.overlapClean[useInd].astype(float) < 1.00))]
            data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

            plt.boxplot(data,sym='')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title(args.site + ' Ground Overlap vs. Residual')
            plt.ylabel('Ground Residual (m)')
            plt.ylim([-30,40])
            plt.xlabel('Ground Overlap')
            plt.xticks(np.arange(1,11,step=1),['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'])
            plt.savefig(args.outRoot+'overlapBox.png',dpi=300)
            plt.close()
            plt.clf()

        if args.hist:
            plt.hist(self.overlapClean[useInd].astype(float),bins=np.arange(0.0,1.0+0.05,0.05),color='grey',edgecolor='white',linewidth=2)
            plt.title(args.site + ' Ground Overlap Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('Ground Overlap')
            plt.savefig(args.outRoot+'overlapHist.png',dpi=300)
            plt.close()
            plt.clf()

        if args.dens:
            nbins=1000
            x=self.overlapClean[useInd].astype(float)
            y=residual[useInd]
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.axhline(y=0.0,color='r')
            plt.title(args.site + ' Ground Overlap vs Residual')
            plt.xlabel('Ground Overlap')
            plt.ylabel('Ground Residual (m)')
            plt.colorbar()
            plt.savefig(args.outRoot+'OverlapVsResDensity.png')
            plt.close()
            plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Beam Sensitivity Plots
        plt.plot(self.gediSensClean[useInd],residual[useInd],'o',markersize=1)
        #plt.plot((580,1100),(580,1100))
        plt.title(args.site + ' Sensitivity vs Residual')
        plt.xlabel('Beam Sensitivity')
        plt.ylabel('Ground Residual (m)')
        #plt.xlim([0.88,1.0])
        #plt.ylim([-5,25])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'BeamSensVsRes.png',dpi=300)
        plt.close()
        plt.clf()

        if args.box:
            x1=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.90) & (self.gediSensClean[useInd].astype(float) < 0.91))]
            x2=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.91) & (self.gediSensClean[useInd].astype(float) < 0.92))]
            x3=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.92) & (self.gediSensClean[useInd].astype(float) < 0.93))]
            x4=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.93) & (self.gediSensClean[useInd].astype(float) < 0.94))]
            x5=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.94) & (self.gediSensClean[useInd].astype(float) < 0.95))]
            x6=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.95) & (self.gediSensClean[useInd].astype(float) < 0.96))]
            x7=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.96) & (self.gediSensClean[useInd].astype(float) < 0.97))]
            x8=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.97) & (self.gediSensClean[useInd].astype(float) < 0.98))]
            x9=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.98) & (self.gediSensClean[useInd].astype(float) < 0.99))]
            x10=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.99) & (self.gediSensClean[useInd].astype(float) < 1.00))]

            x1=residual[np.where((self.gediSensClean[useInd].astype(float) < 0.82))]
            x2=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.82) & (self.gediSensClean[useInd].astype(float) < 0.84))]
            x3=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.84) & (self.gediSensClean[useInd].astype(float) < 0.86))]
            x4=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.86) & (self.gediSensClean[useInd].astype(float) < 0.88))]
            x5=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.88) & (self.gediSensClean[useInd].astype(float) < 0.90))]
            x6=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.90) & (self.gediSensClean[useInd].astype(float) < 0.92))]
            x7=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.92) & (self.gediSensClean[useInd].astype(float) < 0.94))]
            x8=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.94) & (self.gediSensClean[useInd].astype(float) < 0.96))]
            x9=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.96) & (self.gediSensClean[useInd].astype(float) < 0.98))]
            x10=residual[np.where((self.gediSensClean[useInd].astype(float) > 0.98) & (self.gediSensClean[useInd].astype(float) < 1.00))]

            data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

            plt.boxplot(data,sym='')
            plt.axhline(y=0.0,color='r',linestyle='--')
            plt.title(args.site + ' Beam Sensitivity vs. Residual')
            plt.ylabel('Ground Residual (m)')
            plt.ylim([-30,40])
            plt.xlabel('Beam Sensitivity')
            #plt.xticks(np.arange(1,11,step=1),['0.90', '0.91', '0.92', '0.93', '0.94', '0.95', '0.96', '0.97', '0.98', '0.99'])
            plt.xticks(np.arange(1,11,step=1),['0.81', '0.83', '0.85', '0.87', '0.89', '0.91', '0.93', '0.95', '0.97', '0.99'])
            plt.savefig(args.outRoot+'beamSensBox.png',dpi=300)
            plt.close()
            plt.clf()

        if args.hist:
            plt.hist(self.gediSensClean[useInd].astype(float),bins=np.arange(0.7,1.0+0.01,0.01),color='grey',edgecolor='white',linewidth=2)
            plt.title(args.site + ' Beam Sensitivity Distribution')
            plt.ylabel('Frequency')
            plt.xlabel('Beam Sensitivity')
            plt.savefig(args.outRoot+'beamSensHist.png',dpi=300)
            plt.close()
            plt.clf()

        if args.dens:
            nbins=1000
            x=self.gediSensClean[useInd]
            y=residual[useInd]
            k=kde.gaussian_kde([x,y])
            xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(),yi.flatten()]))
            plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
            plt.axhline(y=0.0,color='r')
            plt.title(args.site + ' Sensitivity vs Residual')
            plt.xlabel('Beam Sensitivity')
            plt.ylabel('Ground Residual (m)')
            plt.colorbar()
            plt.savefig(args.outRoot+'BeamSensVsResDensity.png')
            plt.close()
            plt.clf()'''

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def copyPlots(self):
        '''Identify a subset of waveforms by shotnumber and copy those plots'''

        self.alsShotClean=np.empty(0)
        self.alsGroundClean=np.empty(0,dtype=float)
        self.slopeClean=np.empty(0,dtype=float)
        self.alsCoverClean=np.empty(0,dtype=float)
        self.alsRH95Clean=np.empty(0,dtype=float)
        self.overlapClean=np.empty(0,dtype=float)
        self.gediShotClean=np.empty(0)
        self.gediGroundClean=np.empty(0,dtype=float)
        self.gediSensClean=np.empty(0,dtype=float)
        self.gediSolarClean=np.empty(0,dtype=float)
        self.gediQualClean=np.empty(0,dtype=int)

        # Match each record in the ALS metric file to the corresponding L2A GEDI record
        for i in range(alsShotMaster.shape[0]):
            for j in range(gediShotMaster.shape[0]):
                if alsShotMaster[i] == gediShotMaster[j]:
                    self.alsShotClean=np.append(self.alsShotClean,alsShotMaster[i])
                    self.alsGroundClean=np.append(self.alsGroundClean,alsGroundMaster[i])
                    self.slopeClean=np.append(self.slopeClean,slopeMaster[i])
                    self.alsCoverClean=np.append(self.alsCoverClean,alsCoverMaster[i])
                    self.alsRH95Clean=np.append(self.alsRH95Clean,alsRH95Master[i])
                    self.overlapClean=np.append(self.overlapClean,overlapMaster[i])
                    self.gediShotClean=np.append(self.gediShotClean,gediShotMaster[j])
                    self.gediGroundClean=np.append(self.gediGroundClean,gediGroundMaster[j])
                    self.gediSensClean=np.append(self.gediSensClean,gediSensMaster[j])
                    self.gediSolarClean=np.append(self.gediSolarClean,gediSolarMaster[j])
                    self.gediQualClean=np.append(self.gediQualClean,gediQualMaster[j])

        residual=np.subtract(self.gediGroundClean,self.alsGroundClean.astype(float))

        useInd=np.where((self.slopeClean.astype(float) > 45.0) & (self.gediQualClean == 1))
        shotnumbers=self.gediShotClean[useInd]

        print(shotnumbers.shape[0],'matching shotnumbers identified')

        outFile='slope60.sh'
        root='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l1b/pulse/waveforms/'
        with open(outFile,'a') as outf:
            outf.write('#!/usr/bin/env bash\n\n')
            for shot in shotnumbers:
                outf.write('cp {0}plot.beam.{1}.png {2}slope60/\n\n'.format(root,shot,root))




#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    test=compareGround()
    if args.als:
        inputs='sim*.txt'
        path=os.path.join(args.als,inputs)
        fileList=glob.glob(path)

        alsShotMaster=np.empty(0)
        alsGroundMaster=np.empty(0,dtype=float)
        slopeMaster=np.empty(0,dtype=float)
        alsCoverMaster=np.empty(0,dtype=float)
        alsRH95Master=np.empty(0,dtype=float)
        overlapMaster=np.empty(0,dtype=float)

        for file in fileList:
            test.readALS(file)
            alsShotMaster=np.append(alsShotMaster,test.alsShot)
            alsGroundMaster=np.append(alsGroundMaster,test.alsGround)
            slopeMaster=np.append(slopeMaster,test.slope)
            alsCoverMaster=np.append(alsCoverMaster,test.alsCover)
            alsRH95Master=np.append(alsRH95Master,test.alsRH95)
            overlapMaster=np.append(overlapMaster,test.overlap)

    else:
        test.readALS(args.alsFile)

    if args.gedi:
        inputs='*.h5'
        path=os.path.join(args.gedi,inputs)
        fileList=glob.glob(path)

        gediShotMaster=np.empty(0,dtype=str)
        gediGroundMaster=np.empty(0,dtype=float)
        gediSensMaster=np.empty(0,dtype=float)
        gediSolarMaster=np.empty(0,dtype=float)
        gediQualMaster=np.empty(0,dtype=int)
        algorithmMaster=np.empty(0,dtype=int)
        ground1Master=np.empty(0,dtype=float)
        ground2Master=np.empty(0,dtype=float)
        ground3Master=np.empty(0,dtype=float)
        ground4Master=np.empty(0,dtype=float)
        ground5Master=np.empty(0,dtype=float)
        ground6Master=np.empty(0,dtype=float)

        for file in fileList:
            test.readGEDI(file)
            gediShotMaster=np.append(gediShotMaster,test.gediShot)
            gediGroundMaster=np.append(gediGroundMaster,test.gediGround)
            gediSensMaster=np.append(gediSensMaster,test.gediSens)
            gediSolarMaster=np.append(gediSolarMaster,test.gediSolar)
            gediQualMaster=np.append(gediQualMaster,test.gediQual)
            algorithmMaster=np.append(algorithmMaster,test.algorithm)
            ground1Master=np.append(ground1Master,test.ground1)
            ground2Master=np.append(ground2Master,test.ground2)
            ground3Master=np.append(ground3Master,test.ground3)
            ground4Master=np.append(ground4Master,test.ground4)
            ground5Master=np.append(ground5Master,test.ground5)
            ground6Master=np.append(ground6Master,test.ground6)

    else:
        test.readGEDI(args.gediFile)

    if args.plots:
        test.plotData()

    if args.waveforms:
        test.copyPlots()
