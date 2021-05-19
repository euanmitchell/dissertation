#######################################
# Script to plot Version 2.0 preferred
# ground estimate against collocated
# ALS ground estimate.
#######################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
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
    p.add_argument('--als', dest='als', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l1b/pulse/metric/',
        help=('The path to the directory containing the ALS metric files.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l1b/pulse/metric/'))
    p.add_argument('--alsFile', dest='alsFile', type=str, default='',
        help=('The path to a single ALS metric file.\nDefault is not set'))
    p.add_argument('--gedi', dest='gedi', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l2a/',
        help=('The path to the directory containing the GEDI data.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l2a/'))
    p.add_argument('--gediFile', dest='gediFile', type=str, default='',
        help=('The path to a single GEDI data file.\nDefault is not set'))
    p.add_argument('--beams', dest='beams', type=str, default='all',
        help=('Which beams to use: "all", "power", or "coverage".\nDefault is all'))
    p.add_argument('--outRoot', dest='outRoot', type=str, default='../data/figures/',
        help=('Output root for plots.\nDefault is "../data/figures/"'))
    p.add_argument('--plots', dest='plots', action='store_true', default=False,
        help=('Call the plotData() method to make comparision plots'))
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

        '''GET GROUND ELEVATION ARRAYS FOR EACH OF THE SIX ALGORITHMS'''

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
        f=h5py.File(input,'r')
        for beam in beamlist:
            # Need to handle empty beams
            try:
                self.gediShot=np.append(self.gediShot,f[beam]['shot_number'])
                self.gediGround=np.append(self.gediGround,f[beam]['elev_lowestmode'])
                self.gediSens=np.append(self.gediSens,f[beam]['sensitivity'])
                self.gediSolar=np.append(self.gediSolar,f[beam]['solar_elevation'])
                self.gediQual=np.append(self.gediQual,f[beam]['quality_flag'])
                # BEAM SENSITIVITY
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

        # Specify a subset of data to plot here with useInd?
        #useInd=np.where((self.gediQualClean == 1))       # Default is all data, otherwise
        useInd=np.where((self.gediQualClean == 1) & (self.gediSolarClean.astype(float) > 0.0))
        print('useInd length',len(useInd[0]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        plt.plot(self.alsGroundClean[useInd].astype(float),self.gediGroundClean[useInd],'o')
        plt.plot((200,1300),(200,1300))
        plt.title('ALS vs GEDI Ground Elevation')
        plt.xlabel('ALS Ground Elevation (m)')
        plt.ylabel('GEDI Ground Elevation (m)')
        plt.xlim([200,1300])
        plt.ylim([200,1300])
        plt.savefig(args.outRoot+'ALSvsGEDI.png')
        plt.close()
        plt.clf()
        #plt.show()

        nbins=1000
        x=self.alsGroundClean[useInd].astype(float)
        y=self.gediGroundClean[useInd]
        k=kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(),yi.flatten()]))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        plt.title('ALS vs GEDI Ground Elevation')
        plt.xlabel('ALS Ground Elevation (m)')
        plt.ylabel('GEDI Ground Elevation (m)')
        plt.colorbar()
        plt.savefig(args.outRoot+'ALSvsGEDIDensity.png')
        plt.close()
        plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        plt.plot(self.slopeClean[useInd].astype(float),residual[useInd],'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Slope vs Residual')
        plt.xlabel('Slope')
        plt.ylabel('Ground Residual (m)')
        plt.xlim([0,70])
        plt.ylim([-30,60])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'SlopeVsRes.png')
        plt.close()
        plt.clf()

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

        plt.plot(self.alsCoverClean[useInd].astype(float),residual[useInd],'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Canopy Cover vs Residual')
        plt.xlabel('Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        plt.xlim([0.0,1.0])
        plt.ylim([-30,60])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'CoverVsRes.png')
        plt.close()
        plt.clf()

        nbins=1000
        x=self.alsCoverClean[useInd].astype(float)
        y=residual[useInd]
        k=kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(),yi.flatten()]))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        plt.axhline(y=0.0,color='r')
        plt.title('Canopy Cover vs Residual')
        plt.xlabel('Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        plt.colorbar()
        plt.savefig(args.outRoot+'CoverVsResDensity.png')
        plt.close()
        plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        plt.plot(self.alsRH95Clean[useInd].astype(float),residual[useInd],'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Height vs Residual')
        plt.xlabel('RH95 (m)')
        plt.ylabel('Ground Residual (m)')
        plt.xlim([0,80])
        plt.ylim([-30,60])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'HeightVsRes.png')
        plt.close()
        plt.clf()

        nbins=1000
        x=self.alsRH95Clean[useInd].astype(float)
        y=residual[useInd]
        k=kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(),yi.flatten()]))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        plt.axhline(y=0.0,color='r')
        plt.title('Height vs Residual')
        plt.xlabel('RH95 (m)')
        plt.ylabel('Ground Residual (m)')
        plt.colorbar()
        plt.savefig(args.outRoot+'HeightVsResDensity.png')
        plt.close()
        plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        plt.plot(self.overlapClean[useInd].astype(float),residual[useInd],'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Ground Overlap vs Residual')
        plt.xlabel('Ground Overlap')
        plt.ylabel('Ground Residual (m)')
        plt.xlim([0.0,1.0])
        plt.ylim([-30,60])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'OverlapVsRes.png')
        plt.close()
        plt.clf()

        nbins=1000
        x=self.overlapClean[useInd].astype(float)
        y=residual[useInd]
        k=kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(),yi.flatten()]))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        plt.axhline(y=0.0,color='r')
        plt.title('Ground Overlap vs Residual')
        plt.xlabel('Ground Overlap')
        plt.ylabel('Ground Residual (m)')
        plt.colorbar()
        plt.savefig(args.outRoot+'OverlapVsResDensity.png')
        plt.close()
        plt.clf()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        plt.plot(self.gediSensClean[useInd],residual[useInd],'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Sensitivity vs Residual')
        plt.xlabel('Beam Sensitivity')
        plt.ylabel('Ground Residual (m)')
        plt.xlim([0.9,1.0])
        plt.ylim([-30,60])
        plt.axhline(y=0.0,color='r')
        plt.savefig(args.outRoot+'BeamSensVsRes.png')
        plt.close()
        plt.clf()

        nbins=1000
        x=self.gediSensClean[useInd]
        y=residual[useInd]
        k=kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(),yi.flatten()]))
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
        plt.axhline(y=0.0,color='r')
        plt.title('Sensitivity vs Residual')
        plt.xlabel('Beam Sensitivity')
        plt.ylabel('Ground Residual (m)')
        plt.colorbar()
        plt.savefig(args.outRoot+'BeamSensVsResDensity.png')
        plt.close()
        plt.clf()

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
        inputs='*.txt'
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

        for file in fileList:
            test.readGEDI(file)
            gediShotMaster=np.append(gediShotMaster,test.gediShot)
            gediGroundMaster=np.append(gediGroundMaster,test.gediGround)
            gediSensMaster=np.append(gediSensMaster,test.gediSens)
            gediSolarMaster=np.append(gediSolarMaster,test.gediSolar)
            gediQualMaster=np.append(gediQualMaster,test.gediQual)

    else:
        test.readGEDI(args.gediFile)

    if args.plots:
        test.plotData()

    if args.waveforms:
        test.copyPlots()
