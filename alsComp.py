#######################################
# Script to plot Version 2.0 preferred
# ground estimate against collocated
# ALS ground estimate.
#######################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse

#######################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify input ALS and GEDI data files'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--als', dest='als', type=str, default='',
        help=('The path to the directory containing the ALS metric files.\nDefault is not set'))
    p.add_argument('--alsFile', dest='alsFile', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l1b/subset/metric/GEDI01_B_2019262084455_O04355_03_T04290_02_005_01_V002.metric.txt',
        help=('The path to a single ALS metric file.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l1b/subset/metric/GEDI01_B_2019262084455_O04355_03_T04290_02_005_01_V002.metric.txt'))
    p.add_argument('--gedi', dest='gedi', type=str, default='',
        help=('The path to the directory containing the GEDI data.\nDefault is /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l2a/processed_GEDI02_A_2019262084455_O04355_03_T04290_02_003_01_V002.h5'))
    p.add_argument('--gediFile', dest='gediFile', type=str, default='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l2a/processed_GEDI02_A_2019262084455_O04355_03_T04290_02_003_01_V002.h5',
        help=('The path to a single GEDI data file.\nDefault is not set'))
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
        self.beamSens=np.empty(0,dtype=float)       # 'blairSense'

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
                    self.beamSens=np.append(self.beamSens,data[112])


    def readGEDI(self,input):
        '''Read the GEDI L2A file'''

        '''GET GROUND ELEVATION ARRAYS FOR EACH OF THE SIX ALGORITHMS'''

        print('Reading GEDI L2A file',input)
        beamlist=['BEAM0101','BEAM1000','BEAM0010','BEAM0000','BEAM0011','BEAM0110','BEAM0001','BEAM1011']
        # Create empty array for shotnumbers
        self.gediShot=np.empty(0,dtype=str)
        # Create empty array for GEDI ground estimate
        self.gediGround=np.empty(0,dtype=float)
        f=h5py.File(input,'r')
        for beam in beamlist:
            # Need to handle empty beams
            try:
                self.gediShot=np.append(self.gediShot,f[beam]['shot_number'])
                self.gediGround=np.append(self.gediGround,f[beam]['elev_lowestmode'])
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
        self.beamSensClean=np.empty(0,dtype=float)
        self.gediShotClean=np.empty(0)
        self.gediGroundClean=np.empty(0,dtype=float)

        # Match each record in the ALS metric file to the corresponding L2A GEDI record
        for i in range(self.alsShot.shape[0]):
            for j in range(self.gediShot.shape[0]):         # Replace this line with a np.where() query?
                if self.alsShot[i] == self.gediShot[j]:
                    self.alsShotClean=np.append(self.alsShotClean,self.alsShot[i])
                    self.alsGroundClean=np.append(self.alsGroundClean,self.alsGround[i])
                    self.slopeClean=np.append(self.slopeClean,self.slope[i])
                    self.alsCoverClean=np.append(self.alsCoverClean,self.alsCover[i])
                    self.alsRH95Clean=np.append(self.alsRH95Clean,self.alsRH95[i])
                    self.overlapClean=np.append(self.overlapClean,self.overlap[i])
                    self.beamSensClean=np.append(self.beamSensClean,self.beamSens[i])
                    self.gediShotClean=np.append(self.gediShotClean,self.gediShot[j])
                    self.gediGroundClean=np.append(self.gediGroundClean,self.gediGround[j])

        residual=np.subtract(self.gediGroundClean,self.alsGroundClean.astype(float))

        plt.plot(self.alsGroundClean.astype(float),self.gediGroundClean,'o')
        plt.plot((580,1100),(580,1100))
        plt.title('ALS vs GEDI Ground Elevation')
        plt.xlabel('ALS Ground Elevation (m)')
        plt.ylabel('GEDI Ground Elevation (m)')
        plt.savefig('code/ALSvsGEDI.png')
        plt.close()
        plt.clf()
        #plt.show()

        plt.plot(self.slopeClean.astype(float),residual,'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Slope vs Residual')
        plt.xlabel('Slope')
        plt.ylabel('Ground Residual (m)')
        plt.savefig('code/SlopeVsRes.png')
        plt.close()
        plt.clf()

        plt.plot(self.alsCoverClean.astype(float),residual,'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Canopy Cover vs Residual')
        plt.xlabel('Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        plt.savefig('code/CoverVsRes.png')
        plt.close()
        plt.clf()

        plt.plot(self.alsRH95Clean.astype(float),residual,'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Height vs Residual')
        plt.xlabel('RH95 (m)')
        plt.ylabel('Ground Residual (m)')
        plt.savefig('code/HeightVsRes.png')
        plt.close()
        plt.clf()

        plt.plot(self.overlapClean.astype(float),residual,'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Ground Overlap vs Residual')
        plt.xlabel('Ground Overlap')
        plt.ylabel('Ground Residual (m)')
        plt.savefig('code/OverlapVsRes.png')
        plt.close()
        plt.clf()

        plt.plot(self.beamSensClean.astype(float),residual,'o')
        #plt.plot((580,1100),(580,1100))
        plt.title('Sensitivity vs Residual')
        plt.xlabel('Beam Sensitivity')
        plt.ylabel('Ground Residual (m)')
        plt.savefig('code/BeamSensVsRes.png')
        plt.close()
        plt.clf()

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    test=compareGround()
    test.readALS(args.alsFile)
    test.readGEDI(args.gediFile)
    test.plotData()
