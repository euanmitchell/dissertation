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
        self.alsGround=np.empty(0,dtype=float)
        self.alsSlope=np.empty(0,dtype=float)
        self.alsCover=np.empty(0,dtype=float)
        # Do this with numpy.loadtxt() instead?
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

    def readGEDI(self,input):
        '''Read the GEDI L2A file'''

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
        self.gediShotClean=np.empty(0)
        self.gediGroundClean=np.empty(0,dtype=float)

        # Match each record in the ALS metric file to the corresponding L2A GEDI record
        for i in range(self.alsShot.shape[0]):
            for j in range(self.gediShot.shape[0]):         # Replace this line with a np.where() query?
                if self.alsShot[i] == self.gediShot[j]:
                    self.alsShotClean=np.append(self.alsShotClean,self.alsShot[i])
                    self.alsGroundClean=np.append(self.alsGroundClean,self.alsGround[i])
                    self.gediShotClean=np.append(self.gediShotClean,self.gediShot[j])
                    self.gediGroundClean=np.append(self.gediGroundClean,self.gediGround[j])

        residual=(self.gediGroundClean-self.alsGroundClean)

        plt.plot(self.alsGroundClean.astype(float),self.gediGroundClean,'o')
        plt.plot((580,1100),(580,1100))
        plt.title('ALS vs GEDI Ground Elevation')
        plt.xlabel('ALS Ground Elevation (m)')
        plt.ylabel('GEDI Ground Elevation (m)')
        #plt.savefig('name')
        #plt.close()
        #plt.clf()
        plt.show()

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    test=compareGround()
    test.readALS(args.alsFile)
    test.readGEDI(args.gediFile)
    test.plotData()
