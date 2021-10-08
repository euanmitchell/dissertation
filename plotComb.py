#######################################
# Script to read text files output by
# alsComp.py and read them in to make
# histogram and box plots for all
# three sites.
#######################################

import h5py
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import argparse

#######################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify input ALS and GEDI data files and program control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--outRoot', dest='outRoot', type=str, default='../data/newPaperFigs/',
        help=('Output root for plots.\nDefault is "../data/newPaperFigs/"'))
    p.add_argument('--hist', dest='hist', action='store_true', default=False,
        help=('Make histogram plots'))
    p.add_argument('--box', dest='box', action='store_true', default=False,
        help=('Make box plots'))
    cmdargs=p.parse_args()
    return cmdargs

#######################################

class comparisonPlots(object):
    '''A class to handle reading and plotting ALS and GEDI ground estimates for multiple sites'''

    def readData(self):
        '''Read the data arrays from the text files'''

        # The data files
        alsW='alsDataWREF.txt'
        alsG='alsDataGRSM.txt'
        alsL='alsDataLSBS.txt'
        gediW='gediDataWREF.txt'
        gediG='gediDataGRSM.txt'
        gediL='gediDataLSBS.txt'

        # ALS data column structure
        # 0 = ground elevation
        # 1 = ground slope
        # 2 = canopy cover
        # 3 = canopy height (RH95)
        # 4 = residual

        # GEDI data column structure
        # 0 = ground elevation
        # 1 = beam sensitivity
        # 2 = solar angle
        # 3 = quality flag
        # 4 = residual

        self.alsWREF=np.loadtxt(alsW)
        self.alsGRSM=np.loadtxt(alsG)
        self.alsLSBS=np.loadtxt(alsL)
        self.gediWREF=np.loadtxt(gediW)
        self.gediGRSM=np.loadtxt(gediG)
        self.gediLSBS=np.loadtxt(gediL)

        # Will be useful to have residual in ALS arrays - add extra column with column_stack
        self.alsWREF=np.column_stack((self.alsWREF,self.gediWREF[:,4]))
        self.alsGRSM=np.column_stack((self.alsGRSM,self.gediGRSM[:,4]))
        self.alsLSBS=np.column_stack((self.alsLSBS,self.gediLSBS[:,4]))

        print('WREF raw ALS data shape',self.alsWREF.shape)
        print('WREF raw GEDI data shape',self.gediWREF.shape)
        print('GRSM raw ALS data shape',self.alsGRSM.shape)
        print('GRSM raw GEDI data shape',self.gediGRSM.shape)
        print('LSBS raw ALS data shape',self.alsLSBS.shape)
        print('LSBS raw GEDI data shape',self.gediLSBS.shape)

        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50))
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50))
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50))

        print('Mean canopy height at WREF is',np.mean(self.alsWREF[useIndW,3]))
        print('Mean canopy height at GRSM is',np.mean(self.alsGRSM[useIndG,3]))
        print('Mean canopy height at LSBS is',np.mean(self.alsLSBS[useIndL,3]))

    def plotHistograms(self):
        '''Make the histogram plots from ALS data'''

        # Plot parameters
        plt.rcParams['figure.figsize']=(8,6)
        plt.rcParams['xtick.labelsize']=16
        plt.rcParams['ytick.labelsize']=16
        plt.rcParams['axes.labelsize']=18
        plt.rcParams['axes.labelpad']=8.0
        plt.rcParams['axes.titlesize']=20

        # Slope Histogram
        plt.hist(
            [self.alsWREF[:,1],self.alsGRSM[:,1],self.alsLSBS[:,1]],
            bins=np.arange(0.0,60.0+5,5),
            label=['WREF','GRSM','LSBS'],
            edgecolor='white'
        )
        plt.title('Slope Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Slope (degrees)')
        plt.legend(loc='upper right')
        plt.savefig(args.outRoot+'slopeHistAll.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

        # Cover Histogram
        plt.hist(
            [self.alsWREF[:,2],self.alsGRSM[:,2],self.alsLSBS[:,2]],
            bins=np.arange(0.0,1.0+0.05,0.05),
            label=['WREF','GRSM','LSBS'],
            edgecolor='white'
        )
        plt.title('Canopy Cover Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Canopy Cover (fraction)')
        plt.legend(loc='upper left')
        plt.savefig(args.outRoot+'coverHistAll.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

        # Height Histogram
        plt.hist(
            [self.alsWREF[:,3],self.alsGRSM[:,3],self.alsLSBS[:,3]],
            bins=np.arange(0.0,70.0+5,5),
            label=['WREF','GRSM','LSBS'],
            edgecolor='white'
        )
        plt.title('Canopy Height Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Canopy Height (m)')
        plt.legend(loc='upper right')
        plt.savefig(args.outRoot+'heightHistAll.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

        # Beam Sensitivity Histogram
        plt.hist(
            [self.gediWREF[:,1],self.gediGRSM[:,1],self.gediLSBS[:,1]],
            bins=np.arange(0.85,1.0+0.01,0.01),
            label=['WREF','GRSM','LSBS'],
            edgecolor='white'
        )
        plt.title('Beam Sensitivity Distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Beam Sensitivity (fraction)')
        plt.legend(loc='upper left')
        plt.savefig(args.outRoot+'sensHistAll.png',dpi=300)
        plt.close()
        plt.clf()
        #plt.show()

    def plotBoxplots(self):
        '''Make the boxplots of elevation residuals'''

        # Plot parameters
        plt.rcParams['figure.figsize']=(8,6)
        plt.rcParams['xtick.labelsize']=16
        plt.rcParams['ytick.labelsize']=16
        plt.rcParams['axes.labelsize']=18
        plt.rcParams['axes.labelpad']=8.0
        plt.rcParams['axes.titlesize']=20

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Slope Plots

        # Set the number of data bins
        slo=np.arange(0,60,5)
        increment=5
        print('SlopeBins:',slo)

        #WREF
        #useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50))
        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50) & (self.gediWREF[:,3] == 1.0))
        print('WREF useInd length',len(useIndW[0]))
        meansWREF=[]
        stdevsWREF=[]
        for i in slo:
            mean=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > i) & (self.alsWREF[useIndW,1] < (i+increment)))])
            meansWREF.append(mean)
            stdev=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > i) & (self.alsWREF[useIndW,1] < (i+increment)))])
            stdevsWREF.append(stdev)

        #GRSM
        #useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50))
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50) & (self.gediGRSM[:,3] == 1.0))
        print('GRSM useInd length',len(useIndG[0]))
        meansGRSM=[]
        stdevsGRSM=[]
        for i in slo:
            mean=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > i) & (self.alsGRSM[useIndG,1] < (i+increment)))])
            meansGRSM.append(mean)
            stdev=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > i) & (self.alsGRSM[useIndG,1] < (i+increment)))])
            stdevsGRSM.append(stdev)

        #LSBS
        #useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50))
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50) & (self.gediLSBS[:,3] == 1.0))
        print('LSBS useInd length',len(useIndL[0]))
        meansLSBS=[]
        stdevsLSBS=[]
        for i in slo:
            mean=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > i) & (self.alsLSBS[useIndL,1] < (i+increment)))])
            meansLSBS.append(mean)
            stdev=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > i) & (self.alsLSBS[useIndL,1] < (i+increment)))])
            stdevsLSBS.append(stdev)

        # Set the x-axis plot tick locations
        slopes=np.arange(2.5,62,5)
        print('Slope Ticks:',slopes)

        #plt.plot(self.alsWREF[:,1],self.gediWREF[:,4],'o',markersize=1)
        plt.plot(slopes,meansWREF,label='WREF')
        plt.plot(slopes,meansGRSM,label='GRSM')
        plt.plot(slopes,meansLSBS,label='LSBS')
        plt.fill_between(slopes,np.subtract(meansWREF,stdevsWREF),np.add(meansWREF,stdevsWREF),alpha=0.2)
        plt.fill_between(slopes,np.subtract(meansGRSM,stdevsGRSM),np.add(meansGRSM,stdevsGRSM),alpha=0.2)
        plt.fill_between(slopes,np.subtract(meansLSBS,stdevsLSBS),np.add(meansLSBS,stdevsLSBS),alpha=0.2)
        plt.title('Residual vs. Slope')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,30])
        plt.xlabel('Slope (degrees)')
        plt.xlim([0,60])
        plt.legend(loc='upper left')
        #plt.savefig(args.outRoot+'slopeResQual.png',dpi=300)
        #plt.close()
        #plt.clf()
        #print('Slope plot written to file '+args.outRoot+'slopeResQual.png')
        plt.show()

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Cover Plots

        # Set the number of data bins
        covs=np.arange(0,1.0,0.05)
        increment=0.05
        print('CoverBins:',covs)

        #WREF
        #useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50))
        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50) & (self.gediWREF[:,3] == 1.0))
        print('WREF useInd length',len(useIndW[0]))
        meansWREF=[]
        stdevsWREF=[]
        for i in covs:
            mean=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > i) & (self.alsWREF[useIndW,2] < (i+increment)))])
            meansWREF.append(mean)
            stdev=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > i) & (self.alsWREF[useIndW,2] < (i+increment)))])
            stdevsWREF.append(stdev)

        #GRSM
        #useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50))
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50) & (self.gediGRSM[:,3] == 1.0))
        print('GRSM useInd length',len(useIndG[0]))
        meansGRSM=[]
        stdevsGRSM=[]
        for i in covs:
            mean=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > i) & (self.alsGRSM[useIndG,2] < (i+increment)))])
            meansGRSM.append(mean)
            stdev=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > i) & (self.alsGRSM[useIndG,2] < (i+increment)))])
            stdevsGRSM.append(stdev)

        #LSBS
        #useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50))
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50) & (self.gediLSBS[:,3] == 1.0))
        print('LSBS useInd length',len(useIndL[0]))
        meansLSBS=[]
        stdevsLSBS=[]
        for i in covs:
            mean=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > i) & (self.alsLSBS[useIndL,2] < (i+increment)))])
            meansLSBS.append(mean)
            stdev=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > i) & (self.alsLSBS[useIndL,2] < (i+increment)))])
            stdevsLSBS.append(stdev)

        # Set the x-axis plot tick locations
        covers=np.arange(0.025,1.025,0.05)
        print('Cover Ticks:',covers)

        #plt.plot(self.alsWREF[:,2],self.gediWREF[:,4],'o',markersize=1)
        plt.plot(covers,meansWREF,label='WREF')
        plt.plot(covers,meansGRSM,label='GRSM')
        plt.plot(covers,meansLSBS,label='LSBS')
        plt.fill_between(covers,np.subtract(meansWREF,stdevsWREF),np.add(meansWREF,stdevsWREF),alpha=0.2)
        plt.fill_between(covers,np.subtract(meansGRSM,stdevsGRSM),np.add(meansGRSM,stdevsGRSM),alpha=0.2)
        plt.fill_between(covers,np.subtract(meansLSBS,stdevsLSBS),np.add(meansLSBS,stdevsLSBS),alpha=0.2)
        plt.title('Residual vs. Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,30])
        plt.xlabel('Canopy Cover (fraction)')
        plt.xlim([0,1])
        plt.legend(loc='upper left')
        #plt.savefig(args.outRoot+'coverResQual.png',dpi=300)
        #plt.close()
        #plt.clf()
        #print('Cover plot written to file '+args.outRoot+'coverResQual.png')
        plt.show()

#######################################

# The main block
if __name__ == '__main__':
    args=readCommands()

    plots=comparisonPlots()

    plots.readData()
    if args.hist:
        plots.plotHistograms()
    if args.box:
        plots.plotBoxplots()
