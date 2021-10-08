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

        #useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50))
        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50) & (self.gediWREF[:,3] == 1.0))
        print('WREF useInd length',len(useIndW[0]))

        #WREF
        m1W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] < 5))])
        m2W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 5) & (self.alsWREF[useIndW,1] < 10))])
        m3W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 10) & (self.alsWREF[useIndW,1] < 15))])
        m4W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 15) & (self.alsWREF[useIndW,1] < 20))])
        m5W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 20) & (self.alsWREF[useIndW,1] < 25))])
        m6W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 25) & (self.alsWREF[useIndW,1] < 30))])
        m7W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 30) & (self.alsWREF[useIndW,1] < 35))])
        m8W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 35) & (self.alsWREF[useIndW,1] < 40))])
        m9W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 40) & (self.alsWREF[useIndW,1] < 45))])
        m10W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 45) & (self.alsWREF[useIndW,1] < 50))])
        m11W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 50))])

        s1W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] < 5))])
        s2W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 5) & (self.alsWREF[useIndW,1] < 10))])
        s3W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 10) & (self.alsWREF[useIndW,1] < 15))])
        s4W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 15) & (self.alsWREF[useIndW,1] < 20))])
        s5W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 20) & (self.alsWREF[useIndW,1] < 25))])
        s6W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 25) & (self.alsWREF[useIndW,1] < 30))])
        s7W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 30) & (self.alsWREF[useIndW,1] < 35))])
        s8W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 35) & (self.alsWREF[useIndW,1] < 40))])
        s9W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 40) & (self.alsWREF[useIndW,1] < 45))])
        s10W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 45) & (self.alsWREF[useIndW,1] < 50))])
        s11W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,1] > 50))])

        #GRSM
        #useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50))
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50) & (self.gediGRSM[:,3] == 1.0))
        print('GRSM useInd length',len(useIndG[0]))

        m1G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] < 5))])
        m2G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 5) & (self.alsGRSM[useIndG,1] < 10))])
        m3G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 10) & (self.alsGRSM[useIndG,1] < 15))])
        m4G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 15) & (self.alsGRSM[useIndG,1] < 20))])
        m5G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 20) & (self.alsGRSM[useIndG,1] < 25))])
        m6G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 25) & (self.alsGRSM[useIndG,1] < 30))])
        m7G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 30) & (self.alsGRSM[useIndG,1] < 35))])
        m8G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 35) & (self.alsGRSM[useIndG,1] < 40))])
        m9G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 40) & (self.alsGRSM[useIndG,1] < 45))])
        m10G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 45) & (self.alsGRSM[useIndG,1] < 50))])
        m11G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 50))])

        s1G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] < 5))])
        s2G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 5) & (self.alsGRSM[useIndG,1] < 10))])
        s3G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 10) & (self.alsGRSM[useIndG,1] < 15))])
        s4G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 15) & (self.alsGRSM[useIndG,1] < 20))])
        s5G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 20) & (self.alsGRSM[useIndG,1] < 25))])
        s6G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 25) & (self.alsGRSM[useIndG,1] < 30))])
        s7G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 30) & (self.alsGRSM[useIndG,1] < 35))])
        s8G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 35) & (self.alsGRSM[useIndG,1] < 40))])
        s9G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 40) & (self.alsGRSM[useIndG,1] < 45))])
        s10G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 45) & (self.alsGRSM[useIndG,1] < 50))])
        s11G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,1] > 50))])

        #LSBS
        #useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50))
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50) & (self.gediLSBS[:,3] == 1.0))
        print('LSBS useInd length',len(useIndL[0]))

        m1L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] < 5))])
        m2L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 5) & (self.alsLSBS[useIndL,1] < 10))])
        m3L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 10) & (self.alsLSBS[useIndL,1] < 15))])
        m4L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 15) & (self.alsLSBS[useIndL,1] < 20))])
        m5L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 20) & (self.alsLSBS[useIndL,1] < 25))])
        m6L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 25) & (self.alsLSBS[useIndL,1] < 30))])
        m7L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 30) & (self.alsLSBS[useIndL,1] < 35))])
        m8L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 35) & (self.alsLSBS[useIndL,1] < 40))])
        m9L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 40) & (self.alsLSBS[useIndL,1] < 45))])
        m10L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 45) & (self.alsLSBS[useIndL,1] < 50))])
        m11L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 50))])

        s1L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] < 5))])
        s2L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 5) & (self.alsLSBS[useIndL,1] < 10))])
        s3L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 10) & (self.alsLSBS[useIndL,1] < 15))])
        s4L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 15) & (self.alsLSBS[useIndL,1] < 20))])
        s5L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 20) & (self.alsLSBS[useIndL,1] < 25))])
        s6L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 25) & (self.alsLSBS[useIndL,1] < 30))])
        s7L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 30) & (self.alsLSBS[useIndL,1] < 35))])
        s8L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 35) & (self.alsLSBS[useIndL,1] < 40))])
        s9L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 40) & (self.alsLSBS[useIndL,1] < 45))])
        s10L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 45) & (self.alsLSBS[useIndL,1] < 50))])
        s11L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,1] > 50))])

        slopes=np.array([2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5])
        meansWREF=[m1W,m2W,m3W,m4W,m5W,m6W,m7W,m8W,m9W,m10W,m11W]
        stdevsWREF=[s1W,s2W,s3W,s4W,s5W,s6W,s7W,s8W,s9W,s10W,s11W]
        meansGRSM=[m1G,m2G,m3G,m4G,m5G,m6G,m7G,m8G,m9G,m10G,m11G]
        stdevsGRSM=[s1G,s2G,s3G,s4G,s5G,s6G,s7G,s8G,s9G,s10G,s11G]
        meansLSBS=[m1L,m2L,m3L,m4L,m5L,m6L,m7L,m8L,m9L,m10L,m11L]
        stdevsLSBS=[s1L,s2L,s3L,s4L,s5L,s6L,s7L,s8L,s9L,s10L,s11L]

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
        plt.xlim([0,55])
        plt.legend(loc='upper left')
        #plt.savefig(args.outRoot+'slopeResQual.png',dpi=300)
        #plt.close()
        #plt.clf()
        #print('Slope plot written to file '+args.outRoot+'slopeResQual.png')
        plt.show()

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Cover Plots

        covs=np.arange(0,1.0,0.1)

        #WREF
        #useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50))
        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50) & (self.gediWREF[:,3] == 1.0))
        print('WREF useInd length',len(useIndW[0]))

        meansWREF=[]
        for i in covs:
            mean=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > i) & (self.alsWREF[useIndW,2] < (i+0.1)))])
            meansWREF.append(mean)
            print(i,i+0.1)
        stdevsWREF=[]
        for i in covs:
            stdev=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > i) & (self.alsWREF[useIndW,2] < (i+0.1)))])
            stdevsWREF.append(stdev)

        '''m1W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] < 0.1))])
        m2W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.1) & (self.alsWREF[useIndW,2] < 0.2))])
        m3W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.2) & (self.alsWREF[useIndW,2] < 0.3))])
        m4W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.3) & (self.alsWREF[useIndW,2] < 0.4))])
        m5W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.4) & (self.alsWREF[useIndW,2] < 0.5))])
        m6W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.5) & (self.alsWREF[useIndW,2] < 0.6))])
        m7W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.6) & (self.alsWREF[useIndW,2] < 0.7))])
        m8W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.7) & (self.alsWREF[useIndW,2] < 0.8))])
        m9W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.8) & (self.alsWREF[useIndW,2] < 0.9))])
        m10W=np.mean(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.9))])'''

        s1W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] < 0.1))])
        s2W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.1) & (self.alsWREF[useIndW,2] < 0.2))])
        s3W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.2) & (self.alsWREF[useIndW,2] < 0.3))])
        s4W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.3) & (self.alsWREF[useIndW,2] < 0.4))])
        s5W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.4) & (self.alsWREF[useIndW,2] < 0.5))])
        s6W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.5) & (self.alsWREF[useIndW,2] < 0.6))])
        s7W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.6) & (self.alsWREF[useIndW,2] < 0.7))])
        s8W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.7) & (self.alsWREF[useIndW,2] < 0.8))])
        s9W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.8) & (self.alsWREF[useIndW,2] < 0.9))])
        s10W=np.std(self.gediWREF[useIndW,4][np.where((self.alsWREF[useIndW,2] > 0.9))])

        #GRSM
        #useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50))
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50) & (self.gediGRSM[:,3] == 1.0))
        print('GRSM useInd length',len(useIndG[0]))

        m1G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] < 0.1))])
        m2G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.1) & (self.alsGRSM[useIndG,2] < 0.2))])
        m3G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.2) & (self.alsGRSM[useIndG,2] < 0.3))])
        m4G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.3) & (self.alsGRSM[useIndG,2] < 0.4))])
        m5G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.4) & (self.alsGRSM[useIndG,2] < 0.5))])
        m6G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.5) & (self.alsGRSM[useIndG,2] < 0.6))])
        m7G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.6) & (self.alsGRSM[useIndG,2] < 0.7))])
        m8G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.7) & (self.alsGRSM[useIndG,2] < 0.8))])
        m9G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.8) & (self.alsGRSM[useIndG,2] < 0.9))])
        m10G=np.mean(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.9))])

        s1G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] < 0.1))])
        s2G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.1) & (self.alsGRSM[useIndG,2] < 0.2))])
        s3G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.2) & (self.alsGRSM[useIndG,2] < 0.3))])
        s4G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.3) & (self.alsGRSM[useIndG,2] < 0.4))])
        s5G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.4) & (self.alsGRSM[useIndG,2] < 0.5))])
        s6G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.5) & (self.alsGRSM[useIndG,2] < 0.6))])
        s7G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.6) & (self.alsGRSM[useIndG,2] < 0.7))])
        s8G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.7) & (self.alsGRSM[useIndG,2] < 0.8))])
        s9G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.8) & (self.alsGRSM[useIndG,2] < 0.9))])
        s10G=np.std(self.gediGRSM[useIndG,4][np.where((self.alsGRSM[useIndG,2] > 0.9))])

        #LSBS
        #useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50))
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50) & (self.gediLSBS[:,3] == 1.0))
        print('LSBS useInd length',len(useIndL[0]))

        m1L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] < 0.1))])
        m2L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.1) & (self.alsLSBS[useIndL,2] < 0.2))])
        m3L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.2) & (self.alsLSBS[useIndL,2] < 0.3))])
        m4L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.3) & (self.alsLSBS[useIndL,2] < 0.4))])
        m5L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.4) & (self.alsLSBS[useIndL,2] < 0.5))])
        m6L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.5) & (self.alsLSBS[useIndL,2] < 0.6))])
        m7L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.6) & (self.alsLSBS[useIndL,2] < 0.7))])
        m8L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.7) & (self.alsLSBS[useIndL,2] < 0.8))])
        m9L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.8) & (self.alsLSBS[useIndL,2] < 0.9))])
        m10L=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.9))])

        s1L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] < 0.1))])
        s2L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.1) & (self.alsLSBS[useIndL,2] < 0.2))])
        s3L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.2) & (self.alsLSBS[useIndL,2] < 0.3))])
        s4L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.3) & (self.alsLSBS[useIndL,2] < 0.4))])
        s5L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.4) & (self.alsLSBS[useIndL,2] < 0.5))])
        s6L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.5) & (self.alsLSBS[useIndL,2] < 0.6))])
        s7L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.6) & (self.alsLSBS[useIndL,2] < 0.7))])
        s8L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.7) & (self.alsLSBS[useIndL,2] < 0.8))])
        s9L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.8) & (self.alsLSBS[useIndL,2] < 0.9))])
        s10L=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > 0.9))])

        covers=np.array([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95])
        #meansWREF=[m1W,m2W,m3W,m4W,m5W,m6W,m7W,m8W,m9W,m10W]
        #stdevsWREF=[s1W,s2W,s3W,s4W,s5W,s6W,s7W,s8W,s9W,s10W]
        meansGRSM=[m1G,m2G,m3G,m4G,m5G,m6G,m7G,m8G,m9G,m10G]
        stdevsGRSM=[s1G,s2G,s3G,s4G,s5G,s6G,s7G,s8G,s9G,s10G]
        meansLSBS=[m1L,m2L,m3L,m4L,m5L,m6L,m7L,m8L,m9L,m10L]
        stdevsLSBS=[s1L,s2L,s3L,s4L,s5L,s6L,s7L,s8L,s9L,s10L]

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

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Moving average solution do this for covers?
        '''def movingAvg(x,n):
            cumsum=np.cumsum(np.insert(x,0,0))
            return (cumsum[n:] - cumsum[:-n])/float(n)

        wrefSorted=self.alsWREF[self.alsWREF[:,2].argsort()]
        n=50
        WREF=movingAvg(wrefSorted[:,4],n)
        plt.plot(wrefSorted[49:,2],WREF)
        plt.show()'''

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
