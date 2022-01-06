#######################################
# Script to read text files output by
# alsComp.py and read them in to make
# histogram and box plots for all
# three sites.
#######################################

import h5py
import numpy as np
from mpl_toolkits import mplot3d
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
    p.add_argument('--stats', dest='stats', action='store_true', default=False,
        help=('Print statistics for the three sites'))
    p.add_argument('--hist', dest='hist', action='store_true', default=False,
        help=('Make histogram plots'))
    p.add_argument('--box', dest='box', action='store_true', default=False,
        help=('Make box plots'))
    cmdargs=p.parse_args()
    return cmdargs

#######################################

class comparisonPlots(object):
    '''A class to handle reading and plotting ALS and GEDI ground estimates for multiple sites'''

    def __init__(self):
        '''Initialise the class'''

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

        if args.stats:
            print('Mean canopy height at WREF is',np.mean(self.alsWREF[useIndW,3]))
            print('Mean canopy height at GRSM is',np.mean(self.alsGRSM[useIndG,3]))
            print('Mean canopy height at LSBS is',np.mean(self.alsLSBS[useIndL,3]))

            print('Median canopy height at WREF is',np.median(self.alsWREF[useIndW,3]))
            print('Median canopy height at GRSM is',np.median(self.alsGRSM[useIndG,3]))
            print('Median canopy height at LSBS is',np.median(self.alsLSBS[useIndL,3]))

            print('Median canopy cover at WREF is',np.median(self.alsWREF[useIndW,2]))
            print('Median canopy cover at GRSM is',np.median(self.alsGRSM[useIndG,2]))
            print('Median canopy cover at LSBS is',np.median(self.alsLSBS[useIndL,2]))

            print('Mean slope at WREF is',np.mean(self.alsWREF[useIndW,1]))
            print('Mean slope at GRSM is',np.mean(self.alsGRSM[useIndG,1]))
            print('Mean slope at LSBS is',np.mean(self.alsLSBS[useIndL,1]))

            print('Median slope at WREF is',np.median(self.alsWREF[useIndW,1]))
            print('Median slope at GRSM is',np.median(self.alsGRSM[useIndG,1]))
            print('Median slope at LSBS is',np.median(self.alsLSBS[useIndL,1]))

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
            edgecolor='white',
            density=True
        )
        plt.title('Slope Distribution')
        plt.ylabel('Normalised Frequency')
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
            edgecolor='white',
            density=True
        )
        plt.title('Canopy Cover Distribution')
        plt.ylabel('Normalised Frequency')
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
            edgecolor='white',
            density=True
        )
        plt.title('Canopy Height Distribution')
        plt.ylabel('Normalised Frequency')
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
            edgecolor='white',
            density=True
        )
        plt.title('Beam Sensitivity Distribution')
        plt.ylabel('Normalised Frequency')
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

        # Slope insert panel histogram
        plt.hist(
            [self.alsWREF[:,1],self.alsGRSM[:,1],self.alsLSBS[:,1]],
            bins=np.arange(0,60+5,5),
            label=['WREF','GRSM','LSBS'],
            edgecolor='white',
            density=True
        )
        #plt.title('Beam Sens. - Cover Distribution')
        #plt.ylabel('Normalised Frequency')
        plt.xlabel('Slope (degrees)')
        plt.ylim([0,0.2])
        plt.xlim([0,60])
        plt.legend(loc='upper right')
        #plt.savefig(args.outRoot+'slopeHistAll.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

        #plt.plot(self.alsWREF[:,1],self.gediWREF[:,4],'o',markersize=1)
        plt.axhline(y=0.0,color='r')
        plt.plot(slopes,meansWREF,label='WREF')
        plt.plot(slopes,meansGRSM,label='GRSM')
        plt.plot(slopes,meansLSBS,label='LSBS')
        plt.fill_between(slopes,np.subtract(meansWREF,stdevsWREF),np.add(meansWREF,stdevsWREF),alpha=0.2)
        plt.fill_between(slopes,np.subtract(meansGRSM,stdevsGRSM),np.add(meansGRSM,stdevsGRSM),alpha=0.2)
        plt.fill_between(slopes,np.subtract(meansLSBS,stdevsLSBS),np.add(meansLSBS,stdevsLSBS),alpha=0.2)
        #plt.title('Residual vs. Slope')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,30])
        plt.yticks([-10, 0, 10, 20, 30])
        #plt.xlabel('Slope (degrees)')
        plt.xlim([0,60])
        plt.xticks(color='w')
        plt.legend(loc='upper left')
        #plt.savefig(args.outRoot+'slopeResQual.png',dpi=300)
        #plt.close()
        #plt.clf()
        #print('Slope plot written to file '+args.outRoot+'slopeResQual.png')
        plt.show()

        ########################################################################

        # Cover Plots

        # Set the number of data bins
        covs=np.arange(0,1.0,0.02)
        increment=0.02
        print('CoverBins:',covs)

        #WREF
        #useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50))
        # Add beam sensitivity threshold for plotting ...
        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50) & (self.gediWREF[:,3] == 1.0) & (self.gediWREF[:,1] < 0.98))
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
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50) & (self.gediGRSM[:,3] == 1.0) & (self.gediGRSM[:,1] < 0.98))
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
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50) & (self.gediLSBS[:,3] == 1.0) & (self.gediLSBS[:,1] < 0.98))
        print('LSBS useInd length',len(useIndL[0]))
        meansLSBS=[]
        stdevsLSBS=[]
        for i in covs:
            mean=np.mean(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > i) & (self.alsLSBS[useIndL,2] < (i+increment)))])
            meansLSBS.append(mean)
            stdev=np.std(self.gediLSBS[useIndL,4][np.where((self.alsLSBS[useIndL,2] > i) & (self.alsLSBS[useIndL,2] < (i+increment)))])
            stdevsLSBS.append(stdev)

        # Set the x-axis plot tick locations
        covers=np.arange(0.01,1.01,0.02)
        print('Cover Ticks:',covers)

        #plt.plot(self.alsWREF[:,2],self.gediWREF[:,4],'o',markersize=1)
        plt.plot(covers,meansWREF,label='WREF')
        plt.plot(covers,meansGRSM,label='GRSM')
        plt.plot(covers,meansLSBS,label='LSBS')
        plt.fill_between(covers,np.subtract(meansWREF,stdevsWREF),np.add(meansWREF,stdevsWREF),alpha=0.2)
        plt.fill_between(covers,np.subtract(meansGRSM,stdevsGRSM),np.add(meansGRSM,stdevsGRSM),alpha=0.2)
        plt.fill_between(covers,np.subtract(meansLSBS,stdevsLSBS),np.add(meansLSBS,stdevsLSBS),alpha=0.2)
        plt.axhline(y=0.0,color='r')
        plt.title('Residual vs. Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,30])
        plt.xlabel('Canopy Cover (fraction)')
        plt.xlim([0,1])
        plt.legend(loc='upper left')
        #plt.savefig(args.outRoot+'coverResBSlt98bin2.png',dpi=300)
        #plt.close()
        #plt.clf()
        #print('Cover plot written to file '+args.outRoot+'coverResBSlt98bin2.png')
        plt.show()

        ########################################################################

        # Beam Sens Plot

        # Set the number of data bins
        bins=np.arange(-0.05,1.0,0.02)
        increment=0.02
        print('Bins:',bins)

        #WREF
        useIndW=np.where((self.alsWREF[:,0] > 0) & (self.gediWREF[:,4] < 50) & (self.gediWREF[:,3] == 1.0))
        sensMinusCovW=np.subtract(self.gediWREF[useIndW,1],self.alsWREF[useIndW,2])
        print('WREF useInd length',len(useIndW[0]))
        meansWREF=[]
        stdevsWREF=[]
        for i in bins:
            mean=np.mean(self.gediWREF[useIndW,4][np.where((sensMinusCovW > i) & (sensMinusCovW < (i+increment)))])
            meansWREF.append(mean)
            stdev=np.std(self.gediWREF[useIndW,4][np.where((sensMinusCovW > i) & (sensMinusCovW < (i+increment)))])
            stdevsWREF.append(stdev)
        print(np.min(sensMinusCovW))
        #plt.plot(sensMinusCovW,self.alsWREF[useIndW,4],'o',markersize=2)
        #plt.title('Wind River Raw')
        #plt.show()

        x1=self.gediWREF[useIndW,4][np.where((sensMinusCovW < 0.0))]
        x2=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.0) & (sensMinusCovW < 0.1))]
        x3=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.1) & (sensMinusCovW < 0.2))]
        x4=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.2) & (sensMinusCovW < 0.3))]
        x5=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.3) & (sensMinusCovW < 0.4))]
        x6=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.4) & (sensMinusCovW < 0.5))]
        x7=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.5) & (sensMinusCovW < 0.6))]
        x8=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.6) & (sensMinusCovW < 0.7))]
        x9=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.7) & (sensMinusCovW < 0.8))]
        x10=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.8) & (sensMinusCovW < 0.9))]
        x11=self.gediWREF[useIndW,4][np.where((sensMinusCovW > 0.9) & (sensMinusCovW < 1.0))]
        data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11]

        plt.boxplot(data,sym='')
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('Wind River Experimental Forest')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,40])
        plt.xlabel('Beam Sensitivity - Canopy Cover')
        plt.xticks(np.arange(1,12,step=1),['-0.05', '0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
        #plt.savefig(args.outRoot+'wrefBox.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #GRSM
        useIndG=np.where((self.alsGRSM[:,0] > 0) & (self.gediGRSM[:,4] < 50) & (self.gediGRSM[:,3] == 1.0))
        sensMinusCovG=np.subtract(self.gediGRSM[useIndG,1],self.alsGRSM[useIndG,2])
        print('GRSM useInd length',len(useIndG[0]))
        meansGRSM=[]
        stdevsGRSM=[]
        for i in bins:
            mean=np.mean(self.gediGRSM[useIndG,4][np.where((sensMinusCovG > i) & (sensMinusCovG < (i+increment)))])
            meansGRSM.append(mean)
            stdev=np.std(self.gediGRSM[useIndG,4][np.where((sensMinusCovG > i) & (sensMinusCovG < (i+increment)))])
            stdevsGRSM.append(stdev)
        print(np.min(sensMinusCovG))
        #plt.plot(sensMinusCovG,self.alsGRSM[useIndG,4],'o',markersize=2)
        #plt.title('Smokies Raw')
        #plt.show()

        x1=self.gediGRSM[useIndG,4][np.where((sensMinusCovG < 0.0))]
        x2=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.0) & (sensMinusCovG < 0.1))]
        x3=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.1) & (sensMinusCovG < 0.2))]
        x4=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.2) & (sensMinusCovG < 0.3))]
        x5=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.3) & (sensMinusCovG < 0.4))]
        x6=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.4) & (sensMinusCovG < 0.5))]
        x7=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.5) & (sensMinusCovG < 0.6))]
        x8=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.6) & (sensMinusCovG < 0.7))]
        x9=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.7) & (sensMinusCovG < 0.8))]
        x10=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.8) & (sensMinusCovG < 0.9))]
        x11=self.gediGRSM[useIndG,4][np.where((sensMinusCovG > 0.9) & (sensMinusCovG < 1.0))]
        data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11]

        plt.boxplot(data,sym='')
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('Great Smoky Mountains NP')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,40])
        plt.xlabel('Beam Sensitivity - Canopy Cover')
        plt.xticks(np.arange(1,12,step=1),['-0.05', '0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
        #plt.savefig(args.outRoot+'grsmBox.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #LSBS
        useIndL=np.where((self.alsLSBS[:,0] > 0) & (self.gediLSBS[:,4] < 50) & (self.gediLSBS[:,3] == 1.0))
        sensMinusCovL=np.subtract(self.gediLSBS[useIndL,1],self.alsLSBS[useIndL,2])
        print('LSBS useInd length',len(useIndL[0]))
        meansLSBS=[]
        stdevsLSBS=[]
        for i in bins:
            mean=np.mean(self.gediLSBS[useIndL,4][np.where((sensMinusCovL > i) & (sensMinusCovL < (i+increment)))])
            meansLSBS.append(mean)
            stdev=np.std(self.gediLSBS[useIndL,4][np.where((sensMinusCovL > i) & (sensMinusCovL < (i+increment)))])
            stdevsLSBS.append(stdev)
        print(np.min(sensMinusCovL))
        #plt.plot(sensMinusCovL,self.alsLSBS[useIndL,4],'o',markersize=2)
        #plt.title('La Selva Raw')
        #plt.show()

        x1=self.gediLSBS[useIndL,4][np.where((sensMinusCovL < 0.0))]
        x2=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.0) & (sensMinusCovL < 0.1))]
        x3=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.1) & (sensMinusCovL < 0.2))]
        x4=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.2) & (sensMinusCovL < 0.3))]
        x5=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.3) & (sensMinusCovL < 0.4))]
        x6=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.4) & (sensMinusCovL < 0.5))]
        x7=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.5) & (sensMinusCovL < 0.6))]
        x8=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.6) & (sensMinusCovL < 0.7))]
        x9=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.7) & (sensMinusCovL < 0.8))]
        x10=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.8) & (sensMinusCovL < 0.9))]
        x11=self.gediLSBS[useIndL,4][np.where((sensMinusCovL > 0.9) & (sensMinusCovL < 1.0))]
        data=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11]

        plt.boxplot(data,sym='')
        plt.axhline(y=0.0,color='r',linestyle='--')
        plt.title('La Selva Biological Station')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-20,40])
        plt.xlabel('Beam Sensitivity - Canopy Cover')
        plt.xticks(np.arange(1,12,step=1),['-0.05', '0.05', '0.15', '0.25', '0.35', '0.45', '0.55', '0.65', '0.75', '0.85', '0.95'])
        #plt.savefig(args.outRoot+'lsbsBox.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Histogram of beam sensitivity - canopy cover parameter
        plt.hist(
            [sensMinusCovW[0],sensMinusCovG[0],sensMinusCovL[0]],
            bins=np.arange(-0.04,1.0+0.02,0.02),
            label=['WREF','GRSM','LSBS'],
            edgecolor='white',
            density=True
        )
        #plt.title('Beam Sens. - Cover Distribution')
        #plt.ylabel('Normalised Frequency')
        plt.xlabel('Beam Sensitivity - Canopy Cover')
        plt.ylim([0,28])
        plt.xlim([-0.06,1])
        plt.legend(loc='upper right')
        #plt.savefig(args.outRoot+'slopeHistAll.png',dpi=300)
        #plt.close()
        #plt.clf()
        plt.show()

        # Set the x-axis plot tick locations
        ticks=np.arange(-0.04,1.01,0.02)
        print('X Ticks:',ticks)

        #plt.plot(self.alsWREF[:,2],self.gediWREF[:,4],'o',markersize=1)
        plt.plot(ticks,meansWREF,label='WREF')
        plt.plot(ticks,meansGRSM,label='GRSM')
        plt.plot(ticks,meansLSBS,label='LSBS')
        plt.fill_between(ticks,np.subtract(meansWREF,stdevsWREF),np.add(meansWREF,stdevsWREF),alpha=0.2)
        plt.fill_between(ticks,np.subtract(meansGRSM,stdevsGRSM),np.add(meansGRSM,stdevsGRSM),alpha=0.2)
        plt.fill_between(ticks,np.subtract(meansLSBS,stdevsLSBS),np.add(meansLSBS,stdevsLSBS),alpha=0.2)
        plt.axhline(y=0.0,color='r')
        #plt.title('Residual vs. Beam Sensitivity minus Canopy Cover')
        plt.ylabel('Ground Residual (m)')
        plt.ylim([-15,30])
        plt.yticks([-10, 0, 10, 20, 30])
        #plt.xlabel('Beam Sensitivity - Canopy Cover')
        plt.xlim([-0.06,1])
        plt.xticks(color='w')
        plt.legend(loc='upper right')
        #plt.savefig(args.outRoot+'resSensMinusCover.png',dpi=300)
        #plt.close()
        #plt.clf()
        #print('Cover plot written to file '+args.outRoot+'resSensMinusCover.png')
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
