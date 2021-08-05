'''
A script to plot a comparison
between real and simulated GEDI
data. Make sure to collocate first
2020
'''


##################################
import numpy as np
import h5py
from sys import exit
import matplotlib.pyplot as plt
if __name__ == '__main__':
  import argparse

# import a GEDI object
import sys,os
sys.path.append(os.environ["GEDIRAT_ROOT"])
from gediHandler import gediData


#######################################
# read the command line

if __name__ == '__main__':
  def readCommands():
    '''
    Read commandline arguments
    '''
    p = argparse.ArgumentParser(description=("Writes out properties of GEDI waveform files"))
    p.add_argument("--real",dest="realName",type=str,help=("Input GEDI HDF5 filename"))
    p.add_argument("--sim", dest ="simName", type=str,help=("Input simulated GEDI HDF5 filename"))
    p.add_argument("--outRoot",dest="outRoot",type=str,default='test',help=("Output graph filename root"))
    p.add_argument("--metric",dest="metName", type=str, default=" ", help=("Optional: Metric filename, for ground estimates"))
    p.add_argument("--l2a",dest="l2aName", type=str,default=" ",help=("Optional: GEDI L2A filename"))
    p.add_argument("--datumDiff",dest="datum", type=float, default=0.0, help=("Optional: Datum difference"))
    cmdargs = p.parse_args()
    return cmdargs


###########################################

class gediMetrics():
  '''Class to hold GEDI metrics'''
  def __init__(self,filename):
    '''Initialiser'''
    if(filename!=" "):
      self.waveID=np.genfromtxt(filename,usecols=(0,), unpack=True, dtype=str,comments='#')
      # Update self.ZG to be the ALS ground estimate
      self.ZG,self.sense=np.loadtxt(filename,usecols=(1,112), unpack=True, dtype=float,comments='#')
      self.shotN=np.empty(self.waveID.shape,dtype=int)
      for i in range(0,self.waveID.shape[0]):
        self.shotN[i]=int(self.waveID[i].split('.')[2])
      self.useMet=1   # use it
    else:
      self.useMet=0   # file is empty. Don't use

###########################################

class gediReal():
    '''My new class to hold GEDI L2A metrics'''
    def __init__(self,filename):
        '''The initialiser'''
        if(filename!=" "):
            f=h5py.File(filename,'r')
            self.shotN=np.empty(0,dtype=np.uint64)
            self.groundZ=np.empty(0,dtype=float)
            self.lat=np.empty(0)
            beamlist=['BEAM0000','BEAM0001','BEAM0010','BEAM0011','BEAM0110','BEAM0101','BEAM1000','BEAM1011']
            for beam in beamlist:
                try:
                    self.shotN=np.append(self.shotN,f[beam]['shot_number'])
                    self.groundZ=np.append(self.groundZ,f[beam]['elev_lowestmode'])
                    self.beamID=np.append(self.beamID,f[beam]['beam'])
                    self.lat=np.append(self.lat,f[beam]['lat_lowestmode'])
                    print('Data in',beam)
                except:
                    print('Empty beam',beam)
            self.useReal=1
        else:
            self.useReal=0


###########################################

def plotComparison(gedi,sim,met,real,outRoot,datum=0):
  '''Function to plot graph comparisons'''

  # loop over waveforms
  for i in range(0,gedi.nWaves):        #gedi.nWaves
    # find corresponding metric
    shotN=gedi.waveID[i]
    beam=gedi.waveID[i]
    simInd=np.where((sim.shotN==shotN)) #*(sim.beamID==gedi.beamID[i]))

    # do we have any?
    if(len(simInd)>0):
      if(len(simInd[0])>0):
        simInd=simInd[0][0]
      else:
        continue
    else:
      continue

    if(met.useMet==1):
      metInd=np.where(met.shotN==shotN)
      if(len(metInd)>0):
        metInd=metInd[0]
      else:
        continue

    if(real.useReal==1):
        realInd=np.where(real.shotN==shotN)
        if(len(realInd)>0):
            realInd=realInd[0]
        else:
            continue

    # make a graph
    output=outRoot+".beam."+str(shotN)+".png"

    # set z arrays
    gedi.nBins=gedi.lenInds[i]
    gedi.res=(gedi.Z0[i]-gedi.ZN[i])/gedi.nBins
    gedi.setOneZ(i)
    sim.setOneZ(simInd)

    # set scalars and bounds
    reflScale,meanN,stdev=gedi.meanNoise(i)
    minX,maxX=gedi.findBounds(meanN,stdev,i)
    #temp=np.sort(np.copy(gedi.wave[i][0:gedi.nBins]))
    #minY=temp[int(temp.shape[0]*0.01)]-5  # to avoid artefacts
    maxY=1.2*np.max(gedi.wave[i][0:gedi.nBins])
    minY=meanN-0.1*(maxY-meanN)

    # arrange the plot
    plt.xlim(minY,maxY)
    plt.ylim((minX+datum,maxX+datum))
    plt.ylabel('Elevation (m)')
    plt.xlabel('DN')
    plt.plot(gedi.wave[i][0:gedi.nBins],gedi.z+datum,label='GEDI')
    plt.plot(sim.wave[simInd][0:sim.nBins]*reflScale+meanN,sim.z,label='Sim')
    plt.plot(sim.gWave[simInd][0:sim.nBins]*reflScale+meanN,sim.z,label='ALS ground')
    if(met.useMet==1):
    # This is now the ALS 'true ground' estimate
      plt.plot([minY,maxY], [met.ZG[metInd],met.ZG[metInd]], color='g', linestyle='--', linewidth=2, label='ALS zGround')
    if(real.useReal==1):
    # This should be the L2A ground estimate
        plt.plot([minY,maxY], [real.groundZ[realInd],real.groundZ[realInd]], color='b', linestyle='--', linewidth=2, label='GEDI zGround')
    plt.legend()
    plt.savefig(output)
    plt.close()
    plt.clf()
    print("Written to",output)

  return

###########################################
# the main block

if __name__ == '__main__':
  # read the command line
  cmd=readCommands()

  # read real GEDI data
  gedi=gediData(filename=cmd.realName)

  # read simulated GEDI data
  sim=gediData(filename=cmd.simName)

  # read metrics, if needed
  met=gediMetrics(cmd.metName)

  # read the L2A file
  real=gediReal(cmd.l2aName)

  # plot them up
  plotComparison(gedi,sim,met,real,cmd.outRoot,datum=cmd.datum)
