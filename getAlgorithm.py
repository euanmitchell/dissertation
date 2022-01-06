########################################
# Script to extract the number of times
# each algorithm setting is used from
# GEDI L2A files using.
########################################

import h5py
import numpy as np
import os
import glob

########################################

def readFile():

    dir='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/grsm/l2a/'
    inputs='*.h5'
    path=os.path.join(dir,inputs)
    fileList=glob.glob(path)

    totalWaveforms=0
    totalPFT=np.empty(0,dtype='int64')

    for file in fileList:
        print('Reading file',file)

        beamList=['BEAM0000','BEAM0001','BEAM0010','BEAM0011','BEAM0110','BEAM0101','BEAM1000','BEAM1011']

        #print('Beam, Algorithm, PFT_Class, Quality_Flag, Surface_Flag, Solar_Elevation, Sensitivity')

        f=h5py.File(file,'r')
        fileWaveforms=0
        fileAlg=np.empty(0,dtype='int64')
        filePFT=np.empty(0,dtype='int64')
        for beam in beamList:
            if ((beam in list(f)) == False):
                #print(beam,'empty')
                continue
            #print(beam)
            try:
                beamID=np.array(f[beam]['beam'])
                algorithm=np.array(f[beam]['selected_algorithm'])
                sensitivity=np.array(f[beam]['sensitivity'])
                qualFlag=np.array(f[beam]['quality_flag'])
                pft=np.array(f[beam]['land_cover_data']['pft_class'])
                solar=np.array(f[beam]['solar_elevation'])
                surfaceFlag=np.array(f[beam]['surface_flag'])
                #for i in range(algorithm.shape[0]):
                    #print(beamID[i], algorithm[i], pft[i], qualFlag[i], surfaceFlag[i], solar[i], sensitivity[i])

                #counts=np.bincount(algorithm)
                #print('The number of occurrences of algorithm 0 - 1 - 2 - 3 - 4 - 5 - 6 are')
                #print(counts)
                fileWaveforms+=int(algorithm.shape[0])
                fileAlg=np.append(fileAlg,algorithm)
                filePFT=np.append(filePFT,pft)

            except:
                print('Selected_algorithm object does not exist in',file,beam)

        #counts=np.bincount(fileMaster)
        #print(counts)
        totalWaveforms+=fileWaveforms
        totalPFT=np.append(totalPFT,filePFT)

    print('The total number of waveforms is',totalWaveforms)
    useInd=np.where(totalPFT==9)
    print('The number of waveforms where PFT is 9 is',len(useInd[0]))
    #print(totalPFT[useInd])

#######################################

# The main block
if __name__ == '__main__':
    readFile()
