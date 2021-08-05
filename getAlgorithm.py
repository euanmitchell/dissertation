########################################
# Script to extract the number of times
# each algorithm setting is used from
# GEDI L2A files using.
########################################

import h5py
import numpy as np

########################################

file='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/processed_GEDI02_A_2019110110221_O01997_02_T03335_02_003_01_V002.h5'

beamList=['BEAM0000','BEAM0001','BEAM0010','BEAM0011','BEAM0110','BEAM0101','BEAM1000','BEAM1011']

print('Beam, Algorithm, PFT_Class, Quality_Flag, Surface_Flag, Solar_Elevation, Sensitivity')

f=h5py.File(file,'r')
for beam in beamList:
    if ((beam in list(f)) == False):
        #print(beam,'empty')
        continue
    #print(beam)
    beamID=np.array(f[beam]['beam'])
    algorithm=np.array(f[beam]['selected_algorithm'])
    sensitivity=np.array(f[beam]['sensitivity'])
    qualFlag=np.array(f[beam]['quality_flag'])
    pft=np.array(f[beam]['land_cover_data']['pft_class'])
    solar=np.array(f[beam]['solar_elevation'])
    surfaceFlag=np.array(f[beam]['surface_flag'])
    for i in range(algorithm.shape[0]):
        print(beamID[i], algorithm[i], pft[i], qualFlag[i], surfaceFlag[i], solar[i], sensitivity[i])

    '''counts=np.bincount(algorithm)
    print('The number of occurrences of algorithm 0 - 1 - 2 - 3 - 4 - 5 - 6 are')
    print(counts)'''
