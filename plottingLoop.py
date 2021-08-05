####################################
# Script to loop through directory
# and make a shell script calling
# 'plotComparison.py' for each file.
####################################

import os
import glob

###################################

root='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/area002/'
simDir='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/area002/sim07/' #72
metricDir='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/area002/metric07/'
subsetDir='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/area002/subset/'
l2aDir='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/'

outFile='shell/plotCompLaSelva_002_07.sh'
outRoot='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/area002/waveforms07/plot'
with open(outFile,'a') as outf:
    outf.write('#!/usr/bin/env bash\n\n')

    inputs='*.h5'
    path=os.path.join(simDir,inputs)
    simList=glob.glob(path)
    i=1

    for simFile in simList:
        #print(simFile)
        outf.write('echo "Processing file ' + str(i) + ' of ' + str(len(simList))+'"\n')
        if len(simFile) == 143:
            name=simFile.split('/')[-1]
            filename=name[7:]
            #print(filename)
            simID=name[:7]

        '''elif len(simFile) == 148:           # EDITED FOR ADDITION OF "newFiles" TO PATH
            name=simFile.split('/')[-1]
            filename=name[6:]
            simID=name[:6]'''

        subset=root+'subset/subset_'+filename
        sim=simFile
        metric=root+'metric07/'+simID+filename[:-2]+'metric.txt'
        list=filename.split('_')
        l2a=l2aDir+'processed_GEDI02_A_'+list[2]+'_'+list[3]+'_'+list[4]+'_'+list[5]+'_'+list[6]+'_003_'+list[8]+'_'+list[9]

        outf.write('python3 myPlotComparison.py --real {0} --sim {1} --metric {2} --l2a {3} --outRoot {4}\n\n'.format(subset,sim,metric,l2a,outRoot))
        i+=1
