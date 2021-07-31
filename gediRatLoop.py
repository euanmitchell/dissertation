import numpy as np

minX=569000
maxX=582000
minY=5075000
maxY=5083000

step=1000

inList='~/Dissertation/data/groundBounds.neon_wref.txt'
outRoot='/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/grid/'
pulse='~/Dissertation/code/pulses/meanPulse.BEAM1000.filt'
outFile='shell/gediRatWref2.sh'

with open(outFile,'a') as outf:
    outf.write('#!/usr/bin/env bash\n\n')

counter=1
for x in np.arange(minX,maxX,step):
    for y in np.arange(minY,maxY,step):
        output=outRoot+str(x)+'.'+str(y)+'.h5'
        with open(outFile,'a') as outf:
            outf.write('echo "Processing subset ' + str(counter) + ' of 182"\n')
            outf.write('gediRat -inList {0} -output {1} -ground -hdf -gridBound {2} {3} {4} {5} -gridStep 20 -readPulse {6} -checkCover -pBuff 2 -countOnly\n\n'.format(inList,output,x,x+step,y,y+step,pulse))
        counter+=1

with open(outFile,'a') as outf:
    outf.write('echo "All done"')
