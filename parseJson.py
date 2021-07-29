####################################
# Attempt to cleanup and reduce size
# of the JSON files output from
# GEDI_Subsetter.py for quicker
# loading and rendering in QGIS
####################################

import json
import os
import glob

dir="/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/subset/2009/output/"
inputs='*.json'
path=os.path.join(dir,inputs)
files=glob.glob(path)

features=['ocean','sea_ice','land_ice','inland_water','noise_mean_corrected',
    'rx_sample_count','stale_return_flag','tx_sample_count','txwaveform','rxwaveform']

for file in files:
    filenames=file.split('/')
    filename=filenames[-1]
    outFile=dir+'parsed'+filename[6:]
    print('Cleaning file:', file)
    with open(outFile, 'w') as w:
        with open(file, 'r') as r:
            data=json.load(r)
            for line in data['features']:
                for item in features:
                    del line['properties'][item]
            json.dump(data,w,indent=6)
            print('Written:',outFile)
