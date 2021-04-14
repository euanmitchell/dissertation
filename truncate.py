#######################################
# Truncate full GEDI granule URLs to
# filename alone for pasting into
# NASA EarthData Search page
#######################################

import csv
import argparse

def readCommands():
    p=argparse.ArgumentParser(description=('Input and output files for truncating'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--inFile', dest='inFile', type=str, default='',
        help=('The input URL.\nDefault is not set'))
    p.add_argument('--outFile', dest='outFile', type=str, default='',
        help=('The output text file.\nDefault is not set'))
    cmdargs=p.parse_args()
    return cmdargs

args=readCommands()

inFile=args.inFile
outFile=args.outFile

filenames = ''

with open (inFile) as data:
    dataReader = csv.reader(data)
    for row in dataReader:
        row=row[0].split('/')
        filename=row[-1]
        filenames += filename + ','

with open (outFile, 'w') as cleanedData:
    cleanedData.write(filenames)
    cleanedData.close()
