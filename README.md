# Dissertation

Code developed during work on my MSc dissertation, spring and summer 2021.

## Files in this repository

#### alsComp.py
* First Python script to read text file output from ```gediMetric``` and GEDI L2A HDF5 files and generate ground residual comparison plots.

#### coords.sh
* A simple shell script to execute ```gediHandler.py``` with ```--writeCoords``` option to write footprint coordinates to text files.

#### gridComp.py
* Python script to read text file output from running ```gediMetric``` over simulated waveforms and extract data from TIFFs for comparison plots.

#### subset.sh
* A simple shell script to loop through a directory and execute the ```subsetGEDI.py``` program and spatially subset full orbit GEDI L1B files.

#### truncate.py
* Short Python script to strip the filenames from the end of URLs obtained from the GEDI finder site and write to a text file
as a comma separated list that can be pasted into the NASA EarthData Search page granule IDs search box.
