# Dissertation

Code developed during work on my MSc dissertation, spring and summer 2021.

## Files in this repository

### In the ```shell``` folder:

#### coords.sh
* A simple shell script to execute ```gediHandler.py``` with ```--writeCoords``` option to write footprint coordinates to text files.

#### metric.sh
* A simple shell script to loop through a directory and execute the ```gediMetric``` program for each output HDF5 file from ```collocateWaves```.

#### subset.sh
* A simple shell script to loop through a directory and execute the ```subsetGEDI.py``` program and spatially subset full orbit GEDI L1B files.

### In the main folder:

#### alsComp.py
* Main Python script to read text file output from ```gediMetric``` and GEDI L2A HDF5 files and generate ground residual comparison plots and statistics. The command line options available are:
  * --als - the path to a directory containing text file output from ```gediMetric```.
  * --alsFile - the path and filename for a single ```gediMetric``` text file.
  * --gedi - the path to a directory containing GEDI L2A data files.
  * --gediFile - the path and filename for a single GEDI L2A data file.
  * --beams - specify whether to use all GEDI beams (the default), or only the power or coverage beams with ```all```, ```power```, or ```coverage```, respectively.
  * --outRoot - the path to the directory to save plots to.
  * --site - provide a site name for plot titles.
  * --plots - call the ```plotData()``` method to produce ALS vs. GEDI comparison scatter plots.
  * --hist - include histogram plots in the call to ```plotData()```.
  * --box - include box plots in the call to ```plotData()```.
  * --dens - include scatter density plots in the call to ```plotData()```.
  * --waveforms - call the ```copyPlots()``` method produce a shell script making copies of simulated ALS and real GEDI waveform plots for waveforms meeting criteria specified in the script.
* An example command to produce the full range of plots available for the WREF site for power beams alone and save output to a subdirectory of the data directory is:

  > ```$ python3 alsComp.py --als /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l1b/pulse/metric/ --gedi /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/wref/l2a/ --beams power --outRoot ../data/wrefFigures/ --site WREF --plots --hist --box --dens```

#### collocateSubset.py
* Script to identify bounds and length of L1B HDF5 files and subset appropriately. Output is a shell script containing commands to run ```collocateWaves``` over each spatial subset of each input HDF5 file.

#### getBounds.py
* A short script to read through a .txt file of individual ALS file bounds and extract the ultimate bounds. Can convert between different CRS too.

#### gridComp.py
* Python script to read text file output from running ```gediMetric``` over simulated waveforms and extract data from TIFFs for comparison plots.

#### parseJson.py
* Python script to read the GeoJSON files produced by the NASA LP DAAC ```GEDISubsetter.py``` script and remove unnecessary datasets.

#### randomForest.py
* First attempt to predict canopy height from Sentinel-2 image using Random Forest machine learning.

#### randomForestLaSelva.py
* Edits to ```randomForest.py``` to apply machine learning to real GEDI data.

#### truncate.py
* Short Python script to strip the filenames from the end of URLs obtained from the GEDI finder site and write to a text file
as a comma separated list that can be pasted into the NASA EarthData Search page granule IDs search box.
