# Dissertation

## Overview

Code developed during work on my MSc dissertation, spring and summer 2021. The various programs in this repository process GEDI and ALS data, make various comparison plots of ground elevation accuracy and associated statistics, and execute machine learning algorithms to improve the accuracy of the GEDI ground elevation estimates.

Essential Python packages that are required for most programmes are:

* *h5py*
* *NumPy*
* *matplotlib*
* *argparse*
* *gdal*
* *pyproj*
* *os*
* *glob*

Other packages may be required for certain programmes. The machine learning scripts rely on the *sklearn* package.

## Files in this repository

### In the ```shell``` folder:

#### coords.sh
* A simple shell script to execute ```gediHandler.py``` with ```--writeCoords``` option to write footprint coordinates to text files.

#### gridMetric.sh
* Shell script to execute ```gediMetric``` for simulated waveforms with addition of noise and signal processing options.

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

#### biomass.py
* Script to calculate and plot biomass estimate bias as a function of mean canopy height for temperate or tropical forests. Uses equations of Lefsky et al. (2002) for temperate/boreal forest and Asner et al. (2009) for tropical forest. The command line options available are:
  * --tempBias - the (uncorrected) bias for temperate forests.
  * --tropBias - the (uncorrected) bias for tropical forests.
  * --tempCorr - the second (corrected) bias for temperate forests.
  * --tropCorr - the second (corrected) bias for tropical forests.
  * --output - path and filename for the plot.
* An example command to plot a temperate forest bias of -1 m and a tropical forest bias of 4 m would be:

  > ```$ python3 biomass.py --tempBias -1 --tropBias 4 --output ../data/test```

#### collocateSubset.py
* Script to identify bounds and length of L1B HDF5 files and subset appropriately. Output is a shell script containing commands to run ```collocateWaves``` over each spatial subset of each input HDF5 file.

#### gediRatLoop.py
* Python script to write shell scripts for bulk execution of the ```gediRat``` programme to create a grid of simulated GEDI waveforms from ALS data.

#### getAlgorithm.py
* A short script to read through a GEDI L2A file and extract the number of occurrences of each algorithm setting, as well as other metrics.

#### getBounds.py
* A short script to read through a .txt file of individual ALS file bounds and extract the ultimate bounds. Can convert between different CRS too.

#### gridComp.py
* Python script to read text file output from running ```gediMetric``` over simulated waveforms and extract data from TIFFs for comparison plots. The command line options available are:
  * --inDir - the path to the directory containing the text file output from ```gediMetric```.
  * --inFile - the path and filename of a single text file output from ```gediMetric```.
  * --inTiff - the path and filename of the GeoTIFF image to read in.
  * --tiff - read the file specified by *--inTiff* and make comparison plots.
  * --simPlots - make histogram and box plots from the simulated data.
  * --outRoot - the output root to add to figures.
* The command to read the output from the La Selva grid and compare to the La Selva NDVI raster would be:

  > ```$ python3 gridComp.py --inFile /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/grid20/cov99wid0.5/metricAll.txt --inTiff ../data/laselva/laselvaNDVI_20m.tif --tiff --outRoot ../data/laselvaGridFigs/```

#### mapLidarSubset.py
* Script to write a shell script to execute ```mapLidar``` and generate tiled DSM and DTM rasters from ALS data.

#### myPlotComparison.py
* Edited version of ```plotComparison.py``` to make comparison plots of GEDI and ALS-simulated waveforms.

#### parseJson.py
* Python script to read the GeoJSON files produced by the NASA LP DAAC ```GEDISubsetter.py``` script and remove unnecessary datasets.

#### plotCloud.py
* Takes .pts file output from ```lasPoints``` and generates vertical section plots.

#### plottingLoop.py
* Writes a shell script executing ```myPlotComparison.py``` for the contents of a directory.

#### randomForest.py
* First attempt to predict canopy height from Sentinel-2 image using Random Forest machine learning.

#### randomForestLaSelva.py
* Edits to ```randomForest.py``` to apply machine learning to real GEDI data. Much of the control of the program must be made by direct edits to the code, however, the following command line options are available:
  * --inDir - path to the directory containing the GEDI L2A data.
  * --inDir2 - path to the directory containing the second GEDI L2A data.
  * --inFile - the path and filename of a single text file output from ```gediMetric```.
  * --inTiff - the path and filename of the GeoTIFF image to read in.
  * --inTiff2 - the path and filename of the second GeoTIFF image to read in.
  * --compare - run the ```compareGround()``` method to make comparison plots.
  * --newGround - run the ```utiliseForest()``` method to generate new ground elevation estimates.
  * --beams - specify whether to use all GEDI beams (the default), or only the power or coverage beams with ```all```, ```power```, or ```coverage```, respectively.
* To run the script with one dataset from La Selva and generate new ground estimates with the output of the RF algorithm the command would be:

  > ```$ python3 randomForestLaSelva.py --inDir /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l2a/2009/ --inFile /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/metric/2009/metricAll.txt --inTiff /exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/laselva20mComposite.tif --newGround```

#### randomForestVariogram.py
* Edits to ```randomForestLaSelva.py``` to compute semi-variogram and investigate impact of distance from training data on RF regression score. Accepts the same *--inFile*, *--inTiff*, and *--beams* options as ```randomForestLaSelva.py```.

#### truncate.py
* Short Python script to strip the filenames from the end of URLs obtained from the GEDI finder site and write to a text file
as a comma separated list that can be pasted into the NASA EarthData Search page granule IDs search box.
