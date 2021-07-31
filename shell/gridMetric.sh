#!/usr/bin/env bash

# Run gediMetric over the output of gediRat grid sims

files=/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/grid20/*.h5

for file in $files
do
  echo "Processing file $file"
  root="/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/grid20/cov99wid0.5/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -linkNoise 0 0.99 -varScale 3 -sWidth 0.5 -minWidth 3 -outRoot $root
done

echo "All done"
