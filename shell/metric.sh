#!/usr/bin/env bash

files=/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/sim/*.h5

for file in $files
do
  echo "Processing file $file"
  fname="/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/test/$(basename -- $file .h5)"
  gediMetric -input $file -readHDFgedi -ground -statsLen 8 -outRoot $fname
done

echo "All done"
