#!/usr/bin/env bash

files=/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/*.h5

for file in $files
do
  echo "Processing file $file"
  fname="/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/footprints/$(basename -- $file .h5).txt"
  python3 $GEDIRAT_ROOT/gediHandler.py --input $file --writeCoords > $fname
done

echo "All done"
