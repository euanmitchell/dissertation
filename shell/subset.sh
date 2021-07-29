#!/usr/bin/env bash

# The directory path containing the files to loop through as input
files=/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/raw/*.h5

# Start bash loop
for file in $files
do
  # Shell script print statement
  echo "Processing file $file"
  # Create new variable for output filename - basename strips directory path (and optionally file suffix)
  fname="/exports/csce/datastore/geos/groups/MSCGIS/s2129010/data/laselva/l1b/subset/2009/subset_$(basename -- $file)"
  # Issue the command
  python3 $GEDIRAT_ROOT/subsetGEDI.py --input $file --bounds \-84.06198 10.38395 \-83.96591 10.47477 --output $fname
# End of the loop for each file
done

echo "All done"
