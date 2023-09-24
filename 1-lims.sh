#!/bin/bash

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project
/home/groups/singlecell/mabdalfttah/projects/scripts/limsq.py -sp $1 | sed 's/;/\t/g' > "lims_info_"$1".txt"

echo "Created LIMS information file: lims_info.txt"
