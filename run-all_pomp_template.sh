#!/bin/bash

# We assume running this from the script directory
lls=(%l)
places=(%c)
# missing: AG AI AR GL NW OW SG SH SO SZ TG UR ZG "LU"

for ll in ${lls[@]}; do
  for place in ${places[@]}; do
    Rscript COVID-pomp/scripts/build_canton_pomp.R -p ${place} -l ${ll} -r %r -j %j -w %w
    Rscript COVID-pomp/scripts/filter_canton_pomp.R -p ${place} -l ${ll} -j %j -n %n -w %w
  done
done
