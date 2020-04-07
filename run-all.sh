#!/bin/bash

# We assume running this from the script directory
lls=("d-c-deltah" "d-deltah")

cantons=("BE" "BL" "BS" "FR" "GE" "GR" "JU" "LU" "NE" "TI" "VD" "VS" "ZH")
# missing: AG AI AR GL NW OW SG SH SO SZ TG UR ZG

for ll in ${lls[@]}; do
  for canton in ${cantons[@]}; do
    Rscript COVID-pomp/scripts/build_canton_pomp.R -p ${canton} -l ${ll} -r 3 -j 8
    Rscript COVID-pomp/scripts/filter_canton_pomp.R -p ${canton} -l ${ll} -j 8
  done
done
