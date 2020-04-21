#!/bin/bash

# We assume running this from the script directory
lls=('d-deltah')
#lls=('d')
#places=('CH')
places=('FR' 'GE' 'GR' 'JU' 'NE' 'TI' 'VD' 'VS' 'ZH' 'BS' 'BL')
places=('AG' 'AI' 'AR' 'GL' 'NW' 'OW' 'SG' 'SH' 'SO' 'SZ' 'TG' 'UR' 'ZG' 'LU' 'BE')
# missing: AG AI AR GL NW OW SG SH SO SZ TG UR ZG "LU" BE
for place in ${places[@]}; do
  for ll in ${lls[@]}; do
    Rscript COVID-pomp/scripts/build_canton_pomp.R -p ${place} -l ${ll} -r 3 -o 8 -j 8
    Rscript COVID-pomp/scripts/filter_canton_pomp.R -p ${place} -l ${ll} -o 8 -n 1000 
  done
done
