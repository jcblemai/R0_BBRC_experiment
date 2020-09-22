#!/bin/bash

# We assume running this from the script directory
lls=('d-deltah')
#lls=('d')
#places=('CH')

#places=('FR' 'GE' 'GR' 'JU' 'NE' 'TI' 'VD' 'VS' 'ZH' 'BS' 'BL')
#places=('AG' 'AI' 'AR' 'GL' 'NW' 'OW' 'SG' 'SH' 'SO' 'SZ' 'TG' 'UR' 'ZG' 'LU' 'BE')
#places=("BE" "BL" "BS" "FR" "GE" "GR" "JU" "NE" "TI" "VS" "VD" "ZH")
places=("CH")
# missing: AG AI AR GL NW OW SG SH SO SZ TG UR ZG "LU"
#'BE' 'BL' 'BS'

for place in ${places[@]}; do
  for ll in ${lls[@]}; do
    Rscript COVID-pomp/scripts/build_canton_pomp.R -p ${place} -l ${ll} -r 3 -o 8 -j 16 -s gamma1-3
    Rscript COVID-pomp/scripts/filter_canton_pomp.R -p ${place} -l ${ll} -o 8 -n 1000 -s gamma1-3
    Rscript COVID-pomp/scripts/plot_filter_canton_pomp.R -p ${place} -l ${ll} -s gamma1-3
  done
done
