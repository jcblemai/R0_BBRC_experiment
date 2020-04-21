#!/bin/bash

# We assume running this from the script directory
lls=('d-deltah')
#lls=('d')
#places=('CH')
<<<<<<< HEAD
places=('FR' 'GE' 'GR' 'JU' 'NE' 'TI' 'VD' 'VS' 'ZH' 'BS' 'BL')
places=('AG' 'AI' 'AR' 'GL' 'NW' 'OW' 'SG' 'SH' 'SO' 'SZ' 'TG' 'UR' 'ZG' 'LU' 'BE')
places=("AG" "BE" "BL" "BS" "FR" "GE" "GR" "JU" "NE" "TI" "UR" "VS" "VD" "ZH")
# missing: AG AI AR GL NW OW SG SH SO SZ TG UR ZG "LU"
#'BE' 'BL' 'BS' 
>>>>>>> b55695fea3831c23dde5c7f8227a0cd2b1ade3e9
for place in ${places[@]}; do
  for ll in ${lls[@]}; do
    Rscript COVID-pomp/scripts/build_canton_pomp.R -p ${place} -l ${ll} -r 3 -o 8 -j 8
    Rscript COVID-pomp/scripts/filter_canton_pomp.R -p ${place} -l ${ll} -o 8 -n 1000 
  done
done
