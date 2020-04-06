#!/bin/bash

# We assume running this from the script directory
lls=("d-c-deltah" "d-deltah")
cantons=("VD" "VS" "ZH")

for ll in ${lls[@]}; do
  for canton in ${cantons[@]}; do
    Rscript COVID-pomp/scripts/covid_pomp_generic.R -p ${canton} -l ${ll}
  done
done
