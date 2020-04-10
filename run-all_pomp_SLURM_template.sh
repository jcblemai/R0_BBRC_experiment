#!/bin/bash

# We assume running this from the script directory
lls=(%l)
places=(%c)

for variant in ${lls[@]}; do

    job_file="./${variant}.job"
    SCRIPT1="~/data/Javier/COVID-19_CH/COVID-pomp/scripts/build_canton_pomp.R"
    SCRIPT2="~/data/Javier/COVID-19_CH/COVID-pomp/scripts/filter_canton_pomp.R"
    BASE="~/data/Javier/COVID-19_CH/"
    echo "#!/bin/bash
#SBATCH --job-name=${variant}.job
#SBATCH --output=logs/${variant}.%A_%a.log
#SBATCH --time=0-01:00
#SBATCH --mem=10000
#SBATCH -c %j
#SBATCH --array=%a
#SBATCH -p shared

ml gcc
ml R

Rscript $SCRIPT ${variant}
Rscript $SCRIPT1  -b $BASE --asindex 1 -l ${variant} -r %r 
Rscript $SCRIPT2  -b $BASE --asindex 1 -l ${variant} -n %n"  > $job_file
    sbatch $job_file

done

