#!/bin/bash

# We assume running this from the script directory
lls=("d-c-deltah", "d-deltah")

for variant in ${lls[@]}; do

    job_file="./${variant}.job"
    SCRIPT="~/data/Javier/COVID-19_CH/COVID-pomp/scripts/covid_pomp_generic_SLURM.R"
     
    echo "#!/bin/bash
#SBATCH --job-name=${variant}.job
#SBATCH --output=logs/${variant}.%A_%a.log
#SBATCH --time=0-02:00
#SBATCH --mem=12000
#SBATCH -c 10
#SBATCH --array=1
#SBATCH -p shared

ml gcc
ml R

Rscript $SCRIPT ${variant}" > $job_file
    sbatch $job_file

done
