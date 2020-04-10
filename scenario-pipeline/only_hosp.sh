rm -Rf hospitalization/


Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Current -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Current -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Current -d high -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Stopped -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Stopped -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Stopped -d high -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate -d high -j 16 -c config.yml

