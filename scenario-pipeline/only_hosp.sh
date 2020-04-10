rm -Rf hospitalization/


Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Current -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Current -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Current -d high -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Stopped -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Stopped -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s Stopped -d high -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate1 -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate1 -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate1 -d high -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate2 -d low -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate2 -d med -j 16 -c config.yml
Rscript COVIDScenarioPipeline/R/scripts/hosp_run.R -s TestIsolate2 -d high -j 16 -c config.yml

Rscript NPI/filter_data_likelihood.R

