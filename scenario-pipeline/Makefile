
.PHONY: rerun rerun_simulations rerun_hospitalization clean_hospitalization clean clean_simulations

RSCRIPT=Rscript
PYTHON=python3.7
NCOREPER=4
PIPELINE=COVIDScenarioPipeline
CONFIG=config.yml


.files:
	mkdir $@
report: .files/100_simulation_Wuhan .files/100_hospitalization_Wuhan_low .files/100_hospitalization_Wuhan_med .files/100_hospitalization_Wuhan_high .files/100_simulation_None .files/100_hospitalization_None_low .files/100_hospitalization_None_med .files/100_hospitalization_None_high
	Rscript compile_Rmd.R
.files/100_simulation_Wuhan: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s Wuhan -n 100 -j $(NCOREPER)
	touch .files/100_simulation_Wuhan
.files/100_hospitalization_Wuhan_low: .files .files/100_simulation_Wuhan
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Wuhan -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/100_hospitalization_Wuhan_low
.files/100_hospitalization_Wuhan_med: .files .files/100_simulation_Wuhan
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Wuhan -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/100_hospitalization_Wuhan_med
.files/100_hospitalization_Wuhan_high: .files .files/100_simulation_Wuhan
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Wuhan -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/100_hospitalization_Wuhan_high
.files/100_simulation_None: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s None -n 100 -j $(NCOREPER)
	touch .files/100_simulation_None
.files/100_hospitalization_None_low: .files .files/100_simulation_None
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s None -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/100_hospitalization_None_low
.files/100_hospitalization_None_med: .files .files/100_simulation_None
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s None -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/100_hospitalization_None_med
.files/100_hospitalization_None_high: .files .files/100_simulation_None
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s None -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/100_hospitalization_None_high


rerun: rerun_simulations rerun_hospitalization
rerun_filter:
	rm -f .files/1*_filter
rerun_importation:
	rm -f .files/1*_importation
rerun_simulations: clean_simulations
	rm -f .files/1*_simulation*
rerun_hospitalization:
	rm -f .files/1*_hospitalization*
clean: clean_simulations clean_hospitalization clean_reports
	rm -rf .files
clean_reports:
	rm -f notebooks/*.html
clean_simulations: rerun_simulations
	rm -rf model_output
clean_hospitalization: rerun_hospitalization
	rm -rf hospitalization
