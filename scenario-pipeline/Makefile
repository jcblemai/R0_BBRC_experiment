
.PHONY: rerun rerun_simulations rerun_hospitalization clean_hospitalization clean clean_simulations

RSCRIPT=Rscript
PYTHON=python3.7
NCOREPER=16
PIPELINE=COVIDScenarioPipeline
CONFIG=config.yml


.files:
	mkdir $@
report: .files/1000_simulation_Current .files/1000_hospitalization_Current_low .files/1000_hospitalization_Current_med .files/1000_hospitalization_Current_high .files/1000_simulation_Stopped .files/1000_hospitalization_Stopped_low .files/1000_hospitalization_Stopped_med .files/1000_hospitalization_Stopped_high .files/1000_simulation_TestIsolate .files/1000_hospitalization_TestIsolate_low .files/1000_hospitalization_TestIsolate_med .files/1000_hospitalization_TestIsolate_high
	Rscript compile_Rmd.R
.files/1000_simulation_Current: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s Current -n 1000 -j $(NCOREPER)
	touch .files/1000_simulation_Current
.files/1000_hospitalization_Current_low: .files .files/1000_simulation_Current
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Current -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_Current_low
.files/1000_hospitalization_Current_med: .files .files/1000_simulation_Current
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Current -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_Current_med
.files/1000_hospitalization_Current_high: .files .files/1000_simulation_Current
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Current -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_Current_high
.files/1000_simulation_Stopped: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s Stopped -n 1000 -j $(NCOREPER)
	touch .files/1000_simulation_Stopped
.files/1000_hospitalization_Stopped_low: .files .files/1000_simulation_Stopped
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Stopped -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_Stopped_low
.files/1000_hospitalization_Stopped_med: .files .files/1000_simulation_Stopped
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Stopped -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_Stopped_med
.files/1000_hospitalization_Stopped_high: .files .files/1000_simulation_Stopped
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Stopped -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_Stopped_high
.files/1000_simulation_TestIsolate: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s TestIsolate -n 1000 -j $(NCOREPER)
	touch .files/1000_simulation_TestIsolate
.files/1000_hospitalization_TestIsolate_low: .files .files/1000_simulation_TestIsolate
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s TestIsolate -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_TestIsolate_low
.files/1000_hospitalization_TestIsolate_med: .files .files/1000_simulation_TestIsolate
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s TestIsolate -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_TestIsolate_med
.files/1000_hospitalization_TestIsolate_high: .files .files/1000_simulation_TestIsolate
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s TestIsolate -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_TestIsolate_high


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
