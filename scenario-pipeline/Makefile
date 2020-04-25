
.PHONY: rerun rerun_simulations rerun_hospitalization clean_hospitalization clean clean_simulations

RSCRIPT=Rscript
PYTHON=python3.7
NCOREPER=12
PIPELINE=COVIDScenarioPipeline
CONFIG=config.yml


.files:
	mkdir $@
report: .files/1000_simulation_Stopped .files/1000_hospitalization_Stopped_low .files/1000_hospitalization_Stopped_med .files/1000_hospitalization_Stopped_high .files/1000_simulation_R0-1dot2 .files/1000_hospitalization_R0-1dot2_low .files/1000_hospitalization_R0-1dot2_med .files/1000_hospitalization_R0-1dot2_high .files/1000_simulation_R0-1dot5 .files/1000_hospitalization_R0-1dot5_low .files/1000_hospitalization_R0-1dot5_med .files/1000_hospitalization_R0-1dot5_high
	Rscript compile_Rmd.R
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
.files/1000_simulation_R0-1dot2: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s R0-1dot2 -n 1000 -j $(NCOREPER)
	touch .files/1000_simulation_R0-1dot2
.files/1000_hospitalization_R0-1dot2_low: .files .files/1000_simulation_R0-1dot2
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot2 -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_R0-1dot2_low
.files/1000_hospitalization_R0-1dot2_med: .files .files/1000_simulation_R0-1dot2
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot2 -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_R0-1dot2_med
.files/1000_hospitalization_R0-1dot2_high: .files .files/1000_simulation_R0-1dot2
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot2 -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_R0-1dot2_high
.files/1000_simulation_R0-1dot5: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s R0-1dot5 -n 1000 -j $(NCOREPER)
	touch .files/1000_simulation_R0-1dot5
.files/1000_hospitalization_R0-1dot5_low: .files .files/1000_simulation_R0-1dot5
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot5 -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_R0-1dot5_low
.files/1000_hospitalization_R0-1dot5_med: .files .files/1000_simulation_R0-1dot5
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot5 -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_R0-1dot5_med
.files/1000_hospitalization_R0-1dot5_high: .files .files/1000_simulation_R0-1dot5
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot5 -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/1000_hospitalization_R0-1dot5_high


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
