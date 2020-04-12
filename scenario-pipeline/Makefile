
.PHONY: rerun rerun_simulations rerun_hospitalization clean_hospitalization clean clean_simulations

RSCRIPT=Rscript
PYTHON=python3.7
NCOREPER=16
PIPELINE=COVIDScenarioPipeline
CONFIG=config.yml


.files:
	mkdir $@
report: .files/10000_simulation_Stopped .files/10000_hospitalization_Stopped_low .files/10000_hospitalization_Stopped_med .files/10000_hospitalization_Stopped_high .files/10000_simulation_R0-1dot2 .files/10000_hospitalization_R0-1dot2_low .files/10000_hospitalization_R0-1dot2_med .files/10000_hospitalization_R0-1dot2_high .files/10000_simulation_R0-0dot9 .files/10000_hospitalization_R0-0dot9_low .files/10000_hospitalization_R0-0dot9_med .files/10000_hospitalization_R0-0dot9_high
	Rscript compile_Rmd.R
.files/10000_simulation_Stopped: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s Stopped -n 10000 -j $(NCOREPER)
	touch .files/10000_simulation_Stopped
.files/10000_hospitalization_Stopped_low: .files .files/10000_simulation_Stopped
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Stopped -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_Stopped_low
.files/10000_hospitalization_Stopped_med: .files .files/10000_simulation_Stopped
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Stopped -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_Stopped_med
.files/10000_hospitalization_Stopped_high: .files .files/10000_simulation_Stopped
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s Stopped -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_Stopped_high
.files/10000_simulation_R0-1dot2: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s R0-1dot2 -n 10000 -j $(NCOREPER)
	touch .files/10000_simulation_R0-1dot2
.files/10000_hospitalization_R0-1dot2_low: .files .files/10000_simulation_R0-1dot2
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot2 -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_R0-1dot2_low
.files/10000_hospitalization_R0-1dot2_med: .files .files/10000_simulation_R0-1dot2
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot2 -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_R0-1dot2_med
.files/10000_hospitalization_R0-1dot2_high: .files .files/10000_simulation_R0-1dot2
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-1dot2 -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_R0-1dot2_high
.files/10000_simulation_R0-0dot9: .files 
	$(PYTHON) $(PIPELINE)/simulate.py -c $(CONFIG) -s R0-0dot9 -n 10000 -j $(NCOREPER)
	touch .files/10000_simulation_R0-0dot9
.files/10000_hospitalization_R0-0dot9_low: .files .files/10000_simulation_R0-0dot9
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-0dot9 -d low -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_R0-0dot9_low
.files/10000_hospitalization_R0-0dot9_med: .files .files/10000_simulation_R0-0dot9
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-0dot9 -d med -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_R0-0dot9_med
.files/10000_hospitalization_R0-0dot9_high: .files .files/10000_simulation_R0-0dot9
	$(RSCRIPT) $(PIPELINE)/R/scripts/hosp_run.R -s R0-0dot9 -d high -j $(NCOREPER) -c $(CONFIG)
	touch .files/10000_hospitalization_R0-0dot9_high


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
