# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(foreach)
library(doSNOW)
library(sf)

source("COVID-pomp/scripts/graph_utils.R")


# Setup ------------------------------------------------------------------------
option_list = list(
  optparse::make_option(c("-c", "--config"), action="store", default='pomp_config.yaml', type='character', help="path to the config file")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- covidcommon::load_config(opt$c)

param_suffix <- ifelse(is.null(config$parameters_to_fit), '', str_c(names(config$parameters_to_fit), collapse = '-'))
places_to_analyze <- config$places
  
# Time windows over which to compute changes in R
tw_left <- as.Date(config$timewindow_R0$left)
tw_right <- as.Date(config$timewindow_R0$right)
                   
# Load results -----------------------------------------------------------------

ffilter <- list.files(path = "COVID-pomp/results/good_run/", 
                      pattern = "filtered_", full.names = TRUE) %>% 
  .[str_detect(., ".rds")] %>% 
  .[str_detect(., str_c(places_to_analyze, collapse = "|"))]

# Get results
states_to_plot <- c("tot_I", "Rt", "D", "H_curr", "U_curr")
sims <- getStates(ffilter, states_to_plot) %>% 
  filter(ll_comp == config$likelihoods$to_plot, !is.na(value))
# Extract statistics
filterstats <- computeFilterStats(sims)  %>% filter( !is.nan(mean))
r0_reduction <- computeR0Reduction(filter(sims, var == "Rt"), 
                                   tw_left = tw_left, 
                                   tw_right = tw_right,
                                   date_start = as.Date("2020-03-15"))

saveRDS(filterstats, file = "COVID-pomp/results/filtered_states_all.rds")
saveRDS(r0_reduction, file = "COVID-pomp/results/R0_reduction.rds")

if (write4rep) {
  filterstats %>% filter(paramset == 1, ShortName == "CH") %>% 
    write_csv("scenario-pipeline/reports/filter_states.csv")
  
  r0_reduction %>% filter(paramset == 1, ShortName == "CH") %>% 
    write_csv("scenario-pipeline/reports/R0_reduction.csv")
}