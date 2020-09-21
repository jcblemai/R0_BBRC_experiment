# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(foreach)
library(doSNOW)
library(sf)

source("COVID-pomp/scripts/graph_utils.R")

write4rep <- T
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

ffilter <- list.files(path = "COVID-pomp/results/tavary 5mars/results/", 
                      pattern = "filtered_", full.names = TRUE) %>% 
  .[str_detect(., ".rds")] %>% 
  # .[str_detect(., param_suffix)] %>% 
  .[str_detect(str_replace_all(., "COVID_CH", ""), str_c(places_to_analyze, collapse = "|"))]

# Get results
states_to_plot <- c("Rt","D", "H_curr", "U_curr", "a_D", "a_deltaH", "tot_I") 
sims <- getStates(ffilter, states_to_plot)

parsets <- tribble(
  ~ShortName, ~parset,
  "CH", 1,
  "BE", 1,
  "BL", 1,
  "BS", 1,
  "FR", 1,
  "GE", 1,
  "GR", 1,
  "JU", 1,
  "LU", 1,
  "NE", 1, 
  "TI", 2,
  # "UR", 2,
  "VD", 1,
  "VS", 1,
  "ZH", 1,
)

sims2 <- inner_join(sims, parsets)

# Extract statistics
filterstats <- computeFilterStats(sims2) %>%
  filter(!is.nan(mean))
r0_reduction <- computeR0Reduction(filter(sims2, var == "Rt"), 
                                   tw_left = as.Date(c("2020-02-25", "2020-03-03")), 
                                   tw_right = as.Date(c("2020-03-28", "2020-04-02")),
                                   date_start = as.Date("2020-03-05"))

saveRDS(filterstats, file = "COVID-pomp/results/filtered_states_all.rds")
saveRDS(r0_reduction, file = "COVID-pomp/results/R0_reduction.rds")