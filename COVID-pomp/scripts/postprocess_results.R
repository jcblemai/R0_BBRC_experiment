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
# Load data that has been produced
geodata <- read.csv("data/ch/geodata.csv")
# files_filter <- list.files(path = "COVID-pomp/results/", pattern = "filtered_*") %>% 
  # grep(param_suffix, ., value = T)

# places <- str_extract(files_filter, "(?<=CH_)[A-Z]{2}(?=_)")
country <- st_read("data/ch/shp/ch.shp")

for (i in seq(nrow(country))) {
  country$ShortName[i] <- lapply(geodata['ShortName'][geodata$CantonNumber == country$KANTONSNUM[i],], as.character) # TODO hardcoded
}

ffilter <- list.files(path = "COVID-pomp/results/", pattern = "filtered_", full.names = TRUE) %>% 
  .[str_detect(., ".rds")] %>% 
  .[str_detect(., str_c(places_to_analyze, collapse = "|"))]

# Get results
states_to_plot <- c("tot_I", "Rt", "D", "H_curr", "U_curr")
sims <- getStates(ffilter, states_to_plot) %>% filter(ll_comp == "c-d-deltah", !is.na(value))
# Extract statistics
filterstats <- computeFilterStats(sims)  %>% filter( !is.nan(mean))
r0_reduction <- computeR0Reduction(filter(sims, var == "Rt"), 
                                   tw_left = tw_left, 
                                   tw_right = tw_right,
                                   date_start = as.Date("2020-03-15"))

saveRDS(filterstats, file = "COVID-pomp/results/filtered_states_all.rds")
saveRDS(r0_reduction, file = "COVID-pomp/results/R0_reduction.rds")
