# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(foreach)
library(doSNOW)
library(sf)

source("COVID-pomp/scripts/graph_utils.R")


# Setup ------------------------------------------------------------------------

cantons_to_analyze <- c("AG", "BE", "FR", "GE", "NE", "VD", "VS", "ZH")

# Time windows over which to compute changes in R
tw_left <- as.Date(c("2020-03-01", "2020-03-10"))
tw_right <- as.Date(c("2020-03-29", "2020-04-03"))
                   
# Load results -----------------------------------------------------------------
# Load data that has been produced
geodata <- read.csv("data/ch/geodata.csv")

files_filter <- list.files(path = "COVID-pomp/results/", pattern = "filtered_*")

cantons <- lapply(files_filter, function(s) unlist(strsplit(s, '_'))[3])
country <- st_read("data/ch/shp/swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.shp")

for (i in seq(nrow(country))) {
  country$ShortName[i] <- lapply(geodata['ShortName'][geodata$CantonNumber == country$KANTONSNUM[i],], as.character) # TODO hardcoded
}

ffilter <- list.files(path = "COVID-pomp/results/", pattern = "filtered_covid_", full.names = TRUE) %>% 
  .[str_detect(., ".rds")] %>% 
  .[str_detect(., str_c(cantons_to_analyze, collapse = "|"))]

# Get results
states_to_plot <- c("tot_I", "Rt", "D", "H_curr", "U_curr")
sims <- getStates(ffilter, states_to_plot)
# Extract statistics
filterstats <- computeFilterStats(sims)
r0_reduction <- computeR0Reduction(filter(sims, var == "Rt"), 
                                   tw_left = tw_left, 
                                   tw_right = tw_right,
                                   date_start = "2020-03-15")

saveRDS(filterstats, file = "COVID-pomp/results/filtered_states_all.rds")
saveRDS(r0_reduction, file = "COVID-pomp/results/R0_reduction.rds")
save(tw_left, tw_right, file = "COVID-pomp/results/computation_time_windows.rda")
