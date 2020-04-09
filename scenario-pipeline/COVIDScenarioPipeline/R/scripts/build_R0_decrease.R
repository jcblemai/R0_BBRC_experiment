
# Preamble ---------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(sf)
library(doSNOW)


option_list <- list(
  optparse::make_option(c("-c", "--config"), action="store", default='config.yml', type='character', help="path to the config file"),
  optparse::make_option(c("-i", "--input"), action="store", default='COVID-pomp/results/archive/', type='character', help="path to the config file")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- covidcommon::load_config(opt$c)
pomp_config <- yaml::read_yaml("../pomp_config.yaml")

source("../COVID-pomp/scripts/graph_utils.R")
# Data -------------------------------------------------------------------------

# What cantons to keep results for
keep_cantons <- c("BE", "FR", "GE", "GR", "JU", "NE", "TI", "VD", "VS", "ZH")

# Load data from canton inference
# places <- str_extract(files_filter, "(?<=CH_)[A-Z]{2}(?=_)")
geom_file <- paste(config$spatial_setup$base_path, config$spatial_setup$shapefile, sep = "/")
country <- st_read(geom_file)
geodata <- read_csv(paste(config$spatial_setup$base_path, config$spatial_setup$geodata, sep = "/"))
for (i in seq(nrow(country))) {
  country$ShortName[i] <- lapply(geodata['ShortName'][geodata$CantonNumber == country$KANTONSNUM[i],], as.character) # TODO hardcoded
}

all_cantons <- sort(unique(country[["ShortName"]]))

ffilter <- list.files(path = "../COVID-pomp/results/", pattern = "filtered_", full.names = TRUE) %>% 
  .[str_detect(., ".rds")] %>% 
  .[str_detect(., "c-")] %>% 
  .[str_detect(., str_c(keep_cantons, collapse = "|"))]

# Filtered trajectories of R0
sims <- getStates(ffilter, "Rt")
# Reduction in R0
r0_reduction <- computeR0Reduction(sims, 
                                   as.Date(pomp_config$timewindow_R0$left), 
                                   as.Date(pomp_config$timewindow_R0$right), as.Date("2020-03-15")) %>% 
  filter(str_detect(var, "r0"))

write_csv(r0_reduction, "../data/ch/r0_reduction.csv")


all_cantons <- sort(unique(country[["ShortName"]]))

# Seed ---------------------------------------------------------------

seed_dates <- foreach(cnt = all_cantons,
                      .combine = rbind
) %do% {
  
  data_file <- glue("../data/ch/cases/covid_19/fallzahlen_kanton_total_csv/COVID19_Fallzahlen_Kanton_{cnt}_total.csv")
  
  cases_data <- read_csv(data_file, col_types = cols()) %>% 
    mutate(cases = c(NA, diff(ncumul_conf)),
           deaths = c(NA, diff(ncumul_deceased)),
           cum_deaths = ncumul_deceased,
           hosp_curr = ncumul_hosp,
           discharged = c(ncumul_released[1], diff(ncumul_released)),
           delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
           delta_ID = delta_hosp + discharged)
  
  
  start_date <- cases_data$date[1] - 5
  
  if (cnt %in% keep_cantons) {
    canton <- cnt
    lik_components <- str_split(pomp_config$likelihoods$to_plot, "-")[[1]]#c("deltah", "c", "d")
    lik_log <- str_c(str_c("ll_", lik_components), collapse = "+")
    lik <- str_c(str_c("ll_", lik_components), collapse = "*")
    # Test for cases in the likelihood
    ll_cases <- "c" %in% lik_components
    suffix <- glue("{pomp_config$name}_{canton}_{str_c(lik_components, collapse = '-')}_{ifelse(is.null(pomp_config$parameters_to_fit), '', str_c(names(pomp_config$parameters_to_fit), collapse = '-'))}")
    liks <- read_csv(glue("../COVID-pomp/results/loglik_exploration_{suffix}.csv"))
    
    best_params <- liks %>%
      filter(loglik > max(loglik) - 4) %>% 
      select(-contains("log")) %>% 
      slice(1)
    
    amount <- best_params$I_0 * geodata$pop2018[geodata$ShortName == cnt]
  } else {
    amount <- max(c(5*cases_data$ncumul_conf[1], cases_data$ncumul_hosp[1]), na.rm = T)
    if( amount == 0 ){
      amount <- 5
    }
  }
  tibble(place = cnt, date = start_date, amount = amount)
}

write_csv(seed_dates, "../data/ch/seeding.csv")

