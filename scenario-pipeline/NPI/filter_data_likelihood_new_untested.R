
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(covidcommon)
library(report.generation)
library(foreach)
library(magrittr)
library(itertools)

## Block loads the config file and geodata
config <- covidcommon:::load_config("config_filter.yml")
scn_dirs <- paste(config$name,config$interventions$scenarios,sep='_')
scenario_labels <- config$report$formatting$scenario_labels

nfiles <- 100#config$nsimulations## set to NULL or the actual number of sim files to include for final report
reportStateUSPS <- "CH" ## e.g. CA

var_to_filter <- config$report$filtering$var_to_filter

sim_data_dict <- c("NincidDeath" = "deaths",
                   "NhospCurr" = "hosp_curr",
                   "NincidHosp" = "hosp_incid")

n_sim_target <- nfiles
data_path <- config$report$filtering$data_path
out_file <- paste0("model_output/sim_resampling_", paste(config$interventions$scenarios, collapse = "-"), ".csv")

# Functions --------------------------------------------------------------------

# Compute weights and resample
logSumExp_R <- function(lx, ...) {
  # from https://github.com/HenrikBengtsson/matrixStats/wiki/logSumExp
  iMax <- which.max(lx)
  log1p(sum(exp(lx[-iMax] - lx[iMax]))) + lx[iMax]
}

# Function to compute log-likelihood of trajectories
logLikObs <- function(obs, sim, alpha = 1, k = 1, dist = "pois", replace_zero = F) {
  if(replace_zero)
    sim[sim==0] <- 1
  
  if (dist == "pois") {
    lls <- dpois(obs, alpha * sim, log = T)
  } else if (dist == "nbinom") {
    lls <- dnbinom(obs, mu = alpha * sim, size = k, log = T)
  }
  sum(lls, na.rm = T)
}


# Load sims and data -----------------------------------------------------------
# Simulation outputs
state_hosp_totals <- list()
for (i in 1:length(config$hospitalization$parameters$p_death_names)) {
  state_hosp_totals[[i]] <- load_hosp_geocombined_totals(scn_dirs[i],
                                                         num_files = nfiles,
                                                         scenariolabels = config$report$formatting$scenario_labels,
                                                         # name_filter= config$hospitalization$parameters$p_death_names[i]
                                                         name_filter = "csv") %>%
    mutate(pdeath=config$hospitalization$parameters$p_death[i],
           usim_num = paste(scenario_num, i, sim_num, sep = "-"))
}

state_hosp_totals <- dplyr::bind_rows(state_hosp_totals)

# Data
data <- read_csv(data_path)

obs <- data %>% 
  select(date, one_of(sim_data_dict[var_to_filter])) %>% 
  set_colnames(c("time", var_to_filter))

# Compute log-likelihood of data for each sim
# This part can be parallelized
ll_data <- foreach(sim = isplit(state_hosp_totals, state_hosp_totals$usim_num),
                   .combine = rbind,
                   .packages = c("dplyr", "magrittr")
) %do% {
  
  combined <- sim$value %>% 
    select(time, pdeath, one_of(var_to_filter)) %>% 
    inner_join(obs, ., by = "time", suffix = c(".obs", ".sim")) 
  
  res <- lapply(seq_along(var_to_filter), function(i)
    sum(logLikObs(combined[[paste0(var_to_filter[i], ".obs")]],
                  combined[[paste0(var_to_filter[i], ".sim")]],
                  dist = config$report$filtering$dist[i],                  # distribution to use for observation model
                  alpha = config$report$filtering$under_reporting,         # under reporting of observations
                  k = config$report$filtering$k[i],                        # aggregations parameter for NB distribution
                  replace_zero = config$report$filtering$replace_zero[i])) # replace zeros in sims to avoid Inf values
    ) %>% 
    as.data.frame() %>% 
    set_colnames(var_to_filter) %>% 
    mutate(usim_num = sim$key[[1]])
  
  # compute total loglik
  res$ ll <- sum(res[1, var_to_filter])
  res
}

# Compute likelihood weights
ll_data <- ll_data %>% 
  inner_join(distinct(state_hosp_totals, usim_num, sim_num, scenario_num, pdeath)) %>% 
  group_by(scenario_num) %>% 
  filter(!is.infinite(ll)) %>% 
  mutate(w = exp(ll - logSumExp_R(ll))) %>% 
  filter(!is.na(w)) %>% 
  select(-usim_num)

# These are the resampled simulations to use
resampled_sims <- bind_rows(lapply(
  unique(ll_data$scenario_num), 
  function(x) 
    sample_n(select(ll_data, scenario_num, pdeath, sim_num, w) %>% 
               filter(scenario_num == x),
             size = n_sim_target,
             replace = T,
             weight = w))) %>% 
  arrange(scenario_num, pdeath, sim_num) %>% 
  select(scenario_num, pdeath, sim_num) 

resampled_sims <- resampled_sims %>% 
  group_by(scenario_num, pdeath) %>% 
  mutate(new_sim_num = row_number())

write_csv(resampled_sims, out_file)

