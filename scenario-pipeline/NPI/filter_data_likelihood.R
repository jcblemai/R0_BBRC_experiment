#Preamble
library(tidyverse)
library(covidcommon)
library(report.generation)
library(kableExtra)
library(gridExtra)
library(foreach)


## Block loads the config file and geodata
config <- covidcommon:::load_config("config.yml")
scn_dirs <- paste(config$name,config$interventions$scenarios,sep='_')
scenario_labels <- config$report$formatting$scenario_labels

nfiles <- config$nsimulations## set to NULL or the actual number of sim files to include for final report
reportStateUSPS <- "CH" ## e.g. CA

col_to_filter_sim <- "NhospCurr"
col_to_filter_data<- c("NincidDeath" = "deaths", "NhospCurr" = "hosp_curr")

n_sim_target <- nfiles

# Functions --------------------------------------------------------------------

alpha <- .3    # Underreporting for cases from Imperial College estimates for CH

# Compute weights and resample
logSumExp_R <- function(lx, ...) {
  # from https://github.com/HenrikBengtsson/matrixStats/wiki/logSumExp
  iMax <- which.max(lx)
  log1p(sum(exp(lx[-iMax] - lx[iMax]))) + lx[iMax]
}

# Function to compute log-likelihood of trajectories
logLikCases <- function(cases, sim_cases, alpha, k = 1, dist = "pois") {
  if (dist == "pois") {
    lls <- dpois(cases, alpha * sim_cases, log = T)
  } else if (dist == "nbinom") {
    lls <- dnbinom(cases, mu = alpha * sim_cases, size = k, log = T)
  }
  sum(lls, na.rm = T)
}

state_hosp_totals <- load_hosp_geocombined_totals(scn_dirs,
                                                  num_files = nfiles,
                                                  scenariolabels = scenario_labels,
                                                  name_filter = "csv") %>% 
  mutate(pdeath = "", scenario_name = factor(scenario_name, levels = scenario_labels),
         usim_num = paste(scenario_num, sim_num, sep = "-"))


national_data <- read_csv("reports/national_data.csv") %>% 
  # mutate(cases = cumsum(cases)) %>%
  filter(date < "2020-04-08")


ll_data <- foreach(s = unique(state_hosp_totals$usim_num),
        .combine = rbind
) %do% {
  sim <- state_hosp_totals[state_hosp_totals$usim_num == s,c("time", col_to_filter_sim)] %>% 
    magrittr::set_colnames(c("date", "sim"))
  obs <- national_data[, c("date", col_to_filter_data[col_to_filter_sim])] %>% 
    magrittr::set_colnames(c("date", "obs"))
  
  inner_join(sim, obs) %>% 
    na.omit() %>% 
    summarise(ll = logLikCases(obs, sim + 1, k = 5, alpha = 1, "nbinom")) %>% 
    mutate(usim_num = s)
}

# We here only use the likelihood associated to hospitalizations
ll_data <- ll_data %>% 
  inner_join(distinct(state_hosp_totals, usim_num, sim_num, scenario_num)) %>% 
  group_by(scenario_num) %>% 
  filter(!is.infinite(ll)) %>% 
  mutate(w = exp(ll - logSumExp_R(ll))) %>% 
  filter(!is.na(w)) %>% 
  select(-usim_num)

# These are the resampled simulations to use
resampled_sims <- bind_rows(lapply(
  unique(ll_data$scenario_num), 
  function(x) 
    sample_n(select(ll_data, scenario_num, sim_num, w) %>% 
               filter(scenario_num == x),
             size = n_sim_target,
             replace = T,
             weight = w))) %>% 
  arrange(scenario_num, sim_num) %>% 
  select(scenario_num, sim_num) 

resampled_sims <- resampled_sims %>% 
  group_by(scenario_num) %>% 
  mutate(new_sim_num = row_number())

write_csv(resampled_sims, "reports/resample_sims.csv")

