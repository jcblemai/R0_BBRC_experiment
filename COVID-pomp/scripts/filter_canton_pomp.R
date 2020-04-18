library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(doSNOW)
library(pomp)
library(magrittr)
library(foreach)
library(itertools)
library(lubridate)
library(parallel)
library(glue)

select <- dplyr::select
option_list = list(
  optparse::make_option(c("-c", "--config"), action="store", default='pomp_config.yaml', type='character', help="path to the config file"),
  optparse::make_option(c("-p", "--place"), action="store", default='VD', type='character', help="name of place to be run, a Canton abbrv. in CH"),
  optparse::make_option(c("-a", "--asindex"), action="store", default=0, type='numeric', help="whether to use the index of a slurm array"),
  optparse::make_option(c("-b", "--basepath"), action="store", default="COVID-pomp/", type='character', help="base path"),
  optparse::make_option(c("-j", "--jobs"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-o", "--cores"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-n", "--nfilter"), action="store", default=10, type='numeric', help="Number of filtering iterations"),
  optparse::make_option(c("-l", "--likelihood"), action="store", default='d-h', type='character', help="likelihood to be used for filtering"),
  optparse::make_option(c("-w", "--downweight"), action="store", default=0, type='numeric', help="downweight ikelihood to be used for filtering")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- yaml::read_yaml(opt$config)

source(glue("{opt$b}scripts/skellam.R"))
source(glue("{opt$b}scripts/mifCooling.R"))
source(glue("{opt$b}scripts/utils.R"))


if (opt$a == 1 & Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
  array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  opt$place <- config$places[array_id]
}

# Setup ------------------------------------------------------------------------
# Read Rscript arguments
canton <- opt$place
# Which likelihood components to use?
# deltah: balances of inputs and outputs from hospitals
# c: cases
# d: total deaths
lik_components <- str_split(opt$likelihood, "-")[[1]]#c("deltah", "c", "d")
lik_log <- str_c(str_c("ll_", lik_components), collapse = "+")
lik <- str_c(str_c("ll_", lik_components), collapse = "*")
downweight <- opt$downweight
# Test for cases in the likelihood
ll_cases <- "c" %in% lik_components
suffix <- buildSuffix(name = config$name, 
                      place = canton,
                      lik_components = lik_components, 
                      sdfrac = config$sdfrac*100, 
                      params_to_fit = config$parameters_to_fit)

filter_filename <- glue("{opt$b}results/filtered_{suffix}.rds")

liks <- read_csv(glue("{opt$b}results/loglik_exploration_{suffix}.csv"))
load(glue("{opt$b}interm/pomp_{suffix}.rda"))

source(glue("{opt$b}scripts/pomp_skeleton.R"))
cl <- makeCluster(opt$cores)
registerDoSNOW(cl)

# Filter --------------------------------------------------------------------
best_params <- liks %>%
  arrange(desc(loglik)) %>% 
  # filter(loglik > max(loglik) - 4) %>% 
  select(-contains("log")) %>%
  slice(1:3)

t3 <- system.time({
  filter_dists <- foreach(pari = iter(best_params, "row"),
                          pit = icount(nrow(best_params)),
                          .packages = c("pomp", "tidyverse", "foreach", "magrittr"),
                          .combine = rbind,
                          .noexport = c("par")) %do% {
                            foreach(it = icount(opt$nfilter),
                                    .combine = rbind,
                                    .packages = c("pomp", "tidyverse", "foreach", "magrittr")
                            ) %dopar% {
                              
                              pf <- pfilter(covid, params = as.vector(pari), Np = 1e4, filter.traj = T)
                              traj <- filter.traj(pf) %>% 
                                as.data.frame() 
                              
                              t(traj) %>%
                                as_tibble() %>% 
                                mutate(time = as.numeric(str_replace(colnames(traj), "1.", ""))) %>% 
                                gather(var, value, -time) %>% 
                                mutate(it = it,
                                       parset = pit)
                            }
                          }
})

saveRDS(filter_dists %>% mutate(ShortName = canton), file = filter_filename)

stopCluster(cl)
closeAllConnections()

cat("----- Done filtering, took", round(t3["elapsed"]/60), "mins \n")

filter_dists <- readRDS(filter_filename)
filter_stats <- filter_dists %>% 
  group_by(time, parset, var) %>% 
  summarise(mean = mean(value, na.rm = T),
            q025 = quantile(value, 0.025, na.rm = T),
            q975 = quantile(value, 0.975, na.rm = T),
            q25 = quantile(value, 0.25, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T)) %>%
  ungroup() %>% 
  mutate(date = yearsToDate(time),
         parset = factor(parset))

state_names <- unique(filter_stats$var)
plot_states <- c("tot_I", "Rt", state_names[str_detect(state_names, "a_|_curr")], "D", "H_curr") 

# Load data
data <- read_csv(glue("{opt$b}interm/data_{suffix}.csv"))
# data <- cases_data %>% dplyr::select(
#                               date, cases, deaths, hosp_incid, cum_deaths, 
#                               hosp_curr, icu_curr, discharged, delta_hosp, delta_ID, 
#                               icu_incid, deaths_icu_incid, deaths_noicu_incid, r_incid)
# data <- data %>% mutate(deaths_nohosp = deaths - deaths_icu_incid - deaths_noicu_incid)

if (canton == "CH") {
  write_csv(filter_stats, "scenario-pipeline/reports/filter_states.csv")
  write_csv(data, "scenario-pipeline/reports/national_data.csv")
}

p <- ggplot(filter_stats %>% 
              filter(var %in% plot_states,
                     time <= dateToYears(max(data$date[!is.na(data$cases)]))), 
            aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
  geom_point(data = data %>%
               rename(
                 a_I = cases,
                 a_H = hosp_incid,
                 U_curr = icu_curr,
                 H_curr = hosp_curr,
                 D = cum_deaths,
                 a_D = deaths,
                 # a_O = r_incid,
                 a_deltaH = delta_hosp,
                 a_deltaID = delta_ID,
                 # a_U = icu_incid,
                 # a_DH = deaths_noicu_incid,
                 # a_DU = deaths_icu_incid,
                 # a_DI = deaths_nohosp
               ) %>% 
               gather(var, value, -date) %>% 
               mutate(time = dateToYears(date)),
             aes(y = value)) +
  facet_wrap(~var, scales = "free")  +
  theme_bw()

p
print('ok)')
ggsave(p, filename = glue("{opt$b}results/figs/plot_{suffix}.png"), width = 9, height = 6)
