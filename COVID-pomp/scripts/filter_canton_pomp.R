library(tidyverse)
library(doSNOW)
library(pomp)
library(magrittr)
library(foreach)
library(itertools)
library(lubridate)
library(parallel)
library(sf)
library(glue)

select <- dplyr::select
source("COVID-pomp/scripts/skellam.R")
source("COVID-pomp/scripts/mifCooling.R")
source("COVID-pomp/scripts/utils.R")

option_list = list(
  optparse::make_option(c("-c", "--config"), action="store", default='config.yaml', type='character', help="path to the config file"),
  optparse::make_option(c("-p", "--place"), action="store", default='BE', type='character', help="name of place to be run, a Canton abbrv. in CH"),
  optparse::make_option(c("-j", "--jobs"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-l", "--likelyhood"), action="store", default='c-d-deltah', type='character', help="likelyhood to be used for filtering")
)
opt = optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- load_config(opt$c)

# Setup ------------------------------------------------------------------------
n_filter <- 1e3
# Read Rscript arguments
canton <- opt$place
# Which likelihood components to use?
# deltah: balances of inputs and outputs from hospitals
# c: cases
# d: total deaths
lik_components <- str_split(opt$likelyhood, "-")[[1]]#c("deltah", "c", "d")
lik_log <- str_c(str_c("ll_", lik_components), collapse = "+")
lik <- str_c(str_c("ll_", lik_components), collapse = "*")
# Test for cases in the likelihood
ll_cases <- "c" %in% lik_components
suffix <- glue("{config$name}_{canton}_{str_c(lik_components, collapse = '-')}")

filter_filename <- glue("COVID-pomp/results/filtered_{suffix}.rds")

liks <- read_csv(glue("COVID-pomp/results/loglik_exploration_{suffix}.csv"))
load(glue("COVID-pomp/interm/pomp_{suffix}.rda"))
source("COVID-pomp/scripts/pomp_skeleton.R")
cl <- makeCluster(opt$jobs)
registerDoSNOW(cl)

# Filter --------------------------------------------------------------------


best_params <- liks %>%
  filter(loglik > max(loglik) - 4) %>% 
  select(-contains("log")) %>% 
  slice(1:3)

t3 <- system.time({
filter_dists <- foreach(pari = iter(best_params, "row"),
                        pit = icount(nrow(best_params)),
                        .packages = c("pomp", "tidyverse", "foreach", "magrittr"),
                        .combine = rbind,
                        .noexport = c("par")) %do% {
                          foreach(it = icount(n_filter),
                                  .combine = rbind,
                                  .packages = c("pomp", "tidyverse", "foreach", "magrittr")
                          ) %dopar% {
                            
                            pf <- pfilter(covid, params = as.vector(pari), Np = 3e3, filter.traj = T)
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

plot_states <- c("tot_I", "Rt", state_names[str_detect(state_names, "a_|_curr")], "D") 

# repeat for plot
data_file <- glue("data/ch/cases/COVID19_Fallzahlen_Kanton_{canton}_total.csv")

cases_data <- read_csv(data_file, col_types = cols()) %>% 
  mutate(cases = c(NA, diff(ncumul_conf)),
         deaths = c(NA, diff(ncumul_deceased)),
         cum_deaths = ncumul_deceased,
         hosp_curr = ncumul_hosp,
         discharged = c(ncumul_released[1], diff(ncumul_released)),
         delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
         delta_ID = delta_hosp + discharged)

data <- select(cases_data, 
               date, cases, deaths, cum_deaths, hosp_curr, discharged, delta_hosp, delta_ID) %>% 
  mutate(hosp_incid = NA)

# Set start date 
start_date <- with(data, 
                   min(c(date[which(!is.na(cases))[1]] - 5, 
                         date[which(!is.na(hosp_curr))[1]] - 8)))#as.Date("2020-02-20")

end_date <- as.Date("2020-04-08")

# Add rows
data <- rbind(tibble(date = seq.Date(start_date, min(data$date), by = "1 days")) %>% 
                cbind(data[1, -1] %>% mutate_all(function(x) x <- NA)) , data)


p <- ggplot(filter_stats %>% filter(var %in% plot_states), aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
  # geom_line(aes(y = mean)) +
  geom_point(data = data %>%
               mutate(time = dateToYears(date),
                      cases = cases / best_params[["epsilon"]][1],
                      deaths = case_when(is.na(deaths) ~ 0,T ~  deaths),
                      cum_deaths = cumsum(deaths)) %>%
               rename(
                 a_I = cases,
                 # a_DU = deaths_icu_incid,
                 # a_DH = deaths_noicu_incid,
                 # a_DI = deaths_nohosp_incid,
                 a_H = hosp_incid,
                 # a_U = icu_incid,
                 # U_curr = icu_curr,
                 H_curr = hosp_curr,
                 D = cum_deaths,
                 a_D = deaths,
                 a_O = discharged,
                 a_deltaH = delta_hosp,
                 a_deltaID = delta_ID
               ) %>% 
               gather(var, value, -time, -date),
             aes(y = value)) +
  facet_wrap(~var, scales = "free")  +
  theme_bw()
print('ok)')
ggsave(p, filename = glue("COVID-pomp/results/figs/plot_{suffix}.png"), width = 9, height = 6)
