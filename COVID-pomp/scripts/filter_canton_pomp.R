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
library(optparse)

select <- dplyr::select

option_list = list(
  make_option(c("-c", "--config"), default='pomp_config.yaml', type='character', help="path to the config file"),
  make_option(c("-p", "--place"), default='VD', type='character', help="name of place to be run, a place abbrv. in CH"),
  make_option(c("-a", "--asindex"), default=0, type='numeric', help="whether to use the index of a slurm array"),
  make_option(c("-b", "--basepath"), default="COVID-pomp/", type='character', help="base path"),
  make_option(c("-j", "--jobs"), default=detectCores(), type='numeric', help="number of cores used"),
  make_option(c("-o", "--cores"), default=detectCores(), type='numeric', help="number of cores used"),
  make_option(c("-r", "--run_level"), default = 1, type = "numeric", help = "run level for MIF"),
  make_option(c("-n", "--nfilter"), default=10, type='numeric', help="Number of filtering iterations"),
  make_option(c("-l", "--likelihood"), default='d-deltah', type='character', help="likelihood to be used for filtering"),
  make_option(c("-s", "--suffix"), default = NULL, type = "character", help = "custom suffix to add")
)

opt <-parse_args(OptionParser(option_list=option_list))
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
place <- opt$place

# Parse likelihood components to use
parsed_lik <- parseLikelihood(opt$likelihood)

# Check cases in the likelihood
use_case_incid <- "c" %in% parsed_lik$components

suffix <- buildSuffix(
  name = config$name,
  place = place,
  lik_components = parsed_lik$components,
  params_to_fit = config$parameters_to_fit,
  other =  c(config$sdfrac * 100, opt$suffix)
)

filter_filename <- glue("{opt$b}results/filtered_{suffix}.rds")
liks <- read_csv(glue("{opt$b}results/loglik_exploration_{suffix}.csv"))

load(glue("{opt$b}interm/pomp_{suffix}.rda"))
source(glue("{opt$b}scripts/{config$model}"))

sir_NpLL <- c(1e2, 1e4, 1e4)

cl <- makeCluster(opt$cores)
registerDoSNOW(cl)

# Filter --------------------------------------------------------------------
best_params <- liks %>%
  arrange(desc(loglik)) %>% 
  select(-contains("log")) %>%
  slice(1)

t3 <- system.time({
  filter_dists <- foreach(pari = iter(best_params, "row"),
                          pit = icount(nrow(best_params)),
                          .combine = rbind,
                          .noexport = c("par")
  ) %do% {
    foreach(it = icount(opt$nfilter),
            .combine = rbind,
            .packages = c("pomp", "tidyverse", "foreach", "magrittr")
    ) %dopar% {
      
      pf <- pfilter(covid, params = as.vector(pari), Np = sir_NpLL[opt$run_level], filter.traj = T)
      traj <- filter.traj(pf) %>% 
        as.data.frame() 
      
      t(traj) %>%
        as.data.frame() %>% 
        mutate(time = as.numeric(str_replace(colnames(traj), "1.", ""))) %>% 
        gather(var, value, -time) %>% 
        mutate(it = it,
               parset = pit)
    }
  }
})

saveRDS(filter_dists %>% mutate(ShortName = place), file = filter_filename)

stopCluster(cl)
closeAllConnections()


# Do plot ----------------------------------------------------------------------

cat("----- Done filtering, took", round(t3["elapsed"]/60), "mins \n")

filter_dists <- readRDS(filter_filename)
filter_stats <- filter_dists %>% 
  group_by(time, parset, var) %>% 
  summarise(mean = mean(value, na.rm = T),
            median = median(value, na.rm = T),
            q025 = quantile(value, 0.025, na.rm = T),
            q975 = quantile(value, 0.975, na.rm = T),
            q25 = quantile(value, 0.25, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T)) %>%
  ungroup() %>% 
  mutate(date = yearsToDate(time),
         parset = factor(parset))

state_names <- unique(filter_stats$var)
plot_states <- c("tot_I", "Rt", state_names[str_detect(state_names, "a_|_curr")], "D", "H_curr") 

comp_dict <- c(a_I = "case_incid",
               a_H = "hosp_incid",
               U_curr = "icu_curr",
               H_curr = "hosp_curr",
               D = "cum_deaths",
               a_D = "death_incid",
               a_O = "discharge_incid",
               a_deltaH = "delta_hosp",
               a_deltaID = "delta_ID",
               a_deltaU = "delta_icu",
               a_U = "icu_incid",
               a_DH = "death_noicu_incid",
               a_DU = "death_icu_incid",
               a_DI = "death_nohosp_incid")

# Load epidata
epidata <- read_csv(glue("{opt$b}interm/data_{suffix}.csv"))%>%
  rename(!!comp_dict) 
epidata4plot <- epidata %>% 
  select(-date)  %>% 
  gather(var, value, -time) 

epidata2 <- read_csv("data/VD_hosp_data.csv") %>% select(date, contains("incid"), contains("curr"))
epidata4plot2 <- epidata2 %>%
  mutate(time = dateToYears(date)) %>% 
  select(-date) %>% 
  rename(a_H = hosp_incid,
         a_U = icu_incid,
         a_DH = deaths_noicu_incid,
         a_DU = deaths_icu_incid,
         a_O = r_incid,
         H_curr = hosp_curr,
         U_curr = icu_curr) %>% 
  gather(var, value, -time) 


epidata3 <- read_csv("data/vd/current_hosp_official.csv") 
epidata4plot3 <- epidata3 %>%
  mutate(time = dateToYears(as.Date(date, "%m/%d/%Y"))) %>% 
  select(time, hosp_curr, icu_curr, deaths) %>% 
  arrange(time) %>% 
  mutate(a_deltaH = c(hosp_curr[1], diff(hosp_curr))) %>% 
  rename(H_curr = hosp_curr,
         U_curr = icu_curr,
         D = deaths) %>% 
  gather(var, value, -time) 


epidata4 <- read_csv("data/vd/current_hosp_VD_2.csv") %>%
  mutate(time = dateToYears(date)) %>% 
  arrange(time) %>% 
  mutate(a_deltaH = c(hosp_curr[1], diff(hosp_curr)),
         a_deltaH = case_when(a_deltaH< -100 ~ as.numeric(NA), T ~ a_deltaH),
         a_deltaU = c(icu_curr[1], diff(icu_curr)),
         a_deltaU = case_when(a_deltaU < -40 ~ as.numeric(NA), T ~ a_deltaU)) %>% 
  rename(H_curr = hosp_curr,
         U_curr = icu_curr)
epidata4plot4 <- epidata4  %>% 
  select(-date) %>% 
  gather(var, value, -time) 

if (place == "CH") {
  write_csv(filter_stats, "scenario-pipeline/reports/filter_states.csv")
  write_csv(epidata, "scenario-pipeline/reports/national_epidata.csv")
}


tstart <- dateToYears(max(epidata$date[!is.na(epidata$case_incid)]))
                      
p <- ggplot(filter_stats %>% 
              filter(var %in% plot_states,
                     time <= tstart), 
            aes(x = time)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
  geom_point(data = epidata4plot,  aes(y = value), size = 1) +
  geom_point(data = epidata4plot2, aes(y = value), col = "blue", size = 1) +
  geom_point(data = epidata4plot3, aes(y = value), col = "red", size = 1) +
  geom_point(data = epidata4plot4, aes(y = value), col = "green", size = 1) +
  facet_wrap(~var, scales = "free")  +
  theme_bw() +
  scale_x_continuous(labels = yearsToDateLabel)
p
print('ok)')
ggsave(p, filename = glue("{opt$b}results/figs/plot_{suffix}.png"), width = 9, height = 6)
