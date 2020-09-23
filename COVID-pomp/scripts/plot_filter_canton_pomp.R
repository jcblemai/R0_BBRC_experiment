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
  make_option(c("-p", "--place"), default='CH', type='character', help="name of place to be run, a place abbrv. in CH"),
  make_option(c("-b", "--basepath"), default="COVID-pomp/", type='character', help="base path"),
  make_option(c("-l", "--likelihood"), default='d-deltah', type='character', help="likelihood to be used for filtering"),
  make_option(c("-E", "--nbE"), default = 1, type = "numeric", help = "number of E compartement"),
  make_option(c("-I", "--nbI"), default = 3, type = "numeric", help = "number of I compartement"),
  make_option(c("-s", "--suffix"), default = NULL, type = "character", help = "custom suffix to add")
)


opt <-parse_args(OptionParser(option_list=option_list))
config <- yaml::read_yaml(opt$config)

nc_E <- opt$nbE
nc_I <- opt$nbI

source(glue("{opt$b}scripts/utils.R"))


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
  other =  c(config$sdfrac * 100, opt$suffix, 'NB', nc_E, nc_I)
)

filter_filename <- glue("{opt$b}results/filtered_{suffix}.rds")
liks <- read_csv(glue("{opt$b}results/loglik_exploration_{suffix}.csv"))

load(glue("{opt$b}interm/pomp_{suffix}.rda"))
source(glue("{opt$b}scripts/{config$model}"))

# Do plot ----------------------------------------------------------------------

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


tstart <- dateToYears(max(epidata$date[!is.na(epidata$a_I)], na.rm = T))
                      
p <- ggplot(filter_stats %>% 
              filter(var %in% plot_states,
                     time <= tstart), 
            aes(x = time)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
  geom_point(data = epidata4plot,  aes(y = value), size = 1) +
  facet_wrap(~var, scales = "free")  +
  theme_bw() +
  scale_x_continuous(labels = yearsToDateLabel)
p
print('ok)')
ggsave(p, filename = glue("{opt$b}results/figs/plot_{suffix}.png"), width = 9, height = 6)
