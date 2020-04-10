# Title: POMP model of COVID19 in CH, cantonal level
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
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

option_list = list(
  optparse::make_option(c("-c", "--config"), action="store", default='pomp_config.yaml', type='character', help="path to the config file"),
  optparse::make_option(c("-a", "--asindex"), action="store", default=0, type='numeric', help="whether to use the index of a slurm array"),
  optparse::make_option(c("-b", "--basepath"), action="store", default="", type='character', help="base path"),
  optparse::make_option(c("-j", "--jobs"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-o", "--cores"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-r", "--run_level"), action="store", default=1, type='numeric', help="run level for MIF"),
  optparse::make_option(c("-p", "--place"), action="store", default='AR', type='character', help="name of place to be run, a Canton abbrv. in CH"),
  optparse::make_option(c("-l", "--likelihood"), action="store", default='c-d-deltah', type='character', help="likelihood to be used for filtering"),
  optparse::make_option(c("-w", "--downweight"), action="store", default=0, type='numeric', help="downweight ikelihood to be used for filtering")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- yaml::read_yaml(opt$config)

source(glue("{opt$b}COVID-pomp/scripts/skellam.R"))
source(glue("{opt$b}COVID-pomp/scripts/mifCooling.R"))
source(glue("{opt$b}COVID-pomp/scripts/utils.R"))

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
suffix <- buildSuffix(config$name, canton, lik_components, config$parameters_to_fit)


# Level of detail on which to run the computations
run_level <- opt$run_level
sir_Np <- c(1e3, 3e3, 3e3)
sir_Nmif <- c(2, 20, 150)
sir_Ninit_param <- c(opt$jobs/2, opt$jobs, opt$jobs)
sir_NpLL <- c(1e3, 1e4, 1e4)
sir_Nreps_global <- c(2, 5, 10)


# Data ---------------------------------------------------------------

data_file <- glue("{opt$b}data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2/COVID19_Fallzahlen_Kanton_{canton}_total.csv")

cases_data <- read_csv(data_file, col_types = cols()) %>% 
  filter(date < "2020-04-08") %>%
  mutate(cases = c(NA, diff(ncumul_conf)),
         deaths = c(NA, diff(ncumul_deceased)),
         cum_deaths = ncumul_deceased,
         hosp_incid = new_hosp,
         icu_curr = current_icu,
         hosp_curr = current_hosp,
         discharged = c(ncumul_released[1], diff(ncumul_released)),
         delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
         delta_ID = delta_hosp + discharged)

cases_data$cases[cases_data$cases < 0] <- NA

if (canton == "CH") {
  cases_data$delta_ID <- NA
  cases_ofsp <- read_csv("{opt$b}data/ch/cases_CH.csv") %>% 
    mutate(date = as.Date(date, format = "%m/%d/%Y"))
  cases_data <- select(cases_data, -cases) %>% 
    left_join(cases_ofsp)
}

data <- select(cases_data, 
               date, cases, deaths, hosp_incid, cum_deaths, 
               hosp_curr, icu_curr, discharged, delta_hosp, delta_ID)

# Set start date 
start_date <- with(data, 
                   min(c(date[which(!is.na(cases))[1]] - 5, 
                         date[which(!is.na(hosp_curr))[1]] - 8)))#as.Date("2020-02-20")
start_date <- cases_data$date[1] - 5

end_date <- as.Date("2020-04-08")

# Add rows
data <- rbind(tibble(date = seq.Date(start_date, min(data$date), by = "1 days")) %>% 
                cbind(data[1, -1] %>% mutate_all(function(x) x <- NA)) , data)

data <- data %>% complete(date = seq.Date(start_date,end_date, by="day"))

write_csv(data, glue("{opt$b}COVID-pomp/interm/data_{suffix}.csv"))

# Build pomp object -------------------------------------------------------
source(glue("{opt$b}COVID-pomp/scripts/pomp_skeleton_v2.R"))

# Start and end dates of epidemic
t_start <- dateToYears(start_date)
t_end <- dateToYears(end_date)

geodata <- read_csv(config$setup, col_types = cols())
# initialize empty paramter vector
params <- set_names(rep(0, length(param_names)), param_names)
input_params <- unlist(yaml::read_yaml(config$parameters))
params[param_fixed_names] <- as.numeric(input_params[param_fixed_names])

if (canton == "CH") {
  params["pop"] <- sum(geodata$pop2018, na.rm = T)
} else {
  params["pop"] <- geodata$pop2018[geodata$ShortName == canton]
}

# Initialize the parameters to estimate (just initial guesses)
params["std_W"] <- 0#1e-4
params["std_X"] <- 0#1e-4
params["epsilon"] <- 1
params["k"] <- 5
params["I_0"] <- 10/params["pop"]

# adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days)
params[param_rates_in_days_names] <- params[param_rates_in_days_names] * 365.25
params["R0_0"] <- 2

# rate of simulation in fractions of years
dt_yrs <- 1/(3*365.25)

# data for pomp
data_pomp <- data %>% 
  mutate(time = dateToYears(date)) %>% 
  select(-date, -cum_deaths)

# Transformation of parameters
log_trans <- c("std_X", "R0_0", ifelse(ll_cases, "k", NA)) %>% .[!is.na(.)]
logit_trans <- c("I_0", ifelse(ll_cases, "epsilon", NA)) %>% .[!is.na(.)]

for (par in seq_along(config$parameters_to_fit)) {
  if (config$parameters_to_fit[[par]]$transform == "log"){
    log_trans <- c(log_trans, names(config$parameters_to_fit)[par])
  } else if (config$parameters_to_fit[[par]]$transform == "logit"){
    logit_trans <- c(logit_trans, names(config$parameters_to_fit)[par])
  } 
}

# Time before which we constrain R0
tvary <- dateToYears(as.Date(config$tvary))
# Hospitalized CFR
hcfr <- (1- params["pi2hs"]) * (params["ph2u"] * params["pu2d"] + 1 - params["ph2u"])
sdfrac <- config$sdfrac

covid <- pomp(
  # set data
  data = data_pomp, 
  # time column
  times = "time",
  # initialization time
  t0 = t_start - dt_yrs,
  # parameter names
  paramnames = param_names,
  # parameter vector
  params = params,
  # names of state variables
  statenames = state_names,
  # process simulator
  rprocess = euler(step.fun = proc.Csnippet, delta.t = dt_yrs),
  # measurement model simulator
  rmeasure =  rmeasure.Csnippet,
  # measurement model density
  dmeasure = dmeasure.Csnippet,
  # names of accumulator variables to be re-initalized at each observation timestep 
  accumvars = state_names[str_detect(state_names, "a_")],
  # initializer
  rinit = init.Csnippet,
  partrans = parameter_trans(log = log_trans, logit = logit_trans),
  # Add density and generation fuctions for Skellam dist
  globals = str_c(dskellam.C, rskellam.C, 
                  glue("double tvary = {tvary}; \n
                        double hcfr = {hcfr};\n
                       double sdfrac = {sdfrac};"), sep = "\n")
)

save(covid, file = glue("COVID-pomp/interm/pomp_{suffix}.rda"))

cat("************ \n SAVED:", canton, "pomp object \n************")


cat("************ \nFITTING:", canton, "with likelihood:", opt$likelihood, "\n************")


# files for results
ll_filename <- glue("{opt$b}COVID-pomp/results/loglik_exploration_{suffix}.csv")
mif_filename <- glue("{opt$b}COVID-pomp/results/mif_sir_sde_{suffix}.rda")


# parallel computations
cl <- makeCluster(opt$cores)
registerDoSNOW(cl)

# Setup MIF parameters -----------------------------------------------------

# values of the random walks standard deviations
rw.sd_rp <- 0.02  # for the regular (process and measurement model) parameters
rw.sd_rp_prior <- 0.005  # for the more uncertain random walk on std_X
rw.sd_ivp <- 0.2  # for the initial value parameters
rw.sd_ivp_prior <- 0.05  # for the initial value parameters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_rp_prior, rw.sd_ivp, rw.sd_ivp_prior), c("regular", "prior", "ivp", "ivpprior"))

# lower bound for positive parameter values
min_param_val <- 1e-5 
# define the bounds for the parameters to estimate
parameter_bounds <- tribble(
  ~param, ~lower, ~upper,
  # Process noise
  "std_X", -3, -1, #in log-scale
  # Measurement model
  "k", .1, 10,
  "epsilon", 0.2, 0.5,
  # Initial conditions
  "I_0", 10/params["pop"], 100/params["pop"],
  "R0_0", 1.5, 3
)

if (!is.null(config$parameters_to_fit)) {
  # additional params to fit
  other_bounds <- mapply(x = names(config$parameters_to_fit),
                         y = config$parameters_to_fit, 
                         function(x, y) tibble(param = x, lower = y$lower, upper = y$upper),
                         SIMPLIFY = F) %>% 
    bind_rows()
  parameter_bounds <- rbind(parameter_bounds, other_bounds)
}

if (!ll_cases) {
  parameter_bounds <- parameter_bounds %>% filter(!(param %in% c("k", "epsilon")))
}

# convert to matrix for ease
parameter_bounds <- set_rownames(as.matrix(parameter_bounds[, -1]), parameter_bounds[["param"]])

# create random vectors of initial parameters given the bounds
init_params <- sobolDesign(lower = parameter_bounds[, "lower"],
                           upper = parameter_bounds[, "upper"], 
                           nseq =  sir_Ninit_param[run_level]) 

# convert certain parameters back to log-scale 
param_logbound_names <- c("std_X")
init_params[, param_logbound_names] <- exp(init_params[, param_logbound_names])

# bind with the fixed valued parameters
init_params <- cbind(init_params, 
                     matrix(rep(params[param_fixed_names],
                                each = sir_Ninit_param[run_level]),
                            nrow = sir_Ninit_param[run_level]) %>% 
                       set_colnames(param_fixed_names))



if (!is.null(config$parameters_to_fit)) {
  other_rw <- lapply(names(config$parameters_to_fit), 
                     function(x) glue(", {x} = {rw.sd_param['regular']}")) %>% 
    unlist() %>% 
    str_c(sep = " ")
} else {
  other_rw <- ""
}

rw_text <- glue("rw.sd( std_X  = {rw.sd_param['regular']}
                , k  = {ifelse(ll_cases, rw.sd_param['regular'], 0)}
                , epsilon   = {ifelse(ll_cases, rw.sd_param['regular'], 0)}
                 , I_0  = ivp({rw.sd_param['ivp']})
                 , R0_0  = ivp({rw.sd_param['ivp']})
                {other_rw})")

job_rw.sd <- eval(parse(text = rw_text))

# MIF --------------------------------------------------------------------------

# MIF it! 
t1 <- system.time({
  mf <- foreach(parstart = iter(init_params, by = "row"),
                .inorder = F, 
                .packages = "pomp",
                .errorhandling = "remove") %dopar% {
                  mif2(covid,
                       params = parstart,
                       Np = sir_Np[run_level],
                       Nmif = sir_Nmif[run_level],
                       cooling.type = "geometric",
                       cooling.fraction.50 = findAlpha(sir_Nmif[run_level], nrow(data_pomp)),
                       rw.sd = job_rw.sd,
                       verbose = F)
                }})

save(t1, mf, file = mif_filename)

cat("----- Done MIF, took", round(t1["elapsed"]/60), "mins \n")
# Log-lik ------------------------------------------------------------------

t2 <- system.time({
  liks <- foreach(mfit = mf,
                  .inorder = T, 
                  .packages = "pomp",
                  .combine = rbind,
                  .errorhandling = "remove"
  ) %dopar% {
    # compute log-likelihood estimate by repeating filtering with the given param vect
    ll <-  logmeanexp(
      replicate(sir_Nreps_global[run_level],
                logLik(
                  pfilter(covid,
                          params = pomp::coef(mfit),
                          Np = sir_NpLL[run_level])
                )
      ), se = TRUE)
    # save to dataframe
    data.frame(loglik = ll[1], loglik_se = ll[2], t(coef(mfit)))
  }
})

write_csv(liks, path = ll_filename, append = file.exists(ll_filename))

cat("----- Done LL, took", round(t2["elapsed"]/60), "mins \n")

