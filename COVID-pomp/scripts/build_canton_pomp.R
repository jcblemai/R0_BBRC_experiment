# Title: POMP model of COVID19 in CH, placeal level
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
library(parallel)
library(glue)
library(optparse)

select <- dplyr::select

option_list <- list(
  make_option(c("-c", "--config"), default = "pomp_config.yaml", type = "character", help = "path to the config file"),
  make_option(c("-a", "--asindex"), default = 0, type = "numeric", help = "whether to use the index of a slurm array"),
  make_option(c("-b", "--basepath"), default = "COVID-pomp/", type = "character", help = "base path"),
  make_option(c("-j", "--jobs"), default = detectCores(), type = "numeric", help = "number of cores used"),
  make_option(c("-o", "--cores"), default = detectCores(), type = "numeric", help = "number of cores used"),
  make_option(c("-r", "--run_level"), default = 1, type = "numeric", help = "run level for MIF"),
  make_option(c("-p", "--place"), default = "CH", type = "character", help = "name of place to be run, a place abbrv. in CH"),
  make_option(c("-l", "--likelihood"), default = "d-deltah", type = "character", help = "likelihood to be used for filtering"),
  make_option(c("-s", "--suffix"), default = NULL, type = "character", help = "custom suffix to add")
  )

opt <- parse_args(OptionParser(option_list = option_list))
config <- yaml::read_yaml(opt$config)

source(glue("{opt$b}scripts/skellam.R"))
source(glue("{opt$b}scripts/mifCooling.R"))
source(glue("{opt$b}scripts/utils.R"))
source(glue("{opt$b}scripts/utils_to_custom.R"))

if (opt$a == 1 & Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
  array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  opt$place <- config$places[array_id]
}

# Level of detail on which to run the computations
run_level <- opt$run_level
sir_Np <- c(1e2, 3e3, 5e3)
sir_Nmif <- c(2, 20, 150)
sir_Ninit_param <- c(2, opt$jobs, opt$jobs)
sir_NpLL <- c(1e2, 1e4, 1e4)
sir_Nreps_global <- c(2, 5, 20)

simtest <- F  # flag to test simulations from pomp object

# parallel computations
cl <- makeCluster(opt$cores)
registerDoSNOW(cl)

# Setup ------------------------------------------------------------------------
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

# Data ---------------------------------------------------------------
data_file <- grep(place, list.files(config$data_path, full.names = T), value = T)

epidata <- ingestData(fdata = data_file,
                      likelihood = parsed_lik,
                      dt = "day",
                      params = list(place = place),
                      endate = as.Date('2020-04-24'))

# Write data to csv for reproducibility
write_csv(epidata, glue("{opt$b}interm/data_{suffix}.csv"))

# Build pomp object -------------------------------------------------------
source(glue("{opt$b}scripts/{config$model}"))

# Start and end dates of epidemic
t_start <- dateToYears(min(epidata$date))
t_end <- dateToYears(max(epidata$date))

# Set parameters
params <- setParameters(config, opt)

# rate of simulation in fractions of years
dt_yrs <- 1 / (3 * 365.25)

# Transformation of parameters
log_trans <- c("std_X", ifelse(use_case_incid, "k", NA)) %>% .[!is.na(.)]
logit_trans <- c(ifelse(use_case_incid, "epsilon", NA)) %>% .[!is.na(.)]

for (par in seq_along(config$parameters_to_fit)) {
  if (config$parameters_to_fit[[par]]$transform == "log") {
    log_trans <- c(log_trans, names(config$parameters_to_fit)[par])
  } else if (config$parameters_to_fit[[par]]$transform == "logit") {
    logit_trans <- c(logit_trans, names(config$parameters_to_fit)[par])
  }
}

# Time before which we constrain R0
tvary <- dateToYears(as.Date(config$tvary))
# Hospitalized CFR
sdfrac <- config$sdfrac
globals <-  glue("double tvary = {tvary}; 
                 double sdfrac = {sdfrac};")
covid <- pomp(
  data = select(epidata, -date),  # set data
  times = "time",  # time column
  t0 = t_start - dt_yrs,     # initialization time
  paramnames = param_names,  # parameter names
  params = params,           # parameter vector
  statenames = state_names,  # names of state variables
  rprocess = euler(step.fun = proc.Csnippet, delta.t = dt_yrs),   # process simulator
  rmeasure = rmeasure.Csnippet,  # measurement model simulator
  dmeasure = dmeasure.Csnippet,  # measurement model density
  # names of accumulator variables to be re-initalized at each observation timestep
  accumvars = state_names[str_detect(state_names, "a_")],
  rinit = init.Csnippet,   # initializer
  partrans = parameter_trans(log = log_trans, logit = logit_trans),
  # Add density and generation fuctions for Skellam dist
  globals = str_c(dskellam.C, rskellam.C, globals, sep = "\n")
)

if (simtest) {
  simulate(covid, nsim = 10, format = "data.frame") %>%
    gather(var, value, -time, -.id) %>%
    filter(time <= dateToYears(as.Date(c("2020-03-20")))) %>%
    ggplot(aes(x = time, y = value, color = .id, group = .id)) +
    geom_line() +
    scale_color_viridis_d() +
    facet_wrap(~var, scales = "free_y")
}

# Save pomp object 
save(covid, file = glue("{opt$b}interm/pomp_{suffix}.rda"))


# Setup MIF parameters -----------------------------------------------------
cat("************ \nFITTING:", place, "with likelihood:", opt$likelihood, "\n************")

# files for results
ll_filename <- glue("{opt$b}results/loglik_exploration_{suffix}.csv")
mif_filename <- glue("{opt$b}results/mif_sir_sde_{suffix}.rda")

# values of the random walks standard deviations
rw.sd_rp <- 0.02 # for the regular (process and measurement model) parameters
rw.sd_ivp <- 0.2 # for the initial value parameters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_ivp), c("regular", "ivp"))

init_params <- getInitParams(config = config, 
                             params = params,
                             param_fixed_names = param_fixed_names,
                             use_case_incid = use_case_incid, 
                             npar = sir_Ninit_param[run_level])

job_rw.sd <- getRWSD(config = config, rw.sd_param = rw.sd_param, tvary = tvary)

# MIF --------------------------------------------------------------------------

# MIF it!
t1 <- system.time({
  mf <- foreach(
    parstart = iter(init_params, by = "row"),
    .inorder = F,
    .packages = "pomp",
    .errorhandling = "stop"
  ) %dopar% {
    mif2(covid,
         params = parstart,
         Np = sir_Np[run_level],
         Nmif = sir_Nmif[run_level],
         cooling.type = "geometric",
         cooling.fraction.50 = findAlpha(sir_Nmif[run_level], nrow(epidata), 0.05),
         rw.sd = job_rw.sd,
         verbose = T
    )
  }
})

save(t1, mf, file = mif_filename)

cat("----- Done MIF, took", round(t1["elapsed"] / 60), "mins \n")

# Log-lik ------------------------------------------------------------------

t2 <- system.time({
  liks <- foreach(
    mfit = mf,
    .inorder = T,
    .packages = "pomp",
    .combine = rbind,
    .errorhandling = "remove"
  ) %dopar% {
    # compute log-likelihood estimate by repeating filtering with the given param vect
    ll <- logmeanexp(
      replicate(
        sir_Nreps_global[run_level],
        logLik(
          pfilter(covid,
                  params = pomp::coef(mfit),
                  Np = sir_NpLL[run_level]
          )
        )
      ), se = TRUE
    )
    # save to dataframe
    data.frame(loglik = ll[1], loglik_se = ll[2], t(coef(mfit)))
  }
})

write_csv(liks, path = ll_filename, append = file.exists(ll_filename))

cat("----- Done LL, took", round(t2["elapsed"] / 60), "mins \n")

stopCluster(cl)
