# Title: POMP analysis of COVID epidemic in CH
# Description: POMP analysis

# Preamble ---------------------------------------------------------------
library(tidyverse)
library(doSNOW)
library(pomp)
library(magrittr)
library(foreach)
library(itertools)
library(lubridate)
library(sf)
library(glue)
# rm(list = ls())
select <- dplyr::select

source("COVID-pomp/scripts/skellam.R")
source("COVID-pomp/scripts/mifCooling.R")
# Setup ------------------------------------------------------------------------
# Read Rscript arguments
args <- commandArgs(trailingOnly = T)
# args <- c("VD", "d")
canton <- args[1]

# Which likelihood components to use?
# deltah: balances of inputs and outputs from hospitals
# c: cases
# d: total deaths
lik_components <- str_split(args[2], "-")[[1]]#c("deltah", "c", "d")
lik_log <- str_c(str_c("ll_", lik_components), collapse = "+")
lik <- str_c(str_c("ll_", lik_components), collapse = "*")

# Test for cases in the likelihood
ll_cases <- "c" %in% lik_components

cat("************ \nFITITNG:", canton, "with likelihood:", args[2], "\n************")

suffix <- glue("covid_{canton}_{str_c(lik_components, collapse = '-')}")

# files for results
ll_filename <- glue("COVID-pomp/results/loglik_exploration_{suffix}.csv")
mif_filename <- glue("COVID-pomp/results/mif_sir_sde_{suffix}.rda")
filter_filename <- glue("COVID-pomp/results/filtered_{suffix}.rds")

ncpus <- 8

# Level of detail on which to run the computations
run_level <- 3
sir_Np <- c(1e3, 3e3, 4e3)
sir_Nmif <- c(2, 20, 200)
sir_Ninit_param <- c(ncpus, 8, ncpus)
sir_NpLL <- c(1e3, 1e4, 1e4)
sir_Nreps_global <- c(2, 5, 10)

n_filter <- 1e3

# parallel computations
cl <- makeCluster(ncpus)
registerDoSNOW(cl)

# pomp analysis -----------------------------------------------------------

# function to convert dates to fractions of years for model
dateToYears <- function(date, origin = as.Date("2020-01-01"), yr_offset = 2020) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}


# Data ---------------------------------------------------------------

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

end_date <- as.Date("2020-05-31")

# Add rows
data <- rbind(tibble(date = seq.Date(start_date, min(data$date), by = "1 days")) %>% 
                cbind(data[1, -1] %>% mutate_all(function(x) x <- NA)) , data)


# Model specification -----------------------------------------------------

## state variable names:
# S:  Susceptibles
# E:  Exposed
# L:  Latent (non-symptomatic but infectious)
# I:  Infected
# R:  Recovered

## data names:
# hosp: reported number of hospitalized people

## parameter names:
### Pop dynamics
# N: total population
### Infection dynamics (suffix _s|d indicate survival/death outcome)
# l2i: rate L -> I
# i2d: rate I -> D
# i2h: rate I -> H/H_s
# h2hd: rate H -> H_d
# hs2r: rate H_s -> R
# hd2d: rate H_d -> D
# h2u: rate H -> U_s/U_d
# us2r: rate U_s -> R
# ud2d: rate U_d -> D
#### Extra-demographic stochasticity
# tau: generation time
# std_W:    standard deviation of the weiner process to perturb the foi       
### Measurement model
# epsilon:  under-reporting fraction

# Set variables -----------------------------------------------------------

# Number of compartements for each variable to represent Erlang distributions
n_compartments <- list(
  S = 1,
  E = 1,
  L = 1,
  I = 3,
  I_d = 1,
  H_s = 2,
  H = 1,
  U_s = 2,
  H_d = 2,
  U_d = 2,
  R = 1
)

# define stat variable names for each district
state_names <- mapply(
  function(n, comp) {if(n == 1) {comp} else {str_c(comp, 1:n)}},
  n = n_compartments, 
  comp = names(n_compartments)) %>% 
  unlist() %>% 
  c("H_curr", "U_curr", "D", "X", "N",
    "a_I", "a_H", "a_U", "a_D", "a_DH", "a_DI", "a_DU", "a_O", "a_deltaH",
    "Rt", "tot_I") # prefix a_ represent accumulator variables for incidences

# define parameter names for pomp
## process model parameters names to estimate
#param_proc_est_names <- c("std_X", "k", "epsilon")

## initial value parameters to estimate
param_iv_est_names <- c("I_0", "R0_0")

## fixed process model parameters 
rate_names <- c("e2i", "l2i", "id2o", "i2h", "i2o", "hs2r", "hd2d", "h2u", "us2r", "ud2d")
prob_names <- c("psevere", "pi2d", "pi2h", "pi2hs", "ph2u", "pu2d")
#param_proc_fixed_names <- c("pop", rate_names, prob_names,  "alpha", "std_W")
# define parameter names for pomp
## process model parameters names to estimate
if (ll_cases) {
  param_proc_est_names <- c("std_X", "k", "epsilon")
  param_proc_fixed_names <- c("pop", rate_names, prob_names, "std_W")
} else {
  param_proc_est_names <- c("std_X")
  param_proc_fixed_names <- c("pop", rate_names, prob_names, "std_W", "k", "epsilon")
}

# all parameter names to estimate
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters
param_fixed_names <- c(param_proc_fixed_names)

# all param names
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates
param_rates_in_days_names <- rate_names

# Measurment model  -------------------------------------------------------

# measurement model
dmeasure.Csnippet <- Csnippet(glue("
  double ll_na;
  double ll_c, ll_h, ll_hc, ll_d, ll_deltah, ll_dh, ll_du, ll_di;
  ll_na = (give_log) ? 0 : 1;
  if (!ISNA(cases)) {
    ll_c = dnbinom_mu(cases, k, epsilon * a_I, give_log);
  } else {
    ll_c = ll_na;
  }
  if (!ISNA(delta_ID)) {
    ll_deltah = dskellam(delta_ID, a_H, a_DH + a_DU, give_log);
  } else {
      if (!ISNA(delta_hosp)) { 
      // if info on releases not available, give hosp curr
      ll_deltah = dskellam(delta_hosp, a_H, a_DH + a_DU + a_O, give_log);
    } else {
      ll_deltah = ll_na;
    }
  }
  if (!ISNA(hosp_curr)) {
    ll_hc = dpois(hosp_curr, H_curr, give_log);
  } else {
    ll_hc = ll_na;
  }
  if (!ISNA(hosp_incid)) {
    ll_h = dpois(hosp_incid, a_H, give_log);
  } else {
    ll_h = ll_na;
  }
  if (!ISNA(deaths)) {
   ll_d = dpois(deaths, a_D, give_log);
  } else {
    ll_d = ll_na;
  }
  
  if (give_log) {
    lik =  {{lik_log}};
  } else {
    lik =  {{lik}};
  }
                              ", .open = "{{", .close = "}}"))

## NegBinomial simulator
rmeasure.Csnippet <- Csnippet("
                      int deltah = rskellam(a_H, a_DH + a_DU);
                      deaths = rpois(a_D);
                      hosp_curr = rpois(H_curr);
                      hosp_incid = rpois(a_H);
                      cases = rnbinom_mu(k, epsilon * a_I);
                      delta_hosp = deltah;
                      delta_ID = deltah + rpois(a_O);
                      ")

# Process model -----------------------------------------------------------------


# create C code for each district
proc.Csnippet <- Csnippet("
                          double foi, foi_stoc; // force of infection and its stochastic version
                          double dw, dWX;            // extra-demographic stochasticity on foi
                          double rate[21];      // vector of all rates in model
                          double dN[21];        // vector of transitions between classes during integration timestep

                          // force of infection
                          //foi = (t <= tlockdown) ? exp(X) * (I1)/N : alpha * exp(X) * (I1)/N;
                          foi =  exp(X) * (I1 + I2 + I3)/N;

                          if(std_W > 0.0) {
                          // white noise (extra-demographic stochasticity)
                            dw = rgammawn(std_W, dt);
                          // apply stochasticity
                            foi_stoc = foi * dw/dt;
                          } else {
                            foi_stoc = foi;
                          }
                          
                          // define transition rates for each type of event
                          // S compartment
                          rate[0] = foi_stoc;   // infections
                          // E compartment
                          rate[1] = e2i;    // transition to L
                          // L compartment
                          //rate[4] = l2i;   // become symptomatic
                          // I1 compartment
                          rate[2] = i2o * 3;  // transition to I2
                           // I2 compartment
                          rate[3] = i2o * 3;  // transition to I3
                          // I3 compartment
                          rate[4] = psevere * (1-pi2h) * i2o * 3;    // deaths without hospitalization
                          rate[5] = psevere * pi2h * pi2hs * i2o * 3;   // hospitalization that WILL NOT go in ICU nor die
                          rate[6] = psevere * pi2h * (1-pi2hs) * (1-ph2u) * i2o * 3;   // hospitalization that WILL NOT go in ICU and die
                          rate[7] = psevere * pi2h * (1-pi2hs) * ph2u * i2o * 3;   // hospitalization that WILL go in ICU or die
                          rate[8] = (1-psevere) * i2o * 3;      // recovery from infection
                          
                          // I_d compartment
                          rate[9] = (1-pi2d) * id2o; // severe infected that recover
                          rate[10] = pi2d * id2o;    // severe infected that die
                          
                          // H_s1 compartment
                          rate[11] = hs2r;   // transition to H_S2
                          // H_s2 compartment
                          rate[12] = hs2r;  // recover
                          
                          // H compartment
                          rate[13] = h2u * pu2d; // ICU admission that die
                          rate[14] = h2u * (1-pu2d); // ICU admission that survive
                          
                          // H_d1 compartment
                          rate[15] = hd2d; // transition to H_d2
                          // H_d2 compartment
                          rate[16] = hd2d; // death
                          
                          // U_s1 compartment
                          rate[17] = us2r;    // recovery
                          // U_s2 compartment
                          rate[18] = us2r;    // recovery
                          // U_d1 compartment
                          rate[19] = ud2d;   // transition to U_d2
                          // U_d2 compartment
                          rate[20] = ud2d;    // death
                          
                          // simulate all transitions
                          reulermultinom(1, S,  &rate[0], dt, &dN[0]);
                          reulermultinom(1, E, &rate[1], dt, &dN[1]);
                          //reulermultinom(1, L,  &rate[4], dt, &dN[4]);
                          reulermultinom(1, I1, &rate[2], dt, &dN[2]);
                          reulermultinom(1, I2, &rate[3], dt, &dN[3]);
                          reulermultinom(5, I3, &rate[4], dt, &dN[4]);
                          reulermultinom(2, I_d, &rate[9], dt, &dN[9]);
                          reulermultinom(1, H_s1, &rate[11], dt, &dN[11]);
                          reulermultinom(1, H_s2, &rate[12], dt, &dN[12]);
                          reulermultinom(2, H,    &rate[13], dt, &dN[13]);
                          reulermultinom(1, H_d1, &rate[15], dt, &dN[15]);
                          reulermultinom(1, H_d2, &rate[16], dt, &dN[16]);
                          reulermultinom(1, U_s1,  &rate[17], dt, &dN[17]);
                          reulermultinom(1, U_s2,  &rate[18], dt, &dN[18]);
                          reulermultinom(1, U_d1, &rate[19], dt, &dN[19]);
                          reulermultinom(1, U_d2, &rate[20], dt, &dN[20]);

                          // update state variables
                          S    += -dN[0];
                          E    += dN[0] - dN[1];
                          //L    += dN[3] - dN[4];
                          I1   += dN[1] - dN[2];
                          I2   += dN[2] - dN[3];
                          I3   += dN[3] - dN[4] - dN[5] - dN[6] - dN[7] - dN[8];
                          I_d  += dN[4] - dN[9] - dN[10];
                          H_s1 += dN[5] - dN[11];
                          H_s2 += dN[11] - dN[12];
                          H += dN[7] - dN[13] - dN[14];
                          H_d1 += dN[6] - dN[15];
                          H_d2 += dN[15] - dN[16];
                          U_s1 += dN[14] - dN[17];
                          U_s2 += dN[17] - dN[18];
                          U_d1 += dN[13] - dN[19];
                          U_d2 += dN[19] - dN[20];
                          R += dN[8] + dN[9] + dN[12] + dN[18];
                          D += dN[10] + dN[16] + dN[20];
                          
                          // Accumulators
                          a_I += dN[1];
                          a_H += dN[5] + dN[6] + dN[7];
                          a_U += dN[13] + dN[14];
                          a_DI += dN[10];
                          a_DH += dN[16];
                          a_DU += dN[20];
                          a_D  = a_DH + a_DI + a_DU;
                          a_O  += dN[12] + dN[18]; // discharged from hospital
                          a_deltaH = a_H - a_DH - a_DU - a_O;
                          // Current
                          U_curr = U_s1 + U_s2 + U_d1 + U_d2;
                          H_curr = H_s1 + H_s2 + H + H_d1 + H_d2 + U_s1 + U_s2 + U_d1 + U_d2;
                          // Total infected
                          tot_I += dN[1];

                          // random walk of beta
                          dWX = rnorm(0, sqrt(dt));
                          X  +=  std_X * dWX;  // 
                          
                          
                          
                          // susceptibles so as to match total population
                          N = pop - D - H_curr - U_curr;
                          Rt = exp(X) / (i2o); 
                          ")

# Initializer 
init.Csnippet <- Csnippet("X = log(R0_0 * i2o);
                          E   = 0;
                          L   = 0;
                          I1  =  nearbyint(I_0 * pop);
                          I2  = 0;
                          I3 = 0;
                          I_d = 0;
                          H_s1 = 0;
                          H_s2 = 0;
                          H = 0;
                          U_s1 = 0;
                          U_s2 = 0;
                          H_d1 = 0;
                          H_d2 = 0;
                          U_d1 = 0;
                          U_d2 = 0;
                          H_curr = 0;
                          U_curr = 0;
                          D = 0;
                          a_I = 0;
                          a_H = 0;
                          a_U = 0;
                          a_DU = 0;
                          a_DH = 0;
                          a_DI = 0;
                          a_D = 0;
                          a_O = 0;
                          a_deltaH = 0;
                          R   = 0;
                          S   = nearbyint(pop - I1);
                          N   = pop;
                          tot_I = 0;")

# Build pomp object -------------------------------------------------------

# Start and end dates of epidemic
t_start <- dateToYears(start_date)
t_end <- dateToYears(end_date)

geodata <- read_csv("data/ch/geodata.csv", col_types = cols())
# initialize empty paramter vector
params <- set_names(rep(0, length(param_names)), param_names)
input_params <- unlist(yaml::read_yaml("COVID-pomp/data/parameters.yaml"))
params[param_fixed_names] <- as.numeric(input_params[param_fixed_names])

params["pop"] <- geodata$pop2018[geodata$ShortName == canton]

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
  partrans = parameter_trans(
    log = c("std_X", "R0_0", ifelse(ll_cases, "k", NA)) %>% .[!is.na(.)],
    logit = c("I_0", ifelse(ll_cases, "epsilon", NA)) %>% .[!is.na(.)]
  ),
  # Add density and generation fuctions for Skellam dist
  globals = str_c(dskellam.C, rskellam.C, sep = "\n")
)

save(covid, file = glue("COVID-pomp/interm/pomp_object_{canton}.rda"))

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
  "std_X", -2, 2, #in log-scale
  # Measurement model
  "k", .1, 10,
  "epsilon", 0.2, 0.5,
  # Initial conditions
  "I_0", 10/params["pop"], 100/params["pop"],
  "R0_0", 1.5, 3
)

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

job_rw.sd <- eval(
  parse(
    text = str_c("rw.sd(",
                 " std_X  = ", rw.sd_param["regular"],
                 ", k  = ", ifelse(ll_cases, rw.sd_param["regular"], 0),
                 ", epsilon  = ", ifelse(ll_cases, rw.sd_param["regular"], 0),
                 ", I_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ", R0_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ")")
  )
)

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

ggsave(p, filename = glue("COVID-pomp/results/figs/plot_{suffix}.png"), width = 9, height = 6)
