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

rm(list = ls())

select <- dplyr::select


# Setup ------------------------------------------------------------------------
suffix <- "covidVD_noalpha_noW"
ncpus <- 8

# files for results
ll_filename <- str_c("results/loglik_exploration_", suffix, ".csv")
mif_filename <- str_c("results/mif_sir_sde_", suffix, ".rda")
filter_filename <- str_c("results/filtered_", suffix, ".rda")

# Level of detail on which to run the computations
run_level <- 3
sir_Np <- c(1e3, 3e3, 3e3)
sir_Nmif <- c(2, 20, 100)
sir_Ninit_param <- c(ncpus, 8, ncpus)
sir_NpLL <- c(1e3, 1e4, 1e4)
sir_Nreps_global <- c(2, 5, 10)

# parallel computations
cl <- makeCluster(ncpus)
registerDoSNOW(cl)

# pomp analysis -----------------------------------------------------------

Sys.setlocale("LC_ALL","C")

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


# load data ---------------------------------------------------------------
start_date <- as.Date("2020-02-20")
end_date <- as.Date("2020-05-31")

cases_data <- read_csv("../covid_19/fallzahlen_kanton_total_csv/COVID19_Fallzahlen_Kanton_VD_total.csv") %>% 
  mutate(cases = c(NA, diff(ncumul_conf)),
         deaths = c(NA, diff(ncumul_deceased)))

# load cases from the March-April epidemic
hosp_data <- read_csv("data/VD_hosp_data.csv")

data <- rbind(select(cases_data, date, cases, deaths)) %>% 
  full_join(hosp_data) %>% 
  mutate(deaths_nohosp_incid = pmax(0, deaths - deaths_icu_incid - deaths_noicu_incid))

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
### Infection dynamics (suffix _s|d indicate survival/death ouctome)
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
# epsilon:  under-rerpoting fraction

# Set variables -----------------------------------------------------------

# Number of compartements for each variable to represent Erlang distributions
n_comparements <- list(
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
  n = n_comparements, 
  comp = names(n_comparements)) %>% 
  unlist() %>% 
  c("H_curr", "U_curr", "D", "X", "N", "a_I", "a_H", "a_U", "a_DH", "a_DI", "a_DU", "Rt", "tot_I", "lockflag") # prefix a_ represent accumulator variables for incidences

# define parameter names for pomp
## process model parameters names to estimate
param_proc_est_names <- c("std_X", "k", "epsilon")

## initial value parameters to estimate
param_iv_est_names <- c("I_0", "R0_0")

## fixed process model parameters 
rate_names <- c("e2i", "l2i", "id2o", "i2h", "i2o", "hs2r", "hd2d", "h2u", "us2r", "ud2d")
prob_names <- c("psevere", "pi2d", "pi2h", "pi2hs", "ph2u", "pu2d")
param_proc_fixed_names <- c("pop", rate_names, prob_names,  "alpha", "std_W")

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
## density

## NegBinomial density (if k -> inf then becomes Poisson)
# sum over districts
dmeasure.Csnippet <- Csnippet(" double ll_c, ll_h, ll_dh, ll_du, ll_di;
  if (!ISNA(cases)) {
    ll_c = dnbinom_mu(cases, k, epsilon * a_I, give_log);
  }
  if (!ISNA(hosp_incid)) {
    ll_h = dpois(hosp_incid, a_H, give_log);
  }
  if(!ISNA(deaths_noicu_incid)) {
   //ll_dh = dpois(deaths_noicu_incid, a_DH, give_log);
  }
  if(!ISNA(deaths_icu_incid)) {
   //ll_du = dpois(deaths_icu_incid, a_DU, give_log);
  }
  if(!ISNA(deaths_nohosp_incid)) {
   //ll_di = dpois(deaths_nohosp_incid, a_DI, give_log);
  }
  
  if (give_log) {
    lik = ll_c + ll_h;
  } else {
    lik = ll_c * ll_h;
  }
                              ")


## NegBinomial simulator
rmeasure.Csnippet <- Csnippet("
                      //hosp = rnbinom_mu(k, mean_hosp);
                      hosp_incid = rpois(a_H);
                      deaths_noicu_incid = rpois(a_DH);
                      deaths_icu_incid = rpois(a_DU);
                      deaths_nohosp_incid = rpois(a_DI);
                      cases = rnbinom_mu(k, epsilon * a_I);
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
                          H    += dN[7] - dN[13] - dN[14];
                          H_d1 += dN[6] - dN[15];
                          H_d2 += dN[15] - dN[16];
                          U_s1 += dN[14] - dN[17];
                          U_s2 += dN[17] - dN[18];
                          U_d1 += dN[13] - dN[19];
                          U_d2 += dN[19] - dN[20];
                          R    += dN[8] + dN[9] + dN[12] + dN[18];
                          D    += dN[10] + dN[16] + dN[20];
                          
                          // Accumulators
                          a_I += dN[1];
                          a_H += dN[5] + dN[6] + dN[7];
                          a_U += dN[13] + dN[14];
                          a_DI += dN[6];
                          a_DH += dN[16];
                          a_DU += dN[20];
                          // Current
                          U_curr = U_s1 + U_s2 + U_d1 + U_d2;
                          H_curr = H_s1 + H_s2 + H + H_d1 + H_d2 + U_s1 + U_s2 + U_d1 + U_d2;
                          // Total infected
                          tot_I += dN[1];

                          // random walk of beta
                          dWX = rnorm(0, sqrt(dt));
                          X  +=  std_X * dWX;  // 
                          
                          if (t >= tlockdown & lockflag == 0) {
                            X  +=  log(alpha);  // 
                            lockflag = 1;
                          }
                          
                          
                          // susceptibles so as to match total population
                          N = pop - D - H_curr - U_curr;
                          //Rt = (t<tlockdown) ? exp(X) / (i2r) :  alpha * exp(X) / (i2r); 
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
                          R   = 0;
                          S   = nearbyint(pop - I1);
                          N   = pop;
                          tot_I = 0;
                          lockflag = 0;")

# Build pomp object -------------------------------------------------------

# Start and end dates of epidemic
t_start <- dateToYears(start_date)
t_end <- dateToYears(end_date)

# initialize empty paramter vector
params <- set_names(rep(0, length(param_names)), param_names)
input_params <- unlist(yaml::read_yaml("scripts/parameter.yaml"))
params[param_fixed_names] <- as.numeric(input_params[param_fixed_names])

# Initialize the parameters to estimate (just initial guesses)
params["std_W"] <- 0#1e-4
params["std_X"] <- 0#1e-4
params["epsilon"] <- .3
params["k"] <- 5
params["I_0"] <- 10/params["pop"]
params["alpha"] <- 1

tlockdown <- dateToYears(as.Date(input_params["date_lockdown"]))
# R0 <- exp(params["X_0"])/params["e2l"] 

# adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days)
params[param_rates_in_days_names] <- params[param_rates_in_days_names] * 365.25
# R0 <- 1
params["R0_0"] <- 2

# rate of simulation in fractions of years
dt_yrs <- 1/(3*365.25)

# data for pomp
data_pomp <- data %>% 
  mutate(time = dateToYears(date)) %>% 
  select(time, hosp_incid, cases, deaths_icu_incid, deaths_noicu_incid, deaths_nohosp_incid) 

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
  # (C for cases, W for the white noise just for plotting)
  accumvars = state_names[str_detect(state_names, "a_")],
  # initializer
  rinit = init.Csnippet,
  partrans = parameter_trans(
    log = c("std_W", "std_X", "R0_0", "k"),
    logit = c("I_0", "epsilon", "alpha")
  ),
  globals = glue::glue("double tlockdown = {tlockdown};")
)

save(covid, file = "interm/pomp_object.rda")

sim <- simulate(covid, 10, format = "data.frame", include.data = F) 

sim %>%
  group_by(.id) %>%
  # slice(1:12) %>%
  gather(variable, value, -time, -.id) %>%
  ggplot(aes(x = yearsToDate(time), y = value, color = .id)) +
  geom_line() +
  facet_wrap(~variable, scales = "free") +
  theme_bw()


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
  # "std_W",  -6, -3,  # in log-scale
  # Measurement model
  "k", .1, 1,
  "epsilon", 0.2, 0.5,
  # "alpha", 0.1, .5,
  # Initial conditions
  "I_0", 10/params["pop"], 100/params["pop"],
  "R0_0", 1.5, 3
)

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
                 # ",std_W  = ", rw.sd_param["regular"],
                 " std_X  = ", rw.sd_param["regular"],
                 ", k  = ", rw.sd_param["regular"],
                 ", epsilon  = ", rw.sd_param["regular"],
                 # ", alpha  = ", rw.sd_param["regular"],
                 ", I_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ", R0_0  = ivp(",  rw.sd_param["ivp"], ")",
                 ")")
  )
)

# MIF --------------------------------------------------------------------------

# MIF it! 
# parstart <- init_params[4, ]

t1 <- system.time({
  mf <- foreach(parstart = iter(init_params, by = "row"),
                .inorder = F, 
                .packages = "pomp2",
                .errorhandling = "remove") %dopar% {
                  # system.time({
                  mif2(covid,
                       params = parstart,
                       Np = sir_Np[run_level],
                       Nmif = sir_Nmif[run_level],
                       cooling.type = "geometric",
                       cooling.fraction.50 = 0.6,
                       rw.sd = job_rw.sd,
                       verbose = F)
                  # })
                }})

save(mf, file = mif_filename)

# Log-lik ------------------------------------------------------------------

t2 <- system.time({
  liks <- foreach(mfit = mf,
                  .inorder = T, 
                  # i = icount(1),
                  .packages = "pomp2",
                  .combine = rbind,
                  .errorhandling = "remove"
  ) %dopar% {
    # compute log-likelihood estimate by repeating filtering with the given param vect
    # system.time({
    ll <-  logmeanexp(
      replicate(sir_Nreps_global[run_level],
                logLik(
                  pfilter(covid,
                          params = pomp2::coef(mfit),
                          Np = sir_NpLL[run_level])
                )
      ), se = TRUE)
    # })
    # save to dataframe
    data.frame(loglik = ll[1], loglik_se = ll[2], t(coef(mfit)))
  }
})

write_csv(liks, path = ll_filename, append = file.exists(ll_filename))


# Filter --------------------------------------------------------------------
liks <- read_csv(ll_filename)

best_params <- liks %>%
  arrange(desc(loglik)) %>%
  slice(1:3) %>% 
  select(-contains("log"))

n_filter <- 500

filter_dists <- foreach(par = iter(best_params, "row"),
                        pit = icount(nrow(best_params)),
                        .packages = c("pomp2", "tidyverse", "foreach", "magrittr"),
                        .combine = rbind) %do% {
                          foreach(it = icount(n_filter),
                                  .combine = rbind,
                                  .packages = c("pomp2", "tidyverse", "foreach", "magrittr")
                          ) %dopar% {
                            
                            pf <- pfilter(covid, params = as.vector(par), Np = 3e3, filter.traj = T)
                            
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
save(filter_dists, file = filter_filename)

stopCluster(cl)
closeAllConnections()


# Plot -------------------------------------------------------------------------

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


p <- ggplot(filter_stats, aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
  # geom_line(aes(y = mean)) +
  geom_point(data = data %>%
               mutate(time = dateToYears(date),
                      cases = cases / best_params[["epsilon"]]) %>%
               rename(a_I = cases,
                      a_DU = deaths_icu_incid,
                      a_DH = deaths_noicu_incid,
                      a_DI = deaths_nohosp_incid,
                      a_H = hosp_incid,
                      a_U = icu_incid,
                      U_curr = icu_curr,
                      H_curr = hosp_curr
               ) %>% 
               gather(var, value, -time, -date),
             aes(y = value)) +
  facet_wrap(~var, scales = "free")  #+

p
ggsave(p, filename = str_c("results/plot_",  suffix, ".png"), width = 10, height = 10)

