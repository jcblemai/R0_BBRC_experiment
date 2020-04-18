
# Set variables -----------------------------------------------------------

# Number of compartements for each variable to represent Erlang distributions
n_compartments <- list(
  S = 1,
  E = 1,
  I = 3,
  I_d = 1,
  I_h = 1,
  H_sa = 1,
  H_sb = 2,
  H = 1,
  U_sa = 1,
  U_sb = 1,
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

## initial value parameters to estimate
param_iv_est_names <- c("I_0", "R0_0")

## fixed process model parameters 
rate_names <- c("e2i", "i2h", "id2o", "i2o", "hs2ra", "hs2rb", "hd2d", "h2u", "us2ra", "us2rb", "ud2d")
prob_names <- c("psevere", "pi2d", "pi2h", "pi2hs", "ph2u", "pu2d", "pusa", "phsa")
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

if (!is.null(config$parameters_to_fit)) {
  allin <- unlist(lapply(names(config$parameters_to_fit), function(x) !(x %in% c(rate_names, prob_names))))
  if (sum(allin) == 0) {
  param_proc_est_names <- c(param_proc_est_names, names(config$parameters_to_fit))
  param_proc_fixed_names <- param_proc_fixed_names[!(param_proc_fixed_names %in% names(config$parameters_to_fit))]
  }
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
    if ({{downweight}} == 1) {
      ll_c = (give_log) ? ll_c + log(0.5) : ll_c * 0.5;
    }
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
  if (!ISNA(delta_ID)) {
    ll_h = dskellam(delta_ID, a_H, a_DH + a_DU, give_log);
  } else {
  if (!ISNA(delta_hosp)) { 
    // if info on releases not available, give hosp curr
    ll_h = dskellam(delta_hosp, a_H, a_DH + a_DU + a_O, give_log);
  } else {
    ll_h = ll_na;
  }
  }
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
                          double dw, dWX;       // extra-demographic stochasticity on foi
                          double rate[25];      // vector of all rates in model
                          double dN[25];        // vector of transitions between classes during integration timestep

                          // force of infection
                          foi =  exp(X) * (I1 + I2 + I3)/N;
                          
                          // define transition rates for each type of event
                          // S compartment
                          rate[0] = foi;   // infections
                          // E compartment
                          rate[1] = e2i;    // transition to L
                          // I1 compartment
                          rate[2] = i2o * 3;  // transition to I2
                           // I2 compartment
                          rate[3] = i2o * 3;  // transition to I3
                          // I3 compartment
                          rate[4] = psevere * (1-pi2h) * i2o * 3;    // deaths without hospitalization
                          rate[5] = psevere * pi2h * i2o * 3;    // hospitalizations
                          rate[6] = (1-psevere) * i2o * 3;      // recovery from infection
                          
                          // I_d compartment
                          rate[7] = (1-pi2d) * id2o; // severe infected that recover
                          rate[8] = pi2d * id2o;    // severe infected that die
                          
                          // I_h compartment
                          rate[9] =  i2h * pi2hs * i2o * phsa;   // hospitalization that WILL NOT go in ICU nor die
                          rate[10] = i2h * pi2hs * i2o * (1-phsa);   // hospitalization that WILL NOT go in ICU nor die
                          rate[11] = i2h * (1-pi2hs) * (1-ph2u) * i2o;   // hospitalization that WILL NOT go in ICU and die
                          rate[12] = i2h * (1-pi2hs) * ph2u * i2o;   // hospitalization that WILL go in ICU or die
                          
                          // H_sa compartment
                          rate[13] = hs2ra;   // recover
                          // H_sb1 compartment
                          rate[14] = hs2rb;   // transition to H_Sb2
                          // H_sb2 compartment
                          rate[15] = hs2rb;  // recover
                          
                          // H compartment
                          rate[16] = h2u * pu2d; // ICU admission that die
                          rate[17] = h2u * (1-pu2d) * pusa; // ICU admission that survive
                          rate[18] = h2u * (1-pu2d) * (1-pusa); // ICU admission that survive
                          
                          // H_d1 compartment
                          rate[19] = hd2d; // transition to H_d2
                          // H_d2 compartment
                          rate[20] = hd2d; // death
                          
                          // U_sa compartment
                          rate[21] = us2ra;    // recovery
                          // U_sb compartment
                          rate[22] = us2rb;    // recovery
                          
                          // U_d1 compartment
                          rate[23] = ud2d;   // transition to U_d2
                          // U_d2 compartment
                          rate[24] = ud2d;    // death
                          
                          // simulate all transitions
                          reulermultinom(1, S,  &rate[0], dt, &dN[0]);
                          reulermultinom(1, E, &rate[1], dt, &dN[1]);
                          reulermultinom(1, I1, &rate[2], dt, &dN[2]);
                          reulermultinom(1, I2, &rate[3], dt, &dN[3]);
                          reulermultinom(3, I3, &rate[4], dt, &dN[4]);
                          reulermultinom(2, I_d, &rate[7], dt, &dN[7]);
                          reulermultinom(4, I_h, &rate[9], dt, &dN[9]);
                          reulermultinom(1, H_sa, &rate[13], dt, &dN[13]);
                          reulermultinom(1, H_sb1, &rate[14], dt, &dN[14]);
                          reulermultinom(1, H_sb2, &rate[15], dt, &dN[15]);
                          reulermultinom(3, H,    &rate[16], dt, &dN[16]);
                          reulermultinom(1, H_d1, &rate[19], dt, &dN[19]);
                          reulermultinom(1, H_d2, &rate[20], dt, &dN[20]);
                          reulermultinom(1, U_sa,  &rate[21], dt, &dN[21]);
                          reulermultinom(1, U_sb,  &rate[22], dt, &dN[22]);
                          reulermultinom(1, U_d1, &rate[23], dt, &dN[23]);
                          reulermultinom(1, U_d2, &rate[24], dt, &dN[24]);

                          // update state variables
                          S    += -dN[0];
                          E    += dN[0] - dN[1];
                          I1   += dN[1] - dN[2];
                          I2   += dN[2] - dN[3];
                          I3   += dN[3] - dN[4] - dN[5] - dN[6];
                          I_d  += dN[4] - dN[7] - dN[8];
                          I_h  += dN[5] - dN[9] - dN[10] - dN[11] - dN[12];
                          H_sa += dN[9] - dN[13];
                          H_sb1 += dN[10] - dN[14];
                          H_sb2 += dN[14] - dN[15];
                          H += dN[12] - dN[16] - dN[17] - dN[18];
                          H_d1 += dN[11] - dN[19];
                          H_d2 += dN[19] - dN[20];
                          U_sa += dN[17] - dN[21];
                          U_sb += dN[18] - dN[22];
                          U_d1 += dN[16] - dN[23];
                          U_d2 += dN[23] - dN[24];
                          R += dN[13] + dN[15] + dN[21] + dN[22];
                          D += dN[8] + dN[20] + dN[24];
                          
                          // Accumulators
                          a_I += dN[1];
                          a_H += dN[9] + dN[10] + dN[11] + dN[12];
                          a_U += dN[16] + dN[17] + dN[18];
                          a_DI += dN[8];
                          a_DH += dN[20];
                          a_DU += dN[24];
                          a_D  = a_DH + a_DI + a_DU;
                          a_O  += dN[13] + dN[15] + dN[21] + dN[23]; // discharged from hospital
                          a_deltaH = a_H - a_DH - a_DU - a_O;
                          // Current
                          U_curr = U_sa + U_sb + U_d1 + U_d2;
                          H_curr = H_sa + H_sb1 + H_sb2 + H + H_d1 + H_d2 + U_curr;
                          // Total infected
                          tot_I += dN[1];

                          // random walk of beta
                          dWX = rnorm(0, sqrt(dt));
                          
                          X  += (t <= tvary) ? sdfrac * std_X * dWX : std_X * dWX;  // 
                          
                          // susceptibles so as to match total population
                          N = pop - D - H_curr - U_curr;
                          Rt = exp(X) / (i2o); 
                          ")

# Initializer 
init.Csnippet <- Csnippet("X = log(R0_0 * i2o);
                          E   = 0;
                          I1  =  nearbyint(I_0 * pop);
                          I2  = 0;
                          I3 = 0;
                          I_h = 0;
                          I_d = 0;
                          H_sa = 0;
                          H_sb1 = 0;
                          H_sb2 = 0;
                          H = 0;
                          U_sa = 0;
                          U_sb = 0;
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
