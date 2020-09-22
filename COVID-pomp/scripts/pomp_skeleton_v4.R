
# Set variables -----------------------------------------------------------


# Number of compartements for each variable to represent Erlang distributions
n_compartments <- list(
  E = nc_E,
  I = nc_I
)

# define stat variable names for each district
state_names <- mapply(
  function(n, comp) {str_c(comp, 1:n)},
  n = n_compartments, 
  comp = names(n_compartments)) %>% 
  unlist() %>% 
  c("S","I_d","I_h","H_s", "H", "U_s", "H_d","U_d","R" ,"H_curr", "U_curr", "D", "X", "N",
    "a_I", "a_H", "a_U", "a_D", "a_DH", "a_DI", "a_DU", "a_O", "a_OU", "a_deltaH",
    "a_deltaU","a_deltaID", "Rt", "tot_I") # prefix a_ represent accumulator variables for incidences

# define parameter names for pomp
## process model parameters names to estimate

## initial value parameters to estimate
param_iv_est_names <- c()

## fixed process model parameters 
rate_names <- c("e2i", "i2h", "id2o", "i2o", "hs2r", "hd2d", "h2u", "us2r", "ud2d")
prob_names <- c("psevere", "pi2d", "pi2h", "pi2hs", "ph2u", "pu2d")
regular_names <- c(rate_names, prob_names) # regular parameters (as opposed to initial value parameters)

# define parameter names for pomp
## process model parameters names to estimate
if (use_case_incid) {
  param_proc_est_names <- c("std_X", "k", "epsilon")
  param_proc_fixed_names <- c("pop", rate_names, prob_names)
} else {
  param_proc_est_names <- c("std_X")
  param_proc_fixed_names <- c("pop", rate_names, prob_names, "k", "epsilon")
}

if (!is.null(config$parameters_to_fit)) {
  regular_to_fit <- names(config$parameters_to_fit)[names(config$parameters_to_fit) %in% regular_names]
  iv_to_fit <- names(config$parameters_to_fit)[!(names(config$parameters_to_fit) %in% regular_to_fit)]
  
  allin <- unlist(lapply(regular_to_fit, function(x) !(x %in% regular_names)))
  if (sum(allin) == 0) {
    param_proc_est_names <- c(param_proc_est_names, regular_to_fit)
    param_proc_fixed_names <- param_proc_fixed_names[!(param_proc_fixed_names %in% names(config$parameters_to_fit))]
  }
  
  if (length(iv_to_fit) > 0) {
    param_iv_est_names <- iv_to_fit
    param_iv_fixed_names <- c()
  } else {
    param_iv_fixed_names <- c("I_0", "R0_0")
  }
}

# all parameter names to estimate
param_est_names <- c(param_proc_est_names, param_iv_est_names)
# all fixed parameters
param_fixed_names <- c(param_proc_fixed_names, param_iv_fixed_names)

# all param names
param_names <- c(param_est_names, param_fixed_names)

# names of parameters that are rates
param_rates_in_days_names <- rate_names

# Measurment model  -------------------------------------------------------

# measurement model
# dmeasure.Csnippet <- Csnippet("
# double ll_na;
# ll_na = (give_log) ? 0 : 1;
# double ll_c, ll_h, ll_hc, ll_d, ll_deltah, ll_deltau, ll_dh, ll_du, ll_di;
#                                if (!ISNA(delta_hosp)) { 
#                                // if info on releases not available, give hosp curr
#      ll_deltah = dskellam(delta_hosp, a_H, a_DH + a_DU + a_O, give_log);
#      } else {
#        ll_deltah = ll_na;
#      } lik = ll_deltah;
#      ")
dmeasure.Csnippet <- Csnippet(glue("
  double ll_na;
  double ll_c, ll_h, ll_hc, ll_d, ll_deltah, ll_deltau, ll_dh, ll_du, ll_di;
  ll_na = (give_log) ? 0 : 1;
  if (!ISNA(case_incid)) {
    ll_c = dnbinom_mu(case_incid, k, epsilon * a_I, give_log);
  } else {
    ll_c = ll_na;
  }
 // if (!ISNA(ll_deltau)) {
      // if info on releases not available, give hosp curr
  //   ll_deltau = dskellam(delta_icu, a_U, a_DU + a_OU, give_log);
 //  } else {
  //    ll_deltau = ll_na;
 //   }

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

  if (!ISNA(death_incid)) {
   ll_d = dpois(death_incid, a_D, give_log);
  } else {
    ll_d = ll_na;
  }

  if (give_log) {
    lik =  {{parsed_lik$lsum}};
  } else {
    lik =  {{parsed_lik$lprod}};
  }
                              ", .open = "{{", .close = "}}"))

## NegBinomial simulator
rmeasure.Csnippet <- Csnippet("
                      int deltah = rskellam(a_H, a_DH + a_DU);
                      death_incid = rpois(a_D);
                      hosp_curr = rpois(H_curr);
                      hosp_incid = rpois(a_H);
                      case_incid = rnbinom_mu(k, epsilon * a_I);
                      delta_hosp = deltah;
                      delta_ID = deltah + rpois(a_O);
                      ")

# Process model -----------------------------------------------------------------

Icomp <- state_names[state_names %>% grep('^I[0-9]*$', .)]
Ecomp <- state_names[state_names %>% grep('^E[0-9]*$', .)]
Isum <- paste(Icomp, collapse = '+')
Iinit <- paste(append(Icomp, ' '), collapse = '=0;')
Einit <- paste(append(Ecomp, ' '), collapse = '=0;')
Edraw <- ''
Idraw <- ''
Eeqn <- glue('E1 -= dN[19];')
Ieqn <- glue('I1 -= dN[{18+nc_E+1}];')
for (i in seq(nc_E-1)){
  Edraw <- str_c(Edraw,glue('reulermultinom(1, E{i}, &rate[1], dt, &dN[{18+i}]);'))
  if (i > 1){
  Eeqn <- str_c(Eeqn,glue('E{i} += dN[{18+i-1}] - dN[{18+i}];'))
  }
  
}
for (i in seq(nc_I-1)){
  Idraw <- str_c(Idraw,glue('reulermultinom(1, I{i}, &rate[2], dt, &dN[{18+nc_E+i}]);'))
  if (i > 1){
    Ieqn <- str_c(Ieqn,glue('I{i} += dN[{18+nc_E+i-1}] - dN[{18+nc_E+i}];'))
  }
}
if (nc_E == 1)
{
  Edraw <- ''
  Eeqn <- ''
} else {
  Eeqn <- str_c(Eeqn,glue('E{nc_E} += dN[{18+nc_E-1}];'))
}

if (nc_I == 1)
{
  Idraw <- ''
  Ieqn <- ''
} else {
  Ieqn <- str_c(Ieqn,glue('I{nc_I}   += dN[{18+nc_E+nc_I-1}];'))
}

# create C code for each district
proc.Csnippet <- Csnippet(glue("
                          double foi, foi_stoc; // force of infection and its stochastic version
                          double dw, dWX;       // extra-demographic stochasticity on foi
                          double rate[18];      // vector of all rates in model
                          double dN[18+{nc_E}+{nc_I}];        // vector of transitions between classes during integration timestep

                          // force of infection
                          
                          foi =  exp(X) * ({Isum})/N;
                          
                          // define transition rates for each type of event
                          // S compartment
                          rate[0] = foi;   // infections
                          // E compartment
                          rate[1] = e2i * {nc_E};    // transition to I1
                          // I1 compartment
                          rate[2] = i2o * {nc_I};  // transition to I2

                          // I3 compartment
                          rate[4] = psevere * (1-pi2h) * i2o * 3;    // death_incid without hospitalization
                          rate[5] = psevere * pi2h * i2o * 3;    // hospitalizations
                          rate[6] = (1-psevere) * i2o * 3;      // recovery from infection
                          
                          // I_d compartment
                          rate[7] = (1-pi2d) * id2o; // severe infected that recover
                          rate[8] = pi2d * id2o;    // severe infected that die
                          
                          // I_h compartment
                          rate[9] =  i2h * pi2hs;   // hospitalization that WILL NOT go in ICU nor die
                          rate[10] = i2h * (1-pi2hs) * (1-ph2u);   // hospitalization that WILL NOT go in ICU and die
                          rate[11] = i2h * (1-pi2hs) * ph2u;   // hospitalization that WILL go in ICU or die
                          
                          // H_s compartment
                          rate[12] = hs2r;   // transition to recovery
                          
                          // H compartment
                          rate[13] = h2u * pu2d; // ICU admission that die
                          rate[14] = h2u * (1-pu2d); // ICU admission that survive
                          
                          // H_d compartment
                          rate[15] = hd2d; // Hosp to death
                          
                          // U_s compartment
                          rate[16] = us2r;    // recovery
                          
                          // U_d compartment
                          rate[17] = ud2d;   // ICU to death
                          
                          // simulate all transitions
                          reulermultinom(1, S,  &rate[0], dt, &dN[0]);
                          {Edraw}
                          reulermultinom(1, E{nc_E}, &rate[1], dt, &dN[1]);
                          
                          {Idraw}
                          reulermultinom(3, I{nc_I}, &rate[4], dt, &dN[4]);
                          
                          reulermultinom(2, I_d, &rate[7], dt, &dN[7]);
                          reulermultinom(3, I_h, &rate[9], dt, &dN[9]);
                          reulermultinom(1, H_s, &rate[12], dt, &dN[12]);
                          reulermultinom(2, H,    &rate[13], dt, &dN[13]);
                          reulermultinom(1, H_d, &rate[15], dt, &dN[15]);
                          reulermultinom(1, U_s,  &rate[16], dt, &dN[16]);
                          reulermultinom(1, U_d, &rate[17], dt, &dN[17]);


                          // update state variables
                          S    += -dN[0];
                          E1   += dN[0];
                          {Eeqn}
                          E{nc_E} -= dN[1];
                        
                          I1   += dN[1];
                          {Ieqn}
                          I{nc_I}   += - dN[4] - dN[5] - dN[6];

                          I_d  += dN[4] - dN[7] - dN[8];
                          I_h  += dN[5] - dN[9] - dN[10] - dN[11];
                          H_s += dN[9] - dN[12];
                          H += dN[11] - dN[13] - dN[14];
                          H_d += dN[10] - dN[15];
                          U_s += dN[14] - dN[16];
                          U_d += dN[13] - dN[17];
                          R += dN[6] + dN[7] + dN[12] + dN[16];
                          D += dN[8] + dN[15] + dN[17];
                          
                          // Accumulators
                          a_I += dN[1];
                          a_H += dN[9] + dN[10] + dN[11];
                          a_U += dN[13] + dN[14];
                          a_DI += dN[8];
                          a_DH += dN[15];
                          a_DU += dN[17];
                          a_D  = a_DH + a_DI + a_DU;
                          a_OU += dN[16]; // discharged from ICU
                          a_O  += dN[12] + dN[16]; // discharged from hospital
                          a_deltaH = a_H - a_DH - a_DU - a_O;
                          a_deltaID = a_H - a_DH - a_DU;
                          a_deltaU = a_U - a_DU - a_OU;
                          // Current
                          U_curr = U_s + U_d;
                          H_curr = H_s + H + H_d + U_curr;
                          // Total infected
                          tot_I += dN[1];

                          // random walk of beta
                          dWX = rnorm(0, sqrt(dt));
                          
                          X  += (t <= tvary) ? sdfrac * std_X * dWX : std_X * dWX;  // 
                          
                          // susceptibles so as to match total population
                          N = pop - D - H_curr - U_curr;
                          Rt = exp(X) / (i2o); 
                          "))

# Initializer 
init.Csnippet <- Csnippet(glue("X = log(R0_0 * i2o);
                               {Einit}
                               {Iinit}
                          I1  =  nearbyint(I_0 * pop);
                          I_h = 0;
                          I_d = 0;
                          H_s = 0;
                          H = 0;
                          U_s = 0;
                          H_d = 0;
                          U_d = 0;
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
                          a_OU = 0;
                          a_deltaH = 0;
                          a_deltaID = 0;
                          a_deltaU = 0;
                          R   = 0;
                          S   = nearbyint(pop - I1);
                          N   = pop;
                          tot_I = 0;"))
