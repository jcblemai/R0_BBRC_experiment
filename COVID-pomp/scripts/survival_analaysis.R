# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(foreach)
library(iterators)
library(rstan)
library(glue)
library(survival)
options(mc.cores = parallel::detectCores())  

source("scripts/utils.R")

posteriorSurvivals <- function(fit, pars, dist = "lnorm", times = c(.5, 1:90)) {
  post_pars <- rstan::extract(fit, c(pars, "lambda"))
  densities <- lapply(
    1:2, 
    function(x) {
      lapply(times, function(t) {
        p1 <- post_pars[[pars[1]]][, x]
        p2 <- post_pars[[pars[2]]][, x]
        
        if(dist == "lnorm") {
          cdfs <-  plnorm(t, meanlog = p1, sdlog = p2)
        } else if (dist == "gamma") {
          cdfs <-  pgamma(t, shape = p1, rate = p1/p2)
        }
        
        post_pars$lambda[, x]*(cdfs)
      }
      ) %>% 
        do.call(cbind, .) %>% 
        apply(2, function(i) data.frame(median = median(i),
                                        mean = mean(i),
                                        q025 = quantile(i, 0.025),
                                        q975 = quantile(i, 0.975))) %>% 
        do.call(rbind, .) %>% 
        mutate(group = x,
               time = times)
    }) %>% 
    do.call(rbind, .)
  return(densities)
}
meanLognorm <- function(mu, sigma) { exp(mu + sigma^2/2)}
varLognorm <- function(mu, sigma) {exp(sigma^2 - 1)*exp(sigma^2 + 2 * mu)}

# Load data --------------------------------------------------------------------

hosp_data <- read_csv("data/vd/CHUV_hosp_data_cleaned_2020-04-25.csv") %>% 
  mutate(uid2 = str_c(uid, date_in, sep = "_"))

hosp_data$deceased[which(!hosp_data$deceased)] <- NA
right_date <- as.Date("2020-04-25")
event_data <- hosp_data %>%  
  filter(date_in > "2020-02-24") %>% 
  mutate(hosp_event = case_when(is.na(date_out) | !is.na(transfert) ~ "censor",
                                !is.na(date_out) & is.na(deceased) ~ "discharged",
                                T ~ "deceased"),
         icu_event = case_when(!is.na(icu_in) & 
                                 (is.na(icu_out) | 
                                    (date_out == icu_out & !is.na(transfert))) ~ "censor",
                               !is.na(icu_in) & !is.na(icu_out) & !( 
                                 icu_out == date_out& !is.na(deceased)) ~ "discharged",
                               !is.na(icu_in) & !is.na(icu_out) & !is.na(deceased) ~ "deceased",
                               T ~ as.character(NA)),
         hosp_event = factor(hosp_event, levels = c("censor", "deceased", "discharged")),
         icu_event = factor(icu_event, levels = c("censor", "deceased", "discharged")))

event_data <- event_data %>% 
  mutate(icu_out = case_when(is.na(icu_out) ~ right_date, T ~ icu_out),
        date_out = case_when(is.na(date_out) ~ right_date, T  ~ date_out),
        time_icu = as.numeric(difftime(icu_out, icu_in, units = "days")),
        time_hosp = as.numeric(difftime(date_out, date_in, units = "days")))

# Survival analysis ------------------------------------------------------------

surv_out_path <- str_c("COVID-pomp/results/survival_analaysis")

# Compile stan model
survival_lnorm_stan <- stan_model("COVID-pomp/scripts/survival_competing_risk_lnorm.stan")
survival_gamma_stan <- stan_model("COVID-pomp/scripts/survival_competing_risk_gamma.stan")
# survival_lnorm_mix_stan <- stan_model("COVID-pomp/scripts/survival_competing_risk_mixture_exp_lnorm.stan")

# Stan HMC 
n_chains <- 4
nwarmup <- 2000
niter <- nwarmup + 1000
control <- list(adapt_delta = .9, metric = "dense_e", max_treedepth = 12)

namecols <- c("none", "deceased", "discharged")

states <- foreach(
  in_icu = list(c(F, T), c(F), c(T)),
  out_suffix = c("all", "noICU", "ICU"),
  .combine = rbind) %do% {
    
    dat <- event_data %>% 
      filter(as.character(!is.na(icu_in)) %in% as.character(in_icu)) 
    
    if (in_icu[1]) {
      dat$time <- dat$time_icu
      dat$event <- dat$icu_event
    } else {
      dat$time <- dat$time_hosp
      dat$event <- dat$hosp_event
    }
    dat$time[dat$time == 0] <- .5
    dat <- select(dat, time, event)
    
    km_fit <- survfit(Surv(dat$time, dat$event, type = "mstate") ~ 1)
    km_surv <- km_fit$pstate %>% 
      as_tibble() %>% 
      set_colnames(namecols) %>% 
      mutate(time = km_fit$time) %>% 
      gather(var, mean, -time) %>% 
      inner_join(
        km_fit$lower %>% 
          as_tibble() %>% 
          set_colnames(namecols) %>% 
          mutate(time = km_fit$time) %>% 
          gather(var, q025, -time)
      ) %>% 
      inner_join(
        km_fit$upper %>% 
          as_tibble() %>% 
          set_colnames(namecols) %>% 
          mutate(time = km_fit$time) %>% 
          gather(var, q975, -time)
      ) %>% 
      mutate(median = NA,
             group = case_when(var == "deceased" ~ 1, T ~ 2),
             dist = "Aalen-Johansen",
             subset = out_suffix) %>% 
      filter(var != "none") %>% 
      select(-var)
    
    # Bayesian mixture model
    # Prepare data
    times <- dat$time
    N <- nrow(dat)
    M <- length(unique(dat$event[dat$event != "censor"]))
    event_indices <- lapply(levels(dat$event), function(i) which(dat$event == i))
    n_events <- map_dbl(event_indices, ~length(.))
    indices_mat <- lapply(event_indices, function(x) c(x, rep(0, max(n_events) - length(x)))) %>% 
      do.call(rbind, .)
    
    # Combine data
    surv_data <- list(N = N,
                      M = M,
                      times = times,
                      N_events = n_events,
                      indices = indices_mat)
    
    cat("Events mapping:", paste(0:M, levels(dat$event), sep = ":"))
    
    outfile <- glue("{surv_out_path}_{out_suffix}_stanfits.rda")
    
    if(!file.exists(outfile)) {
      # Fit model with one group per state
      survival_stanfit_lnorm <- sampling(survival_lnorm_stan,
                                         data = surv_data,
                                         chains = n_chains,
                                         warmup = nwarmup, iter = niter,
                                         control = control,
                                         save_warmup = FALSE)
      
      survival_stanfit_gamma <- sampling(survival_gamma_stan,
                                         data = surv_data,
                                         chains = n_chains,
                                         warmup = 3000, iter = 4000,
                                         control = control,
                                         save_warmup = FALSE)
      
      
      # # Fit model with one death and two recovered states
      # survival_stanfit_lnorm_mixR <- sampling(survival_lnorm_mix_stan,
      #                                         data = append(surv_data, list(K = c(1, 2))),
      #                                         chains = 4,
      #                                         warmup = nwarmup, iter = niter,
      #                                         control = control,
      #                                         save_warmup = FALSE)
      
      # Save fits
      save(survival_stanfit_lnorm, 
           survival_stanfit_gamma, 
           # survival_stanfit_lnorm_mixR,
           file = outfile)
    } else {
      load(outfile)
    }
    # Write distribution results
    pars <- rbind(
      summary(survival_stanfit_gamma, 
              pars = c("mu", "alpha", "lambda"))$summary %>% 
        as.data.frame() %>% 
        mutate(param = rownames(.),
               model = "gamma_ungrouped"),
      summary(survival_stanfit_lnorm, 
              pars = c("mu", "sigma", "lambda"))$summary %>% 
        as.data.frame() %>% 
        mutate(param = rownames(.),
               model = "lnorm_ungrouped")#,
      # summary(survival_stanfit_lnorm_mixR, 
      #         pars = c("mu", "sigma", "lambda", "theta[2]", "mu_exp[2]"))$summary  %>% 
      #   as.data.frame() %>% 
      #   mutate(param = rownames(.),
      #          model = "lnorm_grouped")
      ) %>% 
      select(-n_eff, -Rhat, -se_mean) %>% 
      mutate(subset = out_suffix)
    
    write_csv(pars, path = glue("{surv_out_path}_{out_suffix}_params.csv"))
    
    # Return survivals only for ungrouped model
    lnorm_surv <- posteriorSurvivals(survival_stanfit_lnorm, c("mu", "sigma")) %>% 
      mutate(dist = "lnorm",
             subset = out_suffix)
    
    gamma_surv <- posteriorSurvivals(survival_stanfit_gamma, dist = "gamma", c("alpha", "mu")) %>% 
      mutate(dist = "gamma",
             subset = out_suffix)
    
    rbind(km_surv, lnorm_surv, gamma_surv)
    # lnorm_surv_mix <- posteriorSurvivals(survival_stanfit_lnorm_mix, c("mu", "sigma")) %>% mutate(dist = "lnorm")
  }

write_csv(states, "results/survival_analaysis_curves.csv")

compare_models <- foreach(f = c("COVID-pomp/results/survival_analaysis_noICU_stanfits.rda", 
                                "COVID-pomp/results/survival_analaysis_ICU_stanfits.rda"),
                          .combine = rbind
                          
) %do% {
  load(f)
  loo_lognorm <- loo::loo(survival_stanfit_lnorm)
  loo_gamma <- loo::loo(survival_stanfit_gamma)
  # loo_lognorm_mix <- loo::loo(survival_stanfit_lnorm_mixR)
  compare <- loo::loo_compare(loo_lognorm, loo_gamma)
  as.data.frame(compare) %>% 
    mutate(subset = str_extract(f, "(?<=sis_)(.*)(?=_stan)"),
           model = rownames(.))
}


# rbind(lnorm_surv, gamma_surv) %>% 
filter(states, dist != "Aalen-Johansen") %>% 
  ggplot(aes(x = time, y = mean, ymin = q025, ymax = q975, fill = dist)) +
  geom_ribbon(alpha = .2) +
  geom_line(aes(color = dist)) +
  pammtools::geom_stepribbon(data = filter(states, dist == "Aalen-Johansen"), alpha = .2) +
  geom_step(data = filter(states, dist == "Aalen-Johansen"), aes(color = dist)) +
  facet_grid(subset~group, labeller = labeller(group = c("1" = "death", "2" = "recovery"))) +
  ylab("1 - Survival function ") +
  theme_bw()


pars_hosp <- read_csv("COVID-pomp/results/survival_analaysis_noICU_params.csv")
pars_icu <- read_csv("COVID-pomp/results/survival_analaysis_ICU_params.csv")

mu_hosp_d <- pars_hosp$mean[pars_hosp$param == "mu[1]" & pars_hosp$model == "lnorm_ungrouped"]
sigma_hosp_d <- pars_hosp$mean[pars_hosp$param == "sigma[1]" & pars_hosp$model == "lnorm_ungrouped"]
mu_icu_d <- pars_icu$mean[pars_icu$param == "mu[1]" & pars_icu$model == "lnorm_ungrouped"]
sigma_icu_d <- pars_icu$mean[pars_icu$param == "sigma[1]" & pars_icu$model == "lnorm_ungrouped"]

mu_hosp_r <- pars_hosp$mean[pars_hosp$param == "mu[2]" & pars_hosp$model == "lnorm_ungrouped"]
sigma_hosp_r <- pars_hosp$mean[pars_hosp$param == "sigma[2]" & pars_hosp$model == "lnorm_ungrouped"]
mu_icu_r <- pars_icu$mean[pars_icu$param == "mu[2]" & pars_icu$model == "lnorm_ungrouped"]
sigma_icu_r <- pars_icu$mean[pars_icu$param == "sigma[2]" & pars_icu$model == "lnorm_ungrouped"]

meanLognorm(mu_icu_r, sigma_icu_r)
varLognorm(mu_icu_r, sigma_icu_r)

meanLognorm(mu_icu_d, sigma_icu_d)
varLognorm(mu_icu_d, sigma_icu_d)

x <- seq(0.5, 100, by = .5)
# Generate observations for comp fitting
d_h2d <- dlnorm(x, mu_hosp_d, sigma_hosp_d)
d_icu2d <- dlnorm(x, mu_icu_d, sigma_icu_d)
d_h2r <- dlnorm(x, mu_hosp_r, sigma_hosp_r)
d_icu2r <- dlnorm(x, mu_icu_r, sigma_icu_r)
