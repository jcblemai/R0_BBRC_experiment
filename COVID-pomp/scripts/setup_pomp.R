
# Preamble ---------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(foreach)
library(iterators)
source("COVID-pomp/scripts/utils.R")

select <- dplyr::select
# Load hospitalization data ----------------------------------------------------
# Date to replace on the right bound
right_date <- as.Date("2020-04-14")

# Hospital data
hosp_data <- read_csv("data/vd/hospitalization_data_full_20200414.csv") %>% 
  mutate_at(c("date_in", "date_out", "icu_in", "icu_out"), 
            function(x) as.Date(x, format = "%m/%d/%Y")) %>%  
  filter(date_in > "2020-02-28") %>% 
  mutate(deceased = as.logical(deceased),
         outcome = case_when(deceased ~ "deceased",
                             !is.na(date_out) ~ "discharged",
                             T ~ "hospitalized"), 
         time_to_death = difftime(date_out, date_in, units = "days"),
         time_icu_to_death = difftime(date_out, icu_in, units = "days"),
         date_out = case_when(is.na(date_out) ~ right_date,
                              T ~ date_out),
         time_hosp = difftime(date_out, date_in, units = "days"),
         icu_state = case_when(!is.na(icu_in) & !is.na(icu_out) ~ "discharged from ICU",
                               !is.na(icu_in) & is.na(icu_out) ~ "in ICU",
                               T ~ ""),
         icu_out = case_when(!is.na(icu_in) & is.na(icu_out) ~ right_date,
                             !is.na(icu_out) ~ icu_out),
         time_icu = difftime(icu_out, icu_in, units = "days"),
         time_to_icu= difftime(icu_in, date_in, units = "days")) %>% 
  mutate_at(vars(contains("time_")), as.numeric)


hosp_cumul <-  bind_rows(
  lapply(seq.Date(min(hosp_data$date_in), max(hosp_data$date_out), by = "1 days"),
         function(x) {
           with(hosp_data,
                tibble(
                  date = x,
                  hosp_incid = sum(date_in == x),
                  icu_incid = sum(icu_in == x, na.rm = T),
                  deaths_icu_incid = sum(date_out == x & !is.na(icu_in) & outcome == "deceased"),
                  deaths_noicu_incid = sum(date_out == x & is.na(icu_in) & outcome == "deceased"),
                  hosp_curr = sum(date_in <= x & (date_out > x | (date_out ==x & outcome == "hospitalized"))),
                  icu_curr = sum(icu_in <= x & (icu_out > x |(icu_out == x & outcome == "hospitalized")), na.rm = T),
                  r_incid =  sum(date_out == x & outcome == "discharged"),
                  r_cumul =  sum(date_out <= x  & outcome == "discharged"),
                  o_cumul =  sum(date_out <= x  & outcome != "hospitalized"),
                  d_cumul =  sum(date_out <= x  & outcome == "deceased"),
                  hosp_cumul =  sum(date_in <= x)
                )
           )
         })) %>% 
  mutate(predd_cumul = hosp_cumul*0.1)

write_csv(hosp_cumul, path = "data/VD_hosp_data.csv")

hosp_cumul %>% 
  gather(var, value, -date) %>%  filter(str_detect(var, "cumul")) %>% 
  ggplot(aes(x = date, y = value, color = var)) +
  geom_line() #+
facet_wrap(~ var, scales = "free_y")

# Times ------------------------------------------------------------------------

obs_times <- list(
  h2icu = filter(hosp_data, !is.na(time_to_icu), icu_state != "")[["time_to_icu"]],
  icu2d = filter(hosp_data, outcome == "deceased", !is.na(icu_in))[["time_icu_to_death"]],
  icu2r = filter(hosp_data, outcome != "deceased", !is.na(icu_in))[["time_icu"]],
  h2d_noicu = filter(hosp_data, outcome == "deceased", is.na(icu_in))[["time_to_death"]],
  h2r = filter(hosp_data, outcome != "deceased", is.na(icu_in))[["time_hosp"]]
)

# Number of compartments
num_comp <- bind_rows(
  mapply(function(x, t) {
    mutate(fitErland(x, zero_replace = .75, ks = 1:15),
           t = t,
           obsmean = mean(x),
           obsvar = var(x))
  },
  x = obs_times,
  t = names(obs_times),
  SIMPLIFY = F)
)


# Compute probabilities --------------------------------------------------------
hosp_known <- filter(hosp_data, outcome != "hospitalized", outcome_ori != "transfert") %>% 
  mutate(icu = !is.na(icu_in))

hosp_known$deceased[is.na(hosp_known$deceased)] <- F
n_known <- nrow(hosp_known)

# Probability of going to ICU
pu <- sum(!is.na(hosp_data$icu_in))/nrow(hosp_data)
# pi2hs <- with(hosp_known, sum(outcome == "discharged" & !icu)/n_known)
# ph2u <- with(hosp_known, sum(icu)/(sum(icu) + sum(!icu & deceased)))
# Compute ICU CFR
ICUcfr <- computeCCFR(obs_times$icu2d,
                      cumsum(hosp_cumul$icu_incid)[-1],
                      cumsum(hosp_cumul$deaths_icu_incid)[-1], hosp_cumul$date[-1], "lnorm")
pu2d <- ICUcfr$mle[nrow(ICUcfr)]  # with(hosp_known, sum(icu & deceased)/sum(icu))
hcfr <- 0.11
ph2u <- 1/(hcfr/pu - pu2d + 1) 
pi2hs <- 1 - pu/ph2u  
hcfr2 <- (1- pi2hs) * (ph2u * pu2d + 1 - ph2u)

# IFR and prob of hospitaliztion

computePHosp <- function(x, pi2d, hfcr) {
x*pi2d/(1-x)/(hcfr + x*pi2d/(1-x))
}

toOptPI2D <- function(pi2d, ifr, psevere, x, hcfr) {
  pi2h <- computePHosp(x, pi2d, hcfr)
  (ifr - psevere*(pi2h * hcfr + (1 - pi2h) * pi2d))^2
}

x <- 0.5
# proportion of deaths
psymptomatic <- 0.5
ifr <- 0.0075
psevere <- psymptomatic * 0.15
pi2d <- optimize(toOptPI2D, ifr, psevere, x, hcfr, interval = c(0,1))$minimum
# pi2d <- 0.1
pi2h <- computePHosp(x, pi2d, hcfr)
# Check proportion dying in hosp
pi2h*hcfr/(pi2h*hcfr + (1-pi2h) * pi2d)
# Check the ifr
psevere*(pi2h * hcfr + (1 - pi2h) * pi2d)

# Rstan for mixture ------------------------------------------------------------
library(rstan)
options(mc.cores = parallel::detectCores())  

mixture_erland_stan <- stan_model("COVID-pomp/scripts/erlang_mixture.stan")
mixture_gamma_stan <- stan_model("COVID-pomp/scripts/gamma_mixture.stan")
n_chains <- 4

data <- with(filter(hosp_data, 
                    outcome != "deceased", 
                    !is.na(icu_in), 
                    is.na(outcome_ori) | outcome_ori != "transfert",
                    time_icu > 0),
             list(N = length(time_icu),
                  C = 8,
                  y = time_icu,
                  C_prior = rep(1/8, 8)
                  # C_prior = c(.25, .25, .15, .1, .1, .05, .05, .05)
                  ))

datah <- with(filter(hosp_data, 
                    outcome != "deceased", 
                    is.na(icu_in), 
                    # is.na(outcome_ori) | outcome_ori != "transfert",
                    time_hosp > 0),
             list(N = length(time_hosp),
                  C = 8,
                  y = time_hosp,
                  C_prior = rep(1/8, 8)
             ))


# Fit the model
control <- list(adapt_delta = .9, metric = "dense_e", max_treedepth = 12)
nwarmup <- 4000
niter <- nwarmup + 1000
singleerlang_stanfit <- sampling(mixture_erland_stan,
                                 data = append(datah, list(K = 1, mu_scale = c(4, 0))),
                                 chains = n_chains,
                                 warmup = nwarmup, iter = niter,
                                 control = control,
                                 save_warmup = FALSE)

mixerlang_stanfit <- sampling(mixture_erland_stan,
                              data = append(data, list(K = 2, mu_scale = c(2, 10))),
                              chains = n_chains,
                              warmup = nwarmup, iter = niter,
                              control = control,
                              save_warmup = FALSE)

loo_2k <- loo::loo(mixerlang_stanfit)
loo_1k <- loo::loo(singleerlang_stanfit)

loo::loo_compare(loo_1k, loo_2k)

a <- summary(mixerlang_stanfit, pars = "lp")$summary 
a %>% as_tibble()  %>% 
  mutate(param = rownames(a),
         obid  = as.numeric(str_extract(param, "(?<=\\[)[0-9]+(?=,)")),
         grp = as.numeric(str_extract(param, "(?<=,)[0-9]+(?=,)")),
         num_comp = as.numeric(str_extract(param, "(?<=,)[0-9]+(?=\\])"))) %>% 
  group_by(num_comp, grp) %>% 
  summarise(ll = sum(mean)) %>% 
  ggplot(aes(x = num_comp, y = ll, color = factor(grp))) +
  geom_line() +
  facet_wrap(~grp, scales ="free_y")

comp_summary <- summary(mixerlang_stanfit, pars = "ll_c")$summary 
comp_post <-comp_summary %>% as_tibble() %>% 
  mutate(param = rownames(comp_summary),
         grp  = as.numeric(str_extract(param, "(?<=\\[)[0-9]+(?=,)")),
         num_comp = as.numeric(str_extract(param, "(?<=,)[0-9]+(?=\\])")))

comp_post %>% select(mean, num_comp, grp) %>% 
  gather(var, value, -grp, -num_comp) %>% 
  ggplot(aes(x = num_comp, y = value, fill = factor(grp))) +
  geom_bar(position = "dodge", stat = "identity")

k_summary <- summary(mixerlang_stanfit, pars = "ll_k")$summary
group_post <- k_summary %>%
  as_tibble() %>%
  mutate(param = rownames(k_summary),
         obsid = as.numeric(str_extract(param, "(?<=\\[)[0-9]+(?=,)")),
         grp = str_extract(param, "(?<=,)[0-9]+(?=\\])")) %>% 
  group_by(obsid) %>%
  summarise(est_grp = grp[which.max(mean)],
            lambda = mean[grp=="1"]) 


icu_data <- filter(hosp_data, 
                   outcome != "deceased", 
                   !is.na(icu_in), 
                   # is.na(outcome_ori) | outcome_ori != "transfert",
                   time_icu > 0) %>% 
  cbind(group_post)

ggplot(icu_data, aes(x=time_icu, y = 1, fill = lambda)) + 
  geom_bar(stat = "identity")  +
  scale_fill_gradient2(midpoint = .5)


ggplot(icu_data, aes(x = lambda, fill = sex)) + 
  geom_histogram() +
  scale_color_viridis_c()



# Plots ------------------------------------------------------------------------

# Add two-group model for ICU
num_comp <- mutate(num_comp, grp = "1")
num_comp2 <- rbind(num_comp, 
                   tibble(k = c(1, 1), ll = NA, 
                          mean = summary(mixerlang_stanfit, pars = "mu")$summary[,1],
                          t = "icu2r", obsmean =NA, obsvar = NA,
                          grp = c("2", "3")) %>% 
                     mutate(lambda = k/mean,
                            variance = mean^2/k) %>% 
                     select(one_of(colnames(num_comp))))

# Plot to check
p_erlangs  <- foreach(erlang = iter(num_comp2, by = "row"),
                      obs = append(obs_times, list(obs_times$icu2r, obs_times$icu2r)),
                      .combine = rbind
) %do% {
  tibble(var = erlang$t, 
         value = obs, 
         density = NA,
         type = "obs") %>% 
    rbind(tibble(var = erlang$t, 
                 value = seq(0.1, 30, by = .1),
                 type = "density") %>% 
            mutate(density = dgamma(value, shape = erlang$k, rate = erlang$lambda))) %>% 
    mutate(grp = erlang$grp)
} %>% 
  {
    ggplot() +
      geom_histogram(data = filter(., type == "obs"), aes(x = value, y = ..density..), fill = "gray") +
      geom_line(data = filter(., type == "density"), aes(x = value, y = density, color = grp)) +
      geom_text(data = mutate(num_comp, var = t, label = str_c("k = ", k)), aes(x = 25, y = .35, label = label)) +
      facet_wrap(~var) +
      theme_bw()
  }

plot(p_erlangs)



# Scraps -----------------------------------------------------------------------



# %>% 
#   mutate(lambda = case_when(est_grp == 1 ~ lambda, T ~ 1-lambda)) %>% 
#   ggplot(aes(x = lambda, y = grp)) +
#   geom_point(alpha = .1)
# 
n_test <- 100
lambdas <- c(.8, .2)
mix <- runif(n_test) < lambdas[1]
test_data <- c(rgamma(n_test, shape = 1, scale = 2)[mix],
               rgamma(n_test, shape = 10, scale = 4)[!mix])
data <- list(N = n_test,
             K = 2,
             C = 10,
             y = test_data)
#              
#              lambda_est <- summary(mixerlang_stanfit, pars = "ll_k")$summary
# res <- lambda_est %>%
#   as_tibble() %>% 
#   mutate(param = rownames(lambda_est),
#          obsid = str_extract(param, "(?<=\\[)[0-9]+(?=,)"),
#          grp = str_extract(param, "(?<=,)[0-9]+(?=\\])")) %>% 
#   group_by(obsid) %>% 
#   summarise(lambda = max(mean),
#             est_grp = grp[which.max(mean)]) %>% 
#   mutate(est_grp = factor(est_grp),
#          grp = mix,
#          grp = case_when(mix ~ 1, T ~ 2),
#          grp = factor(grp)) %>% 
#   select(grp, est_grp) 
# 
# caret::confusionMatrix(res$est_grp, res$grp)
# mixgamma_stanfit <- sampling(mixture_gamma_stan,
#                              data = append(data, list(K = 2, mu_scale = c(2, 10))),
#                              chains = n_chains,
#                              warmup = nwarmup, iter = niter,
#                              control = control,
#                              save_warmup = FALSE)
# 
# singlegamma_stanfit <- sampling(mixture_gamma_stan,
#                                 data =  append(data, list(K = 1, mu_scale = c(4, 0))),
#                                 chains = n_chains,
#                                 warmup = nwarmup, iter = niter,
#                                 control = control,
#                                 save_warmup = FALSE)