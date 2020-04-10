
# Preamble ---------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(foreach)
library(iterators)
source("COVID-pomp/scripts/utils.R")

select <- dplyr::select
# Load hospitalization data ----------------------------------------------------
# Date to replace on the right bound
right_date <- as.Date("2020-04-09")

# Hospital data
hosp_data <- read_csv("data/vd/hospitalization_data_full_20200408.csv") %>% 
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
  mapply(function(x, t) {x[x==0] = .5; mutate(fitErland(x, zero_replace = .75), t = t, obsmean = mean(x))},
         x = obs_times,
         t = names(obs_times),
         SIMPLIFY = F)
)

# Plot to check
p_erlangs  <- foreach(erlang = iter(num_comp, by = "row"),
                      obs = obs_times,
                      .combine = rbind) %do% {
                        tibble(var = erlang$t, 
                               value = obs, 
                               density = NA,
                               type = "obs") %>% 
                          rbind(tibble(var = erlang$t, 
                                       value = seq(0.1, 30, by = .1),
                                       type = "density") %>% 
                                  mutate(density = dgamma(value, shape = erlang$k, rate = erlang$lambda)))
                      } %>% 
  {
    ggplot() +
      geom_histogram(data = filter(., type == "obs"), aes(x = value, y = ..density..), fill = "gray") +
      geom_line(data = filter(., type == "density"), aes(x = value, y = density)) +
      geom_text(data = mutate(num_comp, var = t, label = str_c("k = ", k)), aes(x = 25, y = .35, label = label)) +
      facet_wrap(~var) +
      theme_bw()
  }

plot(p_erlangs)


# Compute probabilities --------------------------------------------------------
hosp_known <- filter(hosp_data, outcome != "hospitalized") %>% 
  mutate(icu = !is.na(icu_in))

hosp_known$deceased[is.na(hosp_known$deceased)] <- F
n_known <- nrow(hosp_known)

pi2hs <- with(hosp_known, sum(outcome == "discharged" & !icu)/n_known)
ph2u <- with(hosp_known, sum(icu)/(sum(icu) + sum(!icu & deceased)))
pu2d <- with(hosp_known, sum(icu & deceased)/sum(icu))




