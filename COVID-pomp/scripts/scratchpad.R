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



state_names <- mapply(
  function(n, comp) {if(n == 1) {comp} else {str_c(comp, 1:n)}},
  n = n_compartments, 
  comp = names(n_compartments)) %>% 
  unlist() %>% 
  c("H_curr", "U_curr", "D", "X", "N",
    "a_I", "a_H", "a_U", "a_D", "a_DH", "a_DI", "a_DU", "a_O", "a_deltaH",
    "Rt", "tot_I") # prefix a_ represent accumulator variables for incidences

states_to_plot <- c("tot_I", "Rt", state_names[str_detect(state_names, "a_|_curr")], "D") 

p <- ggplot(
      filterstats %>% filter(var %in% states_to_plot) %>% filter(ShortName == 'SG'), aes(x = date)) +
       geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
       geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
       geom_point(data = data %>% filter(ShortName == 'SG') %>% select(-ShortName) %>%
                     mutate(time = dateToYears(date),
                            cases = cases ,       # TODO  / best_params[["epsilon"]][1]
                            deaths = case_when(is.na(deaths) ~ 0,T ~  deaths),
                            cum_deaths = cumsum(deaths)
                            ) %>%
                     rename(
                            a_I = cases,
                       # a_DU = deaths_icu_incid,
                       # a_DH = deaths_noicu_incid,
                       # a_DI = deaths_nohosp_incid,
                      a_H = hosp_incid,
                      #  a_U = icu_incid,
                       #U_curr = icu_curr,
                       H_curr = hosp_curr,
                        D = cum_deaths,
                        a_D = deaths,
                         a_O = discharged,
                        a_deltaH = delta_hosp,
                        a_deltaID = delta_ID
                     ) %>% 
                     gather(var, value, -time, -date),# %>%
                   # mutate(value = as.numeric(factor(value))),
                   aes(y = value)) +
       facet_wrap( ~ var, scales = "free")  +
       theme_bw()
print(p)



a <- data %>% filter(ShortName == 'JU') %>% select(-ShortName) %>%
  mutate(time = dateToYears(date),
         cases = cases  ,       # TODObest_params[["epsilon"]][1]
         deaths = case_when(is.na(deaths) ~ 0,T ~  deaths),
         cum_deaths = cumsum(deaths)
  ) %>%
  rename(
    a_I = cases,
    # a_DU = deaths_icu_incid,
    # a_DH = deaths_noicu_incid,
    # a_DI = deaths_nohosp_incid,
    #a_H = hosp_incid,
   #a_U = icu_incid,
    #U_curr = icu_curr,
    H_curr = hosp_curr,
    D = cum_deaths,
    a_D = deaths,
    a_O = discharged,
    a_deltaH = delta_hosp,
    a_deltaID = delta_ID
  ) %>% 
  gather(var, value, -time, -date) %>%
  mutate(value = as.numeric(factor(value)))
p <- ggplot(a) +  geom_point(aes(x = date, y = value)) +
     facet_wrap( ~ var, scales = "free")  +
        theme_bw()
print(p)
