
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

ggsave(p, filename = glue("results/plot_{suffix}.png"), width = 9, height = 6)
