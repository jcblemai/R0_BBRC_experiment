
dateToYears <- function(date, origin = as.Date("2020-01-01"), yr_offset = 2020) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}

##'Function to plot county final sizes
##'
##'@param country_final_dat sf object with county final sizes
##'
##'@return ggplot object, final sizes with CI (y) by county (x)
##'
make_map_plot <- function(country_data){
  p <- ggplot(country_data,
              aes(x=reorder(county.name,-median), y=median,ymin=pi_low, ymax=pi_high )) +
    geom_pointrange() +
    ylab("Infections") +
    xlab("County") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}

compute_filterstats <- function(filter_data) {
  filter_stats <- filter_data %>% 
    group_by(time, parset, var, ShortName) %>% 
    summarise(mean = mean(value, na.rm = T),
            q025 = quantile(value, 0.025, na.rm = T),
            q975 = quantile(value, 0.975, na.rm = T),
            q25 = quantile(value, 0.25, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T)) %>%
    ungroup() %>% 
    mutate(date = yearsToDate(time),
         parset = factor(parset))
  return(filter_stats)

}


plot_states <- function(filter_stats, data, states_to_plot) {
  p <- ggplot(filter_stats %>% filter(var %in% states_to_plot), aes(x = date)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
    # geom_line(aes(y = mean)) +
    # geom_point(data = data %>%
    #              mutate(time = dateToYears(date),
    #                     #cases = cases / best_params[["epsilon"]][1],       # TODO
    #                     deaths = case_when(is.na(deaths) ~ 0,T ~  deaths),
    #                     cum_deaths = cumsum(deaths)) %>%
    #              rename(
    #                a_I = cases,
    #                # a_DU = deaths_icu_incid,
    #                # a_DH = deaths_noicu_incid,
    #                # a_DI = deaths_nohosp_incid,
    #                #a_H = hosp_incid,
    #                # a_U = icu_incid,
    #                # U_curr = icu_curr,
    #                H_curr = hosp_curr,
    #                D = cum_deaths,
    #                a_D = deaths,
    #                a_O = discharged,
    #                a_deltaH = delta_hosp,
    #                a_deltaID = delta_ID
    #              ) %>% 
    #              gather(var, value, -time, -date),
    #            aes(y = value)) +
    facet_wrap(ShortName ~ var, scales = "free")  +
    theme_bw()
  return(p)
}


plot_cnt_all <- function(filterstats, data, ctn) {
  p <- ggplot(
    filterstats %>% filter(var %in% states_to_plot) %>% filter(ShortName == ctn), aes(x = date)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
    geom_point(data = data %>% filter(ShortName == ctn) %>% select(-ShortName) %>%
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
  return(p)
}