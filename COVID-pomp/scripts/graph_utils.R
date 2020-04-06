
dateToYears <- function(date, origin = as.Date("2020-01-01"), yr_offset = 2020) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}

getStates <- function(flist, states) {
  cl <- parallel::makeCluster(8)
  registerDoSNOW(cl)
  
  res <- foreach(fres = flist,
                 .combine = rbind,
                 .packages = c("dplyr", "stringr")) %dopar% {
                   ll_comp <- str_split(fres, "_")[[1]] %>% 
                     .[length(.)] %>% 
                     str_replace("\\.rds", "")
                   res <- readRDS(fres)
                   dplyr::filter(res, res$var %in% states, parset == 1) %>% 
                     mutate(ll_comp = ll_comp)
                 }
  
  stopCluster(cl)
  return(res)
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

# singleR0Reduction <- function(value, dates, dt_left, dt_right) {
#   r0_left <- value[dates < min(dates) + dt_left & !is.na(value)]
#   r0_right <- value[dates > max(dates) - dt_right]
#   return(mean(r0_right)/mean(r0_left)) 
# } 

singleR0value <- function(value, dates, date_left, date_right) {
  r0 <- value[dates >= date_left & dates <= date_right & !is.na(value)]
  return(mean(r0))
} 

computeR0Reduction <- function(filter_data, dt_left, dt_right) {
  filter_data %>% 
    mutate(date = yearsToDate(time)) %>% 
    group_by(ShortName, it, ll_comp) %>% 
    summarise(r0_right = singleR0value(value, date, max(date) - dt_right, max(date)),
              r0_left = singleR0value(value, date, min(date), min(date) + dt_left),
              r0change = r0_right/r0_left) %>%
    ungroup() %>% 
    gather(var, value, -ShortName, -it, -ll_comp) %>% 
    group_by(ShortName, var, ll_comp) %>% 
    summarise(
      mean = mean(value),
      median = median(value),
      q025 = quantile(value, 0.025),
      q975 = quantile(value, 0.975),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75)
    )
}


computeFilterStats <- function(filter_data) {
  filter_stats <- filter_data %>% 
    group_by(time, parset, var, ShortName, ll_comp) %>% 
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


plot_states <- function(filterstats, states_to_plot, cantons = NULL, scales = "fixed") {
  
  if (is.null(cantons)) {
    cantons <- unique(filterstats$ShortName)
  }
  
  p <- filterstats %>% 
    filter(var %in% states_to_plot, ShortName %in% cantons) %>% 
    ggplot(aes(x = date)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = ll_comp), alpha = .2) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = ll_comp), alpha = .2) +
    theme_bw() +
    scale_fill_manual(values = c("#F7A62D", "#9D5EE0")) +
    labs(x = "", y = "Reproduction number") +
    guides(fill = guide_legend(title = "Data used"))
  
  if (length(states_to_plot) == 1 & length(cantons) > 1) {
    if (states_to_plot == "Rt") {
      p <- p + geom_hline(aes(yintercept = 1), lty = 2, size = .3)
    }
    p <- p + facet_wrap(. ~ ShortName, scales = scales)
  } else if (length(states_to_plot) > 1 & length(cantons) > 1) {
    p <- p + facet_grid(var ~ ShortName, scales = scales)
  } else if (length(states_to_plot) >1 & length(cantons) == 1) {
    p <- p + facet_wrap(. ~ var, scales = scales)
  }
  return(p)
}


plot_cnt_all <- function(filterstats, data, states_to_plot, cantons = NULL) {
  
  if (is.null(cantons)) {
    cantons <- unique(data$ShortName)
  }
  
  # Base plot
  p <- plot_states(filterstats, states_to_plot, cantons = cantons, scales = "free_y")
  
  data_names <- c(a_I = "cases",
                  tot_I = "ncumul_conf",
                  a_H = "hosp_incid",
                  H_curr = "hosp_curr",
                  U_curr = "ncumul_ICU",
                  D = "cum_deaths",
                  a_D = "deaths",
                  a_O = "discharged",
                  a_deltaH = "delta_hosp",
                  a_deltaID = "delta_ID")
  
  var_dict <- c("a_I" = "Case incidence",
                "a_H" = "Hospitalization incidence",
                "H_curr" = "Current hospitalizations",
                "D" = "Cumulative deaths",
                "a_D" = "Death incidence",
                "a_O" = "Dicharge from hospital incidence",
                "a_deltaH" = "Daily hospitalization change",
                "a_deltaID" = "Daily hospitalization change with discharges")
  
  # prepare data
  data4plot <- data %>%
    filter(ShortName %in% cantons) %>% 
    rename(!!data_names) %>% 
    gather(var, value, -date, -ShortName) %>% 
    filter(var %in% states_to_plot)
  
  date_breaks <- as.Date(c("2020-03-01", "2020-03-15", "2020-04-01"))
  # Add data
  p <- p +
    geom_point(data = data4plot, aes(y = value), size = .4) +
    scale_x_date(breaks = date_breaks,
                 labels = format(date_breaks, "%B %d")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  return(p)
}
