
customIngestData <- function(fdata, params = NULL) {
  #' @title Custom ingest data
  #' @description Custom pre-processing of data to meet format requirements
  #' @param fdata file to read
  #' @details The output dataframe should have the following column naming convention:
  #' case_incid:  reported case incidence
  #' death_incid: reported death incidence
  #' hosp_incid:  reported hospitalization incidence
  #' discharge_incid: reported discharge from hospital incidence
  #' icu_incid:   reported ICU incidence
  #' hosp_curr:   current hospitalizations
  #' icu_curr:    current ICUs
  #' cum_deaths:  cumulative death incidence
  #' delta_hosp:  differences in current hospitalizations
  #' delta_ID:    difference in current hospitalizations + discharge_incid
  #' delta_icu:   difference in current ICU
  #' death_noicu_incid: reported incidence of hospitalized deaths with no ICU
  #' death_icu_incid:   reported incidence of hospitalized deaths in ICU
  #' death_nohosp_incid: reported incidence of non-hospitalized deaths
  #' 
  #' All columns that are not present will be filled with NAs
  #' @return dataframe with data
  
  
  if (params$place == "CH") {
    # load national-level estimates from Probst's github
    ch_data_dir <- "data/ch/cases/covid19-cases-switzerland/"
    dfiles <- dir(ch_data_dir, pattern = "openzh.csv", full.names = T)
    
    ch_colnames <- list(
      cum_cases = sym("cases"),
      cum_deaths = sym("fatalities"),
      hosp_curr = sym("hospitalized"),
      icu_curr = sym("icu"),
      cum_discharged = sym("released")
    )
    
    epidata <- purrr::map(dfiles, function(x) {
      read_csv(x) %>%
        select(Date, CH) %>%
        magrittr::set_colnames(c("date", str_extract(x, "(?<=19_)[a-z]+(?=_)")))
    }) %>%
      reduce(inner_join) %>%
      rename(!!!ch_colnames) %>%
      arrange(date) %>%
      mutate(
        case_incid = c(cum_cases[1], diff(cum_cases)),
        death_incid = c(cum_deaths[1], diff(cum_deaths)),
        discharge_incid = c(cum_discharged[1], diff(cum_discharged)),
        delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
        delta_icu = c(icu_curr[1], diff(icu_curr)),
        hosp_incid = NA,
        delta_ID = NA
      )
  } else if (params$place == "Vaud") {
    hosp_data <- read_csv("data/VD_hosp_data.csv") %>%
      select(-hosp_curr)
    # epidata <- select(epidata, one_of(setdiff(colnames(epidata), colnames(hosp_data)[-1]))) %>%
    #   left_join(hosp_data)
  } else {
    # REPLACE CODE HERE
    epidata <- read_csv(fdata, col_types = cols()) %>%
      mutate(
        case_incid = c(ncumul_conf[1], diff(ncumul_conf)),
        death_incid = c(ncumul_deceased[1], diff(ncumul_deceased)),
        hosp_incid = new_hosp,
        discharge_incid = c(ncumul_released[1], diff(ncumul_released)),
        icu_curr = current_icu,
        hosp_curr = current_hosp,
        cum_deaths = ncumul_deceased,
        delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
        delta_icu = c(icu_curr[1], diff(icu_curr))
      )
  }
  return(epidata)
}

getPop <- function(config, place, ...) {
  #' @title Get population
  #' @description Custom function to get population
  #' @param Configuration file
  #' @return pop
  
  geodata <- read_csv(config$setup, col_types = cols())
  if (place == "CH") {
    pop <- sum(geodata$pop2018, na.rm = T)
  } else {
    pop <- geodata$pop2018[geodata$ShortName == place]
  }
  return(pop)
}
