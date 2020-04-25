# function to convert dates to fractions of years for model
dateToYears <- function(date, origin = as.Date("2020-01-01"), yr_offset = 2020) {
  julian(date, origin = origin)/365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateLabel <- function(x, format = "%B-%d") {
  format(yearsToDate(x), format= format)
}
yearsToDateTime <- function(year_frac, origin = as.Date("2020-01-01"), yr_offset = 2020.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}


erlangLL <- function(lambda, x, k) {
  -1 * sum(dgamma(x, shape = k , rate = lambda, log = T))
}

fitErland <- function(x, ks = 1:10, zero_replace = .5) {
  x[x == 0] <- zero_replace
  x <- x[x>0]
  bind_rows(
    lapply(ks, 
           function(k) 
             optimize(erlangLL, x = x, k = k, interval = c(0, 1)) %>% 
             as_tibble() %>% 
             mutate(k = k, 
                    ll = -1 * objective,
                    lambda = minimum,
                    mean = k/lambda,
                    variance = k/(lambda^2)) %>%
             select(k, ll, lambda, mean, variance)
           
    )
  ) %>% 
    arrange(desc(ll)) %>% 
    slice(1)
}

##'Returns a map of configuration loaded from the config YAML
##'@param fname Load configuration from fname (optional, otherwise loads from CONFIG_PATH env var)
##'@example config$parameters_seir$gamma
##'
##'@export
load_config <- function(fname) {
  require(yaml)
  
  if (missing(fname)) {
    fname <- Sys.getenv("CONFIG_PATH")
  }
  if (!missing(fname)) {
    handlers <- list(map=function(x) { class(x) <- "config"; return(x) })
    return(tryCatch(yaml.load_file(fname, handlers=handlers), error = function(e) { stop(paste("Could not find file: ", fname)) }))
  } else {
    return(NA)
  }
}

buildSuffix <- function(name, place, lik_components, params_to_fit, other = NULL){
  param_suffix <- ifelse(is.null(params_to_fit), '', str_c(names(params_to_fit), collapse = '-'))
  
  if (!is.null(other)) {
    other_suffix <- str_c("_", str_c(other, collapse = "_"))
  } else {
    other_suffix <- ""
  }
  suffix <- glue("{name}_{place}_{str_c(lik_components, collapse = '-')}_{param_suffix}{other_suffix}")
}

parseLikelihood <- function(lstring) {
  #' @title Extract likelihood components
  #' @description Parses the likelihood string and verifies that all components exits
  #' @param lstring String passed to Rscript
  #' @return list with components as vector, string of sum and of prod
  
  # Valid likelihoods are the ones that can be parsed by pomp's dmeasure snippet
  valid <- c("deaths" = "d", "delta_hosp" = "deltah", "cases" = "c", "hosp_incid" = "h",
             "delta_icu" = "deltau")
  
  # Parse likelihood components to use
  lik_components <- str_split(lstring, "-")[[1]] 
  
  # Check all valid
  for (l in lik_components)
    if (!(l %in% valid))
      stop(glue("Likelihood component '{l}' not in list of valid strings.
                Please choose among {str_c(str_c(names(valid), valid, sep = ':'), collapse = '|')}."))
  
  lik_sum <- str_c(str_c("ll_", lik_components), collapse = "+")
  lik_prod <- str_c(str_c("ll_", lik_components), collapse = "*")
  
  return(list(components = lik_components, lsum = lik_sum, lprod = lik_prod))
}


getStartEndDates <- function(epidata, ...) {
  #' @title Get start and end dates
  #' @description Gest the start and end dates for the simulations
  #' @param epidate
  #' @return list with start and end dates
  
  start_date <- epidata$date[1] - 5
  end_date <- max(epidata$date + 5)
  
  return(list(start = start_date, end = end_date))
}

ingestData <- function(fdata, likelihood, dt = "day", params = NULL, ...) {
  #' @title Ingest data
  #' @description Ingests data, applies pre-processing, and performs checks
  #' @param fdata data file
  #' @param likelihood likelihood components to use
  #' @param dt time resolution of the data used to fill in missing dates with NA
  #' @param params other parameters to pass to customIngestData
  #' @return tibble with processed data
  
  # Run custom ingestion function
  epidata <- customIngestData(fdata, params)
  
  # Check if all columns are present
  epi_cols <- c("date", "case_incid", "death_incid", "hosp_incid", "discharge_incid", "icu_incid",
                "hosp_curr", "icu_curr", "cum_deaths", "delta_hosp", "delta_ID", "delta_icu",
                "death_noicu_incid", "death_icu_incid", "death_nohosp_incid")
  
  for (ec in epi_cols) {
    if (!(ec %in% colnames(epidata))) {
      if (ec == "date")
        stop("Column 'date' is not in the data, cannot proceed.")
      epidata <- mutate(epidata, !!ec := NA)
      warning(glue("Filling column '{ec}' with NAs!"))
    }}
  
  # Sanity checks
  for (ec in epi_cols[epi_cols != "date" & !str_detect(epi_cols, "delta")]) {
    if (sum(epidata[[ec]] < 0, na.rm = T) > 1) {
      epidata[[ec]][epidata[[ec]] < 0] <- 0
      warning(glue("Found negative numbers in column '{ec}', filling with 0s!"))
    }}
  
  # Check if there is data to compute likelihoods
  lik_dict <- c("c" = "case_incid",
                "d" = "death_incid",
                "h" = "hosp_incid",
                "o" = "discharge_incid",
                "deltah" = "delta_hosp",
                "deltau" = "delta_icu",
                "deltaID" = "delta_ID")
  nlik <- 0
  for (lik in likelihood$components) {
    nobs <- sum(!is.na(epidata[[lik_dict[lik]]]))
    cat("\n", glue("* Using {nobs} observations for {lik_dict[lik]}"))
    nlik <- nlik + (nobs>0)
  }
  
  if (nlik == 0)
    stop(glue("Found no data for likelihood based on {str_c(likelihood$components, collapse = ', ')}"))
  
  # Get start and end dates
  start_end_dates <- getStartEndDates(epidata)
  
  # Add missing rows and fill with NAs (here assumes dt is in days in the data)
  epidata <- epidata %>% 
    tidyr::complete(date = seq.Date(start_end_dates$start, start_end_dates$end, by = dt))  
  
  # Set time column in fraction of years
  epidata$time <- dateToYears(epidata$date)
  
  return(select(epidata, !!c( "time", epi_cols)))
}

setParameters <- function(config, opt, manual = NULL, ...) {
  #' @title Set parameters
  #' @description function to set parameters of the model
  #' @param config Configuration file
  #' @param manual Manual values to set
  #' @return named vector of parameters
  
  # initialize empty parameter vector
  params <- set_names(rep(0, length(param_names)), param_names)
  input_params <- suppressWarnings(unlist(yaml::read_yaml(glue("{opt$b}{config$parameters}"))))
  params[param_fixed_names] <- as.numeric(input_params[param_fixed_names])
  
  # Get population
  params["pop"] <- getPop(config, place)
  
  # Initialize the parameters to estimate (just initial guesses)
  params["std_X"] <- 0 # 1e-4
  params["epsilon"] <- 1
  params["k"] <- 5
  params["I_0"] <- 100 / params["pop"]
  
  # adjust the rate parameters depending on the integration delta time in years (some parameter inputs given in days)
  params[param_rates_in_days_names] <- params[param_rates_in_days_names] * 365.25
  params["R0_0"] <- 3
  
  if (!is.null(manual)) {
    for (p in manual) {
      if (p$param %in% names(params)) {
        cat("changing", p$param)
        params[p$param] <- p$value
      } else {
        warning(glue("Did not find parameter '{p$param}', skipping."))
      }
    }
  }
  
  return(params)
}

getParamBounds <- function(config, params, param_fixed_names, use_case_incid) {
  # lower bound for positive parameter values
  min_param_val <- 1e-5
  
  # define the bounds for the parameters to estimate
  parameter_bounds <- tribble(
    ~param, ~lower, ~upper,
    # Process noise
    "std_X", 2, 3.5
  )
  
  if (use_case_incid) {
    parameter_bounds <- rbind(
      parameter_bounds,
      tribble(
        ~param, ~lower, ~upper,
        # Measurement model
        "k", .1, 10,
        "epsilon", 0.2, 0.5
      ))
  }
  
  if (!is.null(config$parameters_to_fit)) {
    # additional params to fit
    other_bounds <- mapply(
      x = names(config$parameters_to_fit),
      y = config$parameters_to_fit,
      function(x, y) tibble(param = x, lower = y$lower, upper = y$upper),
      SIMPLIFY = F
    ) %>%
      bind_rows()
    parameter_bounds <- rbind(parameter_bounds, other_bounds)
  }
  
  # convert to matrix for ease
  parameter_bounds <- set_rownames(as.matrix(parameter_bounds[, -1]), parameter_bounds[["param"]])
  
  return(parameter_bounds)
}

getInitParams <- function(config, params, param_fixed_names, use_case_incid, npar, profile = F, pars_to_profile = NULL) {
  #' @title Get initial parameters
  #' @description Get initial parameters to run MIF
  #' @param config
  #' @param params
  #' @param param_fixed_names
  #' @param use_case_incid
  #' @param npar
  #' @param profile
  #' @param pars_to_profile
  #' @return return
  
  parameter_bounds <- getParamBounds(config, params, param_fixed_names, use_case_incid)
  
  if (profile) {
    if (is.null(pars_to_profile))
      stop("Please add parameters to profile")
    # TODO add check that all paramaters are known
    parameter_bounds <- parameter_bounds[setdiff(rownames(parameter_bounds), names(pars_to_profile)), ]
    
    # create random vectors of initial parameters given the bounds
    init_params <- do.call(profileDesign,
                           c(pars_to_profile,
                             list(lower = parameter_bounds[, "lower"],
                                  upper = parameter_bounds[, "upper"],
                                  nprof = npar,
                                  type = "sobol"))) 
    
    npar <- nrow(init_params)
    free_pars <- c(rownames(parameter_bounds), names(pars_to_profile))
  } else {
    # create random vectors of initial parameters given the bounds
    init_params <- sobolDesign(
      lower = parameter_bounds[, "lower"],
      upper = parameter_bounds[, "upper"],
      nseq = npar
    )
    free_pars <- rownames(parameter_bounds)
  }
  
  # bind with the fixed valued parameters
  init_params <- cbind(
    init_params,
    matrix(rep(params[!(names(params) %in% free_pars)], each = npar), 
           nrow = npar) %>%
      set_colnames(param_fixed_names)
  )
  
  # Set initial number of infected to proportion
  init_params[, "I_0"] <- init_params[, "I_0"]/init_params[, "pop"]
  
  return(init_params)
}

getRWSD <- function(config, rw.sd_param, tvary, profile = F, pars_to_profile = NULL) {
  
  if (!is.null(config$parameters_to_fit)) {
    other_rw <- lapply(
      names(config$parameters_to_fit),
      function(x) { 
        if (profile) {
          if (x %in% names(pars_to_profile)) {
            glue(", {x} = 0")
          } else {
            glue(", {x} = {rw.sd_param[config$parameters_to_fit[[x]]$rw_sd]}")
          }
        } else {
          glue(", {x} = {rw.sd_param[config$parameters_to_fit[[x]]$rw_sd]}")
        }
      }
    ) %>%
      unlist() %>%
      str_c(collapse = " ")
  } else {
    other_rw <- ""
  }
  
  rw_text <- glue("rw.sd(std_X  =ifelse(time>={tvary}, {rw.sd_param['regular']}, {rw.sd_param['regular']/10}) 
                , k  = {ifelse(use_case_incid, rw.sd_param['regular'], 0)}
                , epsilon   = {ifelse(use_case_incid, rw.sd_param['regular'], 0)}
                {other_rw})")
  
  job_rw.sd <- eval(parse(text = rw_text))
  
  return(job_rw.sd)
}
