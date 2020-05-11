# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(glue)
library(foreach)
library(doSNOW)
library(itertools)
library(mcp)
library(knitr)
source("COVID-pomp/scripts/graph_utils.R")
source("COVID-pomp/scripts/utils.R")

config <- load_config("pomp_config.yaml")

# Functions
extractPosteriors <- function(fit) {
  map_df(as.list(fit$mcmc_post), ~select(as_tibble(.), contains("cp")))
}

fitMCAP <- function(df, var, model) {
  df <- df %>%
    filter(var == var)
  fit <- mcp(model, data = df, chains = 4, cores = 1)
  return(fit)
}

getChangepoints <- function(fit) {
  
  cp_post <- extractPosteriors(fit) %>% 
    gather(cp, value)
  
  cptests <- lapply(grep("cp_", fit$pars$population, value = T),
                    function(cp) {
                      map_df(npis$time, ~hypothesis(fit, glue("{cp} < {.}"))) %>% 
                        mutate(npi_id = npis$label) %>% 
                        mutate(var = str_c("prob_", npi_id)) %>% 
                        select(var, p) %>% 
                        spread(var, p) %>% 
                        mutate(cp = cp)}) %>% 
    bind_rows()
  
  res <- cp_post %>% 
    group_by(cp) %>% 
    summarise(q025 = quantile(value, 0.025),
              q975 = quantile(value, 0.975),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75),
              median = quantile(value, 0.5),
              mean = mean(value)) %>% 
    left_join(cptests)
  res
}

computeProbGreater <- function(x, y) {
  lapply(y, function(y1) sum(x < y1)/length(x)) %>% 
    unlist() %>% 
    mean()
}

# Varibales
# Define NPIS
npis <- data.frame(
  date = as.Date(c("2020-02-28", "2020-03-13", "2020-03-16", "2020-03-20")),
  label = seq(1:4)) %>% 
  mutate(time = dateToYears(date))

# Change points R0 ----------------------------------------------------------------

fit_cantons <- "COVID-pomp/results/changepoints_R0.csv"
model_comp_file <- "COVID-pomp/results/model_comp_changepoints_R0.csv"

if(!file.exists(fit_cantons)) {
  
  ffilter <- list.files(path = "COVID-pomp/results/", 
                        pattern = "filtered_COVID", full.names = TRUE)
  
  # Get results
  sims <- getStates(ffilter, "Rt")
  
  # Define the models
  model1 <- list(value ~ 1, ~ 0 + time, ~ 1)    # one slope
  model2 <- list(value ~ 1, ~ 0 + time, ~ 0 + time, ~ 1)  # two slopes
  model3 <- list(value ~ 1, ~ 0 + time, ~ 0 + time, ~ 0 + time, ~ 1)  # three slopes
  
  n_sample <- 250    # number of samples from the smoothing distribution of R0 to use
  nchains <- 4     # number of MCMC chains to run
  niter <- 1000    # number of MCMC posterior draws to do
  
  r0_change_points <- foreach(df = isplit(sims, sims$ShortName),
                              .combine = rbind,
                              .packages = c("rjags", "mcp", "tidyverse")
  ) %do% {
    
    r0_data <- df$value %>% 
      group_by(time) %>% 
      sample_n(n_sample) %>% 
      filter(time > dateToYears(as.Date("2020-02-25")),
             time <  dateToYears(as.Date("2020-04-03")),
             value < 5)
    
    mfile <- glue("COVID-pomp/results/chgpnts/changepoint_models_{df$key[[1]]}.rda")
    
    if(file.exists(mfile)) {
      load(mfile)
    } else {
      fit1 <- mcp(model1,  data = r0_data, chains = nchains, cores = nchains, iter = niter)
      fit2 <- mcp(model2,  data = r0_data, chains = nchains, cores = nchains, iter = niter)
      fit3 <- mcp(model3,  data = r0_data, chains = nchains, cores = nchains, iter = niter)
      save(fit1, fit2, fit3, file = mfile)
    }
    fits <- list(fit1, fit2, fit3)
    loos <- lapply(fits, loo)
    chnpts <- mapply(fit = fits,
                     loo_fit = loos,
                     md = paste0("model", 1:3),
                     function(fit, loo_fit, md) {
                       getChangepoints(fit) %>% 
                         mutate(loo = loo_fit$estimates[3, 1],
                                loo.se = loo_fit$estimates[3, 2],
                                model = md)
                     },
                     SIMPLIFY = F) %>% 
      bind_rows() %>% 
      mutate(ShortName = df$key[[1]])
    
    loo_comp <- loo::loo_compare(loos)
    loo_comp <- as.data.frame(loo_comp) %>% 
      mutate(model = rownames(loo_comp),
             ShortName = df$key[[1]])
    
    write_csv(loo_comp, path = model_comp_file, append = file.exists(model_comp_file))
    
    return(chnpts)
  }
  
  # Write results to file
  write_csv(r0_change_points, path = fit_cantons)
  
} else {
  r0_change_points <- read_csv(fit_cantons)
}

# Write changepoint table to latex format
library(knitr)
model_comp <- read_csv(model_comp_file) %>% 
  mutate(diff_sign = abs(elpd_diff) > 5 * se_diff) %>% 
  group_by(ShortName) %>% 
  mutate(diff_sign = case_when(row_number() == 1 ~ "", T ~ as.character(diff_sign))) %>% 
  select(ShortName, model, elpd_diff, se_diff, diff_sign)

# Write latex table
kable(model_comp,
      format = "latex",
      booktabs = T,
      linesep =  c('', '', '\\addlinespace'),
      col.names = c("Canton", "Model", "LOO difference", "SE", "significant")) %>% 
  write(file = "COVID-pomp/results/changepoint_comp_table.tex")

# Change points mobility --------------------------------------------------------

var_labels <- c("Grocery & Pharmacy", 
                "Parks",
                "Residential",
                "Retail & Recreation",
                "Transit Stations",
                "Workplace")

gdata_CH <- read_csv("data/google_mobility_CH.csv") %>% 
  mutate(var = factor(var, labels = var_labels))


# Define the model
modelmob = list(
  rollmean ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1
)

mobility_change_points <- foreach(cnt = unique(gdata_CH$ShortName),
                                  .combine = rbind,
                                  .packages = c("mcp", "rjags", "tidyverse", "glue")
) %do% {
  mdf <- filter(gdata_CH, ShortName == cnt)
  mdf <- mdf %>%
    mutate(time = dateToYears(date)) %>% 
    filter(date < "2020-04-01", date > "2020-02-28") %>% 
    filter(!is.na(filled))
  # load chp 1
  load(glue("COVID-pomp/results/chgpnts/changepoint_models_{cnt}.rda"))
  
  cp1_r0 <- purrr::map(as.list(fit1$mcmc_post), ~.[,"cp_1"]) %>% 
    unlist()
  
  lapply(unique(mdf$var), 
         function(cvar) {
           fit <- fitMCAP(mdf, cvar, modelmob)
           cp1_var <- purrr::map(as.list(fit$mcmc_post), ~.[,"cp_1"]) %>% 
             unlist()
           
           pvar_after_R0 <- computeProbGreater(cp1_r0, cp1_var)
           
           getChangepoints(fit) %>% 
             mutate(var = cvar,
                    pvar_after_R0 = pvar_after_R0) 
         }) %>% 
    bind_rows() %>% 
    mutate(ShortName = cnt) 
}

write_csv(mobility_change_points, path = "COVID-pomp/results/changepoints_mobility.csv")
