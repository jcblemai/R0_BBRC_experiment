# Summary statistics --------------------------------------------------------------------
library(foreach)
library(itertools)

select <- dplyr::select
computeQuantiles <- function(df) {
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  sim_groups <- cut(df$time, 8)
  
  res <- foreach(sim = isplit(df, sim_groups),
                 .packages = c("tidyverse", "foreach"),
                 .combine = rbind,
                 .inorder = F) %dopar% {
                   
                   sim$value %>% 
                     select(time, scenario_num, scenario_name, NincidInf, NhospCurr, NICUCurr, NincidDeath) %>% 
                     rename(cases = NincidInf,
                            hosp = NhospCurr,
                            icu = NICUCurr,
                            deaths = NincidDeath) %>% 
                     gather(var, value, -time, -scenario_num, -scenario_name) %>% 
                     group_by(time, var, scenario_num, scenario_name) %>% 
                     summarise(mean = mean(value),
                               median = median(value),
                               q025 = quantile(value, 0.025),
                               q975 = quantile(value, 0.975),
                               q25 = quantile(value, 0.25),
                               q75 = quantile(value, 0.75))
                 }
  parallel::stopCluster(cl)
  return(res)
}


computeProbExceedence <- function(df, icu_thresh) {
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  sim_groups <- cut(df$time, 8)
  
  res <- foreach(sim = isplit(df, sim_groups),
                 .packages = c("tidyverse", "foreach"),
                 .combine = rbind,
                 .inorder = F) %dopar% {
                   
                   foreach(thresh = icu_thresh,
                           .combine = rbind) %do% {
                             sim$value %>% 
                               select(time, NICUCurr, scenario_name) %>% 
                               group_by(time, scenario_name) %>% 
                               summarise(n_exceed = sum(NICUCurr > thresh),
                                         exceed_frac = n_exceed/n()) %>% 
                               mutate(icu_thresh = thresh)
                           }
                   
                 }
  parallel::stopCluster(cl)
  return(res)
}


computeStats <- function(df, time_horizons, thLabel) {
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  sim_groups <- cut(df$sim_num, 8)
  
  res <- foreach(sim = isplit(df, sim_groups),
                 .packages = c("tidyverse", "foreach"),
                 .combine = rbind,
                 .inorder = F) %dopar% {
                   
                   foreach(th = time_horizons,
                           .combine = rbind) %do% {
                             sim$value %>%  
                               filter(time < th) %>% 
                               select(time, sim_num, contains("incid"), contains("Curr"), contains("scenario")) %>% 
                               gather(var, value, -time, -sim_num, -scenario_num, -scenario_name) %>% 
                               group_by(var, sim_num, scenario_num, scenario_name) %>% 
                               summarise(max = max(value),
                                         sum = sum(value)) %>% 
                               mutate(th = thLabel(th)) %>% 
                               gather(stat, value, -var, -sim_num, -th, -scenario_num, -scenario_name) %>% 
                               filter((str_detect(var, "Curr") & stat == "max") |
                                        (str_detect(var, "incid") & stat == "sum"))
                           }
                 }
  parallel::stopCluster(cl)
  return(res)
}

thLabel <- function(th) {
  Sys.setlocale("LC_TIME", "C")
  str_replace(format(th, "%B %d"), " ", "")
}

roundValues <- function(x) {
  om <- 10^max(c(1, round(log10(x))-1))
  format(round(x/om)*om, big.mark=",")
} 

aggregateStats <- function(df) {
  df %>% 
    group_by(var, th, stat, scenario_name) %>% 
    summarise(#mean = mean(value),
              median = median(value),
              q025 = quantile(value, 0.025),
              q975 = quantile(value, 0.975)
              #q25 = quantile(value, 0.25),
              #q75 = quantile(value, 0.75)
              ) %>% 
    ungroup() %>% 
    mutate(var = str_c(var, stat, sep = "-")) %>% 
    select(-stat)
}

formaTableContent <- function(df, var_dict, col_order) {
  df %>% 
    mutate(median = map_chr(median, roundValues), 
           iqr = str_c("(", map(q025, roundValues), "-", map(q975, roundValues), ")")) %>% 
    select(-q025, -q975) %>% 
    pivot_wider(names_from = c("th"),
                values_from = c("median", "iqr")) %>% 
    select(var, scenario_name, one_of(col_order)) %>% 
    inner_join(var_dict, .) %>% 
    select(-var)
}

randomDraw <-  function(x, y, n_sim) {
  usims <- sample(unique(y$sim_num), round(n_sim/3), replace = F)
  filter(y, sim_num %in% usims) %>% 
    mutate(sim_num = str_c(sim_num, x))
}


# Function to print results table
printTable <- function(df, scenario, tabnum = 1, col_headers, var_dict, time_horizons, col_order) {
  formaTableContent(df, var_dict, col_order) %>% 
    kable(booktabs = T,
          full_width = T,
          caption = str_c("Table ", tabnum, ". ", scenario),
          align = "c",
          col.names = c("", rep(c("median", "95%CI"), length(time_horizons)))) %>%
    add_header_above(col_headers)   %>%
    kable_styling(latex_options = c("striped", "hold_position"),
                  full_width = T)
  
}



# Time series ----------------------------------------------------------------
make_combined_plot <- function(df, data, timecut = Inf, nrow = 2) {
  cval <- c("#1C86EE", "#EE2C2C", "#8A2BE2", "#EEAD0E")
  lab_dict <- c("cases" = "Case incidence",
                "hosp" = "Current hospitalizations",
                "deaths" = "Death incidence",
                "icu" = "Current ICUs")
  ggplot() +
    geom_ribbon(data = filter(df, time < timecut),
                aes(x = time, ymin = q025, ymax = q975, fill = var), alpha = .2) +
    geom_ribbon(data = filter(df, time < timecut),
                aes(x = time, ymin = q25, ymax = q75, fill = var), alpha = .5) +
    geom_point(data = data, aes(x = date, y = value), size = .5) +
    facet_wrap(var~., scales = "free_y", labeller = labeller(var = lab_dict), nrow = nrow) +
    theme_bw() +
    scale_color_manual(values = cval) +
    scale_fill_manual(values = cval) +
    guides(fill = "none", color = "none") + 
    labs(x = "", y = "")
}



