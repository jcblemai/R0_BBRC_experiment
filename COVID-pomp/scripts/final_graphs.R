
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(glue)
library(foreach)
library(doSNOW)
library(itertools)
source("COVID-pomp/scripts/graph_utils.R")
source("COVID-pomp/scripts/utils.R")

# Load data --------------------------------------------------------------------

fdata <- list.files(path = "data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2/", pattern = "COVID19_Fallzahlen_Kanton_*", full.names = TRUE)
readCantons <- function(f){
  df <- read_csv(f)
  if(str_detect(f, "NE"))
    df <- mutate(df, date = as.Date(date, "%d.%m.%Y"))
  df
}

data <- purrr::map(fdata, readCantons) %>% bind_rows() %>% 
  # mutate(date = as.Date(date, format = "%d.%m.%Y")) %>% 
  rename(ShortName = abbreviation_canton_and_fl) %>% 
  select(date, ShortName, ncumul_conf, current_hosp, current_icu, ncumul_deceased, ncumul_released)

data <- data %>% arrange(date) %>% group_by(ShortName) %>%
  mutate(cases = c(NA, diff(ncumul_conf)),
         deaths = c(NA, diff(ncumul_deceased)),
         cum_deaths = ncumul_deceased,
         hosp_curr = current_hosp,
         icu_curr = current_icu,
         discharged = c(ncumul_released[1], diff(ncumul_released)),
         delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
         delta_ID = delta_hosp + discharged) %>%
  # select(date, cases, deaths, cum_deaths, hosp_curr, discharged, delta_hosp, delta_ID) %>% 
  mutate(hosp_incid = NA)  %>%
  ungroup()


geodata <- read_csv("data/ch/geodata.csv")

# plot -------------------------------------------------------------------------

config <- load_config("pomp_config.yaml")

#  Time windows
tw_left <- as.Date(config$timewindow_R0$left)
tw_right <- as.Date(config$timewindow_R0$right)
cantons_to_analyze <- config$places

filterstats <- readRDS("COVID-pomp/results/filtered_states_all.rds") %>%
  filter(ShortName %in% config$places)

# Loop over cantons and make plot
fadeAlpha <- function(d, d_start, d_end, alpha_end = 0.4) {
  1 - (1-alpha_end) * as.numeric(difftime(d, d_start, units = "days"))/as.numeric(difftime(d_end, d_start, units = "days"))
}

date_lim_conf <- as.Date("2020-04-05")

# Define NPIS
npis <- data.frame(
  date = as.Date(c("2020-02-28", "2020-03-13", "2020-03-16", "2020-03-20")),
  label = seq(1:4))

npis <- npis %>% mutate(time = dateToYears(date))

alpha_npi <- .8

cnt_plots <- lapply(config$places, function(cnt) {
  
  cnt_label <- data.frame(date = as.Date("2020-04-10"), y = ifelse(cnt == "CH", 3.75, 3.5), label = cnt)
  toplot <- filterstats %>% 
    filter(var == "Rt", ShortName %in% cnt) %>% 
    mutate(alpha = case_when(date < date_lim_conf ~ 1,
                             T ~ fadeAlpha(date, date_lim_conf, as.Date("2020-04-15"))),
           alpha = 0.2 * alpha) %>% 
    mutate(date = seq.Date(min(date), by = "1 days", length.out = length(date)))
  
  p <- toplot %>% 
    ggplot(aes(x = date)) +
    geom_vline(aes(xintercept = as.Date("2020-02-28")), lty = 3, size = .4, alpha = alpha_npi) +
    geom_vline(aes(xintercept = as.Date("2020-03-13")), lty = 3, size = .4, alpha = alpha_npi) +
    geom_vline(aes(xintercept = as.Date("2020-03-16")), lty = 3, size = .4, alpha = alpha_npi) +
    geom_vline(aes(xintercept = as.Date("2020-03-20")), lty = 3, size = .4, alpha = alpha_npi) 
  
  alphas <- unique(toplot$alpha)
  
  for (a in 1:(length(alphas)-1)) {
    p <- p +  
      geom_ribbon(data = filter(toplot, alpha <= alphas[a], alpha >= alphas[a+1]), aes(ymin = q025, ymax = q975), alpha = alphas[a]) +
      geom_ribbon(data = filter(toplot, alpha <= alphas[a], alpha >= alphas[a+1]), aes(ymin = q25, ymax = q75),  alpha = alphas[a])
  }
  p <- p +
    geom_hline(aes(yintercept = 1), lty = 2, size = .3) +
    coord_cartesian(ylim = c(0, 4), xlim = as.Date(c(min(filterstats$date), as.Date("2020-04-18")))) +
    geom_label(data = cnt_label, aes(y = y, label = label), size = ifelse(cnt == "CH", 6, 4)) +
    ylab(ifelse(cnt == "CH", "Basic reproduction number", "")) +
    scale_x_date(date_breaks = ifelse(cnt == "CH", "2 weeks", "1 month"), date_labels = "%B-%d") +
    xlab("") +
    theme_bw()  +
    geom_line(aes(y = median), color = "lightgray", lty = 2, size = .7)
  
  if (cnt == "CH") {
    p <- p + 
      geom_point(data = npis, aes(x = date, y = 3.75), shape = 21,
                 color = "black", fill = "white", size = 7, alpha =1) +
      geom_text(data = npis, aes(x = date, y = 3.75, label = label))
  }
  p
})


lay <- do.call(rbind, 
               list(
                 c(1, 1, 1, 2, 3),
                 c(1, 1, 1, 4, 5),
                 c(6, 7, 8, 9, 10),
                 c(11, 12, 13, 14, NA)))
pr <- arrangeGrob(grobs = cnt_plots, layout_matrix = lay)
ggsave(plot = pr, filename = "COVID-pomp/results/figs/all_R0.png", width = 9, height = 6.5)


# Change points ----------------------------------------------------------------

ffilter <- list.files(path = "COVID-pomp/results/", 
                      pattern = "filtered_", full.names = TRUE) %>% 
  .[str_detect(., "id2o")] #%>% 
# .[str_detect(., param_suffix)] 
# .[str_detect(str_replace_all(., "COVID_CH", ""), str_c(places_to_analyze, collapse = "|"))]

# Get results
sims <- getStates(ffilter, "Rt")


library(mcp)

# Define the model
model1 = list(
  value ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1
)

# Define the model
model2 = list(
  value ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1\
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1
)

# Define the model
model3a = list(
  value ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1\
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 0 + time,
  ~ 1
)


model3b = list(
  value ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1\
  ~ 1,    # joined slope (time_2) at cp_1
  ~ 0 + time,
  ~ 1
)

cl <- parallel::makeCluster(6)
doSNOW::registerDoSNOW(cl)
cp_r0 <- foreach(df = isplit(sims, sims$ShortName),
                 i = icount(1),
                .combine = rbind,
                .packages = c("rjags", "mcp", "tidyverse")
) %dopar% {
  r0_data <- df$value %>% 
    group_by(time) %>% 
    sample_n(1) %>% 
    filter(time > dateToYears(as.Date("2020-02-29")),
           time <  dateToYears(as.Date("2020-04-03")),
           value < 5) %>% 
    mutate(t = row_number()) %>% 
    slice(1:10)
  
  # Fit it. The `ex_demo` dataset is included in mcp
  fit1 <- mcp(model1,  data = r0_data, chains = 4, cores = 1)
  fit2 <- mcp(model2, data = r0_data, chains = 4, cores = 1)
  fit3a <- mcp(model3a, data = r0_data, chains = 4, cores = 1)
  fit3b <- mcp(model3b, data = r0_data, chains = 4, cores = 1)
  tibble(ShortName = df$key[[1]], fits = list(list(fit1, fit2, fit3a, fit3b)))
}

parallel::stopCluster(cl)
# ggplot(r0_data, aes(x = time, y = value)) +
#   geom_point(alpha = 0.3)


p <- plot(fit3c)
p +
  geom_point(data = r0_data %>% group_by(time) %>% summarise(median = median(value)),
             aes(x = time, y = median), col = "red") +
  # geom_vline(aes(xintercept = dateToYears(npis$date[1])), lty = 3, size = .4, alpha = alpha_npi) +
  geom_vline(aes(xintercept = dateToYears(npis$date[2]) + 0.5/365), lty = 3, size = .4, alpha = alpha_npi) +
  geom_vline(aes(xintercept = dateToYears(npis$date[3]) + 0.5/365), lty = 3, size = .4, alpha = alpha_npi) +
  geom_vline(aes(xintercept = dateToYears(npis$date[4]) + 0.5/365), lty = 3, size = .4, alpha = alpha_npi)  +
  scale_x_continuous(labels = yearsToDate)

w1 <- waic(fit1)
w2 <- waic(fit2)
w3 <- waic(fit3)
w3a <- waic(fit3a)

loo1 = loo(fit1)
loo3 = loo(fit3)

loo::loo_compare(loo1, loo3)

# Date crossing 1 --------------------------------------------------------------
r0_reduction <- readRDS("COVID-pomp/results/R0_reduction.rds")

pdate <- r0_reduction %>% 
  filter(var == "t1") %>% 
  inner_join(r0_reduction %>% 
               filter(var == "t1frac") %>% 
               ungroup() %>% 
               select(ShortName, mean) %>% 
               rename(frac = mean)) %>% 
  ungroup() %>% 
  mutate(ShortName = factor(ShortName, ShortName[order(median)])) %>% 
  ggplot(aes(x = ShortName)) +
  geom_point(aes(y = median), size  = 2.5) + 
  geom_errorbar(aes(ymin = q25 , ymax = q75), size = 1.5, width  = 0) +
  geom_errorbar(aes(ymin = q025 , ymax = q975), width  = 0) + 
  coord_flip() +
  scale_y_continuous(breaks = seq(2, 12, by = 2), labels = format(as.Date("2020-03-15")  + seq(2, 12, by = 2), "%B %d")) +
  labs(y = bquote('Date at which '*R[0]*' went bellow 1'), x = "canton") +
  theme_bw()

ggsave(pdate, filename = "COVID-pomp/results/figs/date_crossing.png", width = 5, height = 4)



r0_change <- r0_reduction %>% 
  ungroup() %>% 
  filter(var == "r0change") %>%
  select(ShortName, mean, median, q025, q975) %>%
  mutate_at(c("median", "q025", "q975"), function(x) 1-x)
r0_left <- r0_reduction %>% 
  ungroup() %>% 
  filter(var == "r0_left") %>% 
  select(ShortName, mean, median, q025, q975)

p_change <-r0_change  %>%
  inner_join(r0_left,
             by = "ShortName", suffix = c(".change", ".r0")) %>% 
  # inner_join(select(geodata, ShortName, pop2018)) %>% 
  ggplot(aes(x = median.r0, y = median.change)) +
  geom_point()

# Associations -----------------------------------------------------------------

var_labels <- c("Grocery & Pharmacy", "Parks", "Residential", "Retail & Recreation", "Transit Stations", "Workplace")

gdata_CH <- read_csv("data/google_mobility_CH.csv") %>% 
  # filter(var != "residential") %>% 
  mutate(var = factor(var, labels = var_labels))

r_knot <- filter(filterstats, var  == "Rt") %>% 
  select(date, time, ShortName, median) %>% 
  rename(R0 = median) %>% 
  group_by(ShortName) %>% 
  mutate(date = seq.Date(min(date), by = "1 days", length.out = length(date))) %>% 
  arrange(date) %>% 
  mutate(baseline = mean(R0[!is.na(R0)][1:5]),
         relative = (R0/baseline - 1)*100)



mob_plots <- lapply(c(config$places, "end"), function(cnt) {
  
  cnt_label <- data.frame(date = as.Date("2020-04-10"), y =  ifelse(cnt == "CH", 135, 110), label = cnt)
  
  if( cnt != "end") {
    toplot <- gdata_CH %>% 
      filter(ShortName %in% cnt)
    guidec <-"none"
    p <-  ggplot(toplot, aes(x = date))  +
      geom_vline(aes(xintercept = as.Date("2020-02-28")), lty = 3, size = .4, alpha = alpha_npi) +
      geom_vline(aes(xintercept = as.Date("2020-03-13")), lty = 3, size = .4, alpha = alpha_npi) +
      geom_vline(aes(xintercept = as.Date("2020-03-16")), lty = 3, size = .4, alpha = alpha_npi) +
      geom_vline(aes(xintercept = as.Date("2020-03-20")), lty = 3, size = .4, alpha = alpha_npi) +
      geom_hline(aes(yintercept = 0), lty = 2, size = .4) +
      geom_line(aes(y = filled, color = var), alpha = .3) +
      geom_line(aes(y = rollmean, color = var)) +
      geom_line(data = filter(r_knot, ShortName == cnt, date <= "2020-04-19") , aes(y = relative)) +
      theme_bw() +
      coord_cartesian(ylim = c(-100, 150), xlim = as.Date(c(min(filterstats$date), as.Date("2020-04-18")))) +
      geom_label(data = cnt_label, aes(y = y, label = label), size = ifelse(cnt == "CH", 6, 4)) +
      ylab(ifelse(cnt == "CH", "Percent change", "")) +
      guides(color = guidec) +
      xlab("") +
      ggthemes::scale_color_few()+
      scale_x_date(date_breaks = ifelse(cnt == "CH", "2 weeks", "1 month"), date_labels = "%B-%d") 
    
    if (cnt == "CH") {
      p <- p + 
        geom_point(data = npis, aes(x = date, y = 135), shape = 21,
                   color = "black", fill = "white", size = 7, alpha =1) +
        geom_text(data = npis, aes(x = date, y = 135, label = label))
    }
  } else {
    toplot <- gdata_CH %>% 
      filter(ShortName =="NE")
    p <- ggplot(toplot, aes(x = date))  +
      geom_line(aes(y = filled, color = var), alpha = 0) +
      theme_minimal() +
      labs(x = "" ,y = "") +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            legend.position = c(.5, .61),
            legend.background = element_rect(color = c("#424242"), size = .3),
            legend.key.height = unit(.355, "cm")) + 
      ggthemes::scale_color_few()+
      guides(colour = guide_legend(title = "Mobility type", override.aes = list(alpha = 1)))
  }
  p
})
lay_mob <- lay
lay_mob[4, 5] <- 15
prmob <- arrangeGrob(grobs = mob_plots, layout_matrix = lay_mob)
ggsave(plot = prmob, filename = "COVID-pomp/results/figs/all_mobility.png", width = 9, height = 6.5)


crosscorrs <- read_csv("COVID-pomp/results/mobility_cross_correlations.csv")  %>% 
  # filter(var != "residential") %>% 
  # mutate(var = factor(var, labels = var_labels[var_labels != "Residential"])) %>% 
  mutate(var = factor(var, labels = var_labels)) %>%
  filter(!is.na(lag.type))

preg <- ggplot(crosscorrs, aes(x = max.coef, y = ShortName, color = lag.type, shape = lag.type)) + 
  geom_vline(aes(xintercept = 1), lty =3, size = .6) +
  geom_vline(aes(xintercept = 0), lty =2, size = .4) +
  geom_point() +
  geom_errorbarh(aes(xmin = max.coef.low, xmax = max.coef.high), height = 0) +
  facet_wrap(~var, nrow = 2) +
  theme_bw() +
  scale_color_manual(values = c("#920DD9", "#4ACC6F", "#D94F20")) +
  # ggthemes::scale_color_few() +
  labs(x = "Regression coefficient", y = "") +
  guides(color = guide_legend(title = "Lag between changes"),
         shape =  guide_legend(title = "Lag between changes")) #+
# theme(legend.position = c(0.85, 0.25),
#       legend.background = element_blank())

ggsave(preg, filename = "COVID-pomp/results/figs/mobility_reg_coef.png", width = 8.5, height = 4.5)

pcorr <- ggplot(crosscorrs, aes(x = corr, y = ShortName, color = lag.type, shape = lag.type)) + 
  # geom_vline(aes(xintercept = 1), lty =3, size = .6) +
  geom_vline(aes(xintercept = 0), lty =2, size = .4) +
  geom_point(aes()) +
  geom_errorbarh(aes(xmin = corr.low, xmax = corr.high), height = 0) +
  facet_wrap(~var, nrow = 2) +
  theme_bw() +
  scale_color_manual(values = c("#920DD9", "#4ACC6F", "#D94F20")) +
  labs(x = "Pearson correlation", y = "")  +
  guides(color = "none", shape = "none")   +
  guides(color = guide_legend(title = "Lag between changes"),
         shape =  guide_legend(title = "Lag between changes")) #+
# theme(legend.position = c(0.85, 0.25),
#       legend.background = element_blank())

ggsave(pcorr, filename = "COVID-pomp/results/figs/mobility_corr_coef.png", width = 8.5, height = 4.5)


# pc <- arrangeGrob(pcorr, preg, ncol = 2, layout_matrix = rbind(c(1, 1, 2, 2, 2)))
# ggsave(pc, filename = "COVID-pomp/results/figs/reg_coef.png", width = 6, height = 10)
# Function
getSd <- function(mu, q025, q975) {
  sd1 <- (mu - q025)/2
  sd2 <- (q975 - mu)/2
  return(mean(c(sd1, sd2)))
}

change_join <- gdata_CH %>% 
  group_by(ShortName, var) %>% 
  filter((var != "Residential" & rollmean < 0) |
           (var == "Residential" & rollmean > 0)) %>% 
  arrange(-abs(rollmean)) %>% 
  slice(1) %>% 
  select(ShortName, var, rollmean) %>% 
  # mutate(rollmean = -rollmean) %>% 
  spread(var, rollmean) %>% 
  magrittr::set_colnames(c("ShortName", "gp", "pa", "re", "rr", "ts", "wp")) %>% 
  inner_join(r0_change) %>% 
  mutate_at(c("median", "q025", "q975"), function(x) 100*x) %>% 
  group_by(ShortName) %>% 
  mutate(sd = getSd(mean, q975, q025)) %>% 
  left_join(select(geodata, ShortName, pop2018)) %>% 
  ungroup() %>% 
  filter(ShortName != "CH") %>% 
  # filter(!(ShortName %in% c("BL", "GE"))) %>%
  mutate(w = 1/sd^2/sum(1/sd^2))

getCorr <- function(var) {
  ind <- !is.na( change_join[[var]])
  pho <- cor(change_join$median[ind], change_join[[var]][ind])
  ci <- psychometric::CIr(pho, sum(ind))
  data.frame(cor = pho, cor.low = ci[1], cor.high = ci[2], var = var)
}

correlations <- map_df(c("gp", "pa", "re", "rr", "ts", "wp"),  getCorr)
correlations


p_mc <- change_join %>% 
  gather(var, value, -ShortName, -q025, -q975, -median, -mean,-pop2018, -sd, -w)  %>% 
  mutate(var = factor(var, labels = var_labels),
         alpha = case_when(ShortName == "CH" ~ 0.3, T ~ .5)) %>% 
  ggplot(aes(x = value, y = -median, ymin = -q025, ymax = -q975, alpha = I(alpha))) +
  geom_errorbar(width = 0) +
  geom_label(aes(label = ShortName), size = 2) +
  facet_wrap(~var, scales = "free_x") +
  labs(x = "Percent change in activity", y = "Percent change in R0") +
  theme_bw()

ggsave(p_mc, filename = "COVID-pomp/results/figs/mob_cantons.png", width = 6.5, height = 4.5)


# Changepoints mobility --------------------------------------------------------


# Define the model
modelmob = list(
  rollmean ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1
)

# Define the model
model3 = list(
  rollmean ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1\
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1
)


extractPosteriors <- function(fit) {
  map_df(as.list(fit$mcmc_post), ~select(as_tibble(.), contains("cp")))
}

getChangepoints <- function(df, var, model) {
  df <- df %>%
    filter(var == var)
  fit <- mcp(model, data = df, chains = 4, cores = 1)
  cp_post <- extractPosteriors(fit) %>% 
    gather(cp, value)
  
  cptests <- map_df(npis$time, ~hypothesis(fit, glue("cp_1 < {.}"))) %>% 
    mutate(npi_id = npis$label) %>% 
    mutate(var = str_c("prob_", npi_id)) %>% 
    select(var, p) %>% 
    spread(var, p) %>% 
    mutate(cp = "cp_1")
  
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

cl <- parallel::makeCluster(6)
doSNOW::registerDoSNOW(cl)
change_points <- foreach(cnt = unique(gdata_CH$ShortName),
                         .combine = rbind,
                         .packages = c("mcp", "rjags", "tidyverse", "glue")) %dopar% {
                           mdf <- filter(gdata_CH, ShortName == cnt)
                           mdf <- mdf %>%
                             mutate(time = dateToYears(date)) %>% 
                             filter(date < "2020-04-01", date > "2020-02-28") %>% 
                             filter(!is.na(filled))
                           
                           lapply(unique(mdf$var), 
                                  function(cvar) {
                                    getChangepoints(mdf, cvar, modelmob) %>% 
                                      mutate(var = cvar)
                                  }) %>% 
                             bind_rows() %>% 
                             mutate(ShortName = cnt) 
                         }
parallel::stopCluster(cl)

write_csv(change_points, path = "COVID-pomp/results/mob_changepoints.csv")

cp_labels <- c(
  cp_1 = "Start decrease",
  cp_2 = "Stabilize"
)

yearsToFormat <- function(x) {
  yearsToDate(x) %>% format("%B-%d")
}

p_change <- change_points %>% 
  mutate(cp = factor(cp, labels = cp_labels)) %>% 
  ggplot(aes(x = median, xmin = q025, xmax = q975, y = var, color = cp)) +
  geom_errorbarh(height = 0, alpha = .3) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), size = 1.3, height = 0) +
  geom_point(aes(color = cp, shape = cp), size = 3) +
  facet_wrap(~ShortName) +
  scale_x_continuous(labels = yearsToFormat) + 
  theme_bw() +
  scale_color_manual(values = c("#2BD9B0", "#DA33E6")) +
  guides(color = guide_legend(title = "Activity change"),
         shape = guide_legend(title = "Activity change"))
p_change
ggsave(p_change, filename = "COVID-pomp/results/figs/mob_change_dates.png", width = 10, height = 4)


prob_labels <- c(
  prob_1 = as.character(npis$date[1]),
  prob_2 =  as.character(npis$date[2]),
  prob_3 =  as.character(npis$date[3]),
  prob_4 =  as.character(npis$date[4])
)

p_prob <- change_points %>% 
  filter(cp == "cp_1") %>%
  select(ShortName, contains("prob_"), var) %>% 
  gather(prob, value, -ShortName, -var) %>% 
  ggplot(aes(x = ShortName, y = value)) +
  geom_bar(stat = "identity", fill = "darkgray") +
  facet_grid(var ~ prob, labeller = labeller(prob = prob_labels)) +
  coord_flip() +
  theme_bw() +
  labs(y = "", x = "")
ggsave(p_prob, filename = "COVID-pomp/results/figs/mob_change_prob.png", 
       width = 10, height = 10)

change_points %>%
  filter(var %in% c("Retail & Recreation", "Transit Stations", "Workplace"), cp =="cp_1") %>% 
  group_by(ShortName) %>% 
  summarise_at(vars(contains("prob")), mean)

# Seroprev ---------------------------------------------------------------------

seroprev <- filterstats %>% 
  filter(var == "tot_I") %>% 
  group_by(ShortName) %>% 
  arrange(desc(time)) %>% 
  slice(1) %>% 
  left_join(select(geodata, ShortName, pop2018))

seroprev$pop2018[is.na(seroprev$pop2018)] <- sum(geodata$pop2018)

seroprev %>% 
  mutate(median_prev = median/pop2018, 
         q025_prev = q025/pop2018,
         q975_prev = q975/pop2018) %>% 
  select(ShortName, contains("prev"))

ch_sf <- sf::st_read("data/ch/shp/ch.shp") %>% 
  left_join(seroprev %>% 
              mutate(median_prev = median/pop2018, 
                     q025_prev = q025/pop2018,
                     q975_prev = q975/pop2018) %>% 
              select(ShortName, contains("prev")))
whole_ch <-  summarise(ch_sf, n = n())
p_prev <- ggplot(ch_sf) +
  geom_sf(data = filter(ch_sf, ShortName != "TI"), aes(fill = median_prev*100),
          color = "black", size = .2) +
  geom_sf(data = filter(ch_sf, ShortName == "TI"), fill = "darkred", color ="black", size = .2) +
  geom_sf(data =whole_ch, color ="black", size = .3 ,alpha = 0) +
  ggthemes::theme_map() +
  scale_fill_gradient(low = "#FDCFCE", high = "#C21F14") +
  theme(legend.position = "bottom", legend.key.height = unit(.3, "cm")) +
  guides(fill = guide_colorbar(title = "Sero-Prevalence [%]", title.position = "top"))

ggsave(p_prev, filename = "COVID-pomp/results/figs/seroprev.png", width = 6, height = 5)
