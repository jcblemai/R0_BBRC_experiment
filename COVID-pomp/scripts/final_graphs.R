
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(glue)
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


npis <- data.frame(date = as.Date(c("2020-02-28", "2020-03-13", "2020-03-16", "2020-03-20")),
                   label = seq(1:4))

cnt_plots <- lapply(config$places, function(cnt) {
  
  cnt_label <- data.frame(date = as.Date("2020-04-10"), y = ifelse(cnt == "CH", 3.75, 3.5), label = cnt)
  toplot <- filterstats %>% 
    filter(var == "Rt", ShortName %in% cnt) %>% 
    mutate(alpha = case_when(date < date_lim_conf ~ 1,
                             T ~ fadeAlpha(date, date_lim_conf, as.Date("2020-04-15"))),
           alpha = 0.2 * alpha) 
  
  p <- toplot %>% 
    ggplot(aes(x = date)) 
  
  alphas <- unique(toplot$alpha)
  if (cnt == "CH") {
    p <- p + geom_label(data = npis, aes(x = date, y = 3.75, label = label), fill = "white")
  }
  for (a in 1:(length(alphas)-1)) {
    p <- p +  
      geom_ribbon(data = filter(toplot, alpha <= alphas[a], alpha >= alphas[a+1]), aes(ymin = q025, ymax = q975), alpha = alphas[a]) +
      geom_ribbon(data = filter(toplot, alpha <= alphas[a], alpha >= alphas[a+1]), aes(ymin = q25, ymax = q75),  alpha = alphas[a])
  }
  p <- p +
    geom_hline(aes(yintercept = 1), lty = 2, size = .3) +
    geom_vline(aes(xintercept = as.Date("2020-02-28")), lty = 3, size = .4, alpha = .6) +
    geom_vline(aes(xintercept = as.Date("2020-03-13")), lty = 3, size = .4, alpha = .6) +
    geom_vline(aes(xintercept = as.Date("2020-03-16")), lty = 3, size = .4, alpha = .6) +
    geom_vline(aes(xintercept = as.Date("2020-03-20")), lty = 3, size = .4, alpha = .6) +
    coord_cartesian(ylim = c(0, 4), xlim = as.Date(c(min(filterstats$date), as.Date("2020-04-18")))) +
    geom_label(data = cnt_label, aes(y = y, label = label), size = ifelse(cnt == "CH", 6, 4)) +
    ylab(ifelse(cnt == "CH", "Basic reproduction number", "")) +
    scale_x_date(date_breaks = ifelse(cnt == "CH", "2 weeks", "1 month"), date_labels = "%B-%d") +
    xlab("") +
    theme_bw()  +
    geom_line(aes(y = median), color = "lightgray", lty = 2, size = .7)
})

lay <- do.call(rbind, 
               list(
                 c(1, 1, 1, 2, 3),
                 c(1, 1, 1, 4, 5),
                 c(6, 7, 8, 9, 10),
                 c(11, 12, 13, 14, NA)))
pr <- arrangeGrob(grobs = cnt_plots, layout_matrix = lay)
ggsave(plot = pr, filename = "COVID-pomp/results/figs/all_R0.png", width = 9, height = 6.5)

# Associations -----------------------------------------------------------------

var_labels <- c("Grocery & Pharmacy", "Parks", "Retail & Recreation", "Transit Stations", "Workplace")

gdata_CH <- read_csv("data/google_mobility_CH.csv") %>% 
  filter(var != "residential") %>% 
  mutate(var = factor(var, labels = var_labels))

r_knot <- filter(filterstats, var  == "Rt") %>% 
  select(date, time, ShortName, median) %>% 
  rename(R0 = median) %>% 
  group_by(ShortName) %>% 
  mutate(date = seq.Date(min(date), by = "1 days", length.out = length(date))) %>% 
  arrange(date) %>% 
  mutate(baseline = mean(R0[!is.na(R0)][1:5]),
         relative = (R0/baseline - 1)*100)



mob_plots <- lapply(config$places, function(cnt) {
  
  cnt_label <- data.frame(date = as.Date("2020-04-10"), y =  ifelse(cnt == "CH", 135, 110), label = cnt)
  toplot <- gdata_CH %>% 
    filter(ShortName %in% cnt)
  if( cnt != "A") {
    guidec <-"none"
  } else {
    guidec <- guide_legend(title = "Mobility type")
  }
  
  p <-  ggplot(toplot, aes(x = date)) +
    geom_hline(aes(yintercept = 0), lty = 2, size = .4) +
    geom_line(aes(y = filled, color = var), alpha = .3) +
    geom_line(aes(y = rollmean, color = var)) +
    geom_vline(aes(xintercept = as.Date("2020-02-28")), lty = 3, size = .4, alpha = .6) +
    geom_vline(aes(xintercept = as.Date("2020-03-13")), lty = 3, size = .4, alpha = .6) +
    geom_vline(aes(xintercept = as.Date("2020-03-16")), lty = 3, size = .4, alpha = .6) +
    geom_vline(aes(xintercept = as.Date("2020-03-20")), lty = 3, size = .4, alpha = .6) +
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
    p <- p + geom_label(data = npis, aes(x = date, y = 135, label = label), fill = "white")
  }
  p
})

prmob <- arrangeGrob(grobs = mob_plots, layout_matrix = lay)
ggsave(plot = prmob, filename = "COVID-pomp/results/figs/all_mobility.png", width = 9, height = 6.5)


crosscorrs <- read_csv("COVID-pomp/results/mobility_cross_correlations.csv")  %>% 
  filter(var != "residential") %>% 
  mutate(var = factor(var, labels = var_labels))

preg <- ggplot(crosscorrs, aes(x = max.coef, y = ShortName, color = lag.type)) + 
  geom_vline(aes(xintercept = 1), lty =3, size = .6) +
  geom_vline(aes(xintercept = 0), lty =2, size = .4) +
  geom_point() +
  geom_errorbarh(aes(xmin = max.coef.low, xmax = max.coef.high), height = 0) +
  facet_wrap(~var, nrow = 1) +
  theme_bw() +
  scale_color_manual(values = c("#920DD9", "#4ACC6F", "#D94F20")) +
  # ggthemes::scale_color_few() +
  labs(x = "Regression coefficient", y = "") +
  guides(color = "none")


pcorr <- ggplot(crosscorrs, aes(x = corr, y = ShortName, color = lag.type)) + 
  geom_vline(aes(xintercept = 1), lty =3, size = .6) +
  geom_vline(aes(xintercept = 0), lty =2, size = .4) +
  geom_point() +
  geom_errorbarh(aes(xmin = corr.low, xmax = corr.high), height = 0) +
  facet_wrap(~var, nrow = 1) +
  theme_bw() +
  scale_color_manual(values = c("#920DD9", "#4ACC6F", "#D94F20")) +
  # ggthemes::scale_color_few() +
  labs(x = "Pearson correlation", y = "")

pc <- arrangeGrob(grobs = pcorr, preg, ncol = 1)
ggsave(pc, filename = "COVID-pomp/results/figs/reg_coef.png", width = 9, height = 6)
