
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
alpha_npi <- .8
cnt_plots <- lapply(config$places, function(cnt) {
  
  cnt_label <- data.frame(date = as.Date("2020-04-10"), y = ifelse(cnt == "CH", 3.75, 3.5), label = cnt)
  toplot <- filterstats %>% 
    filter(var == "Rt", ShortName %in% cnt) %>% 
    mutate(alpha = case_when(date < date_lim_conf ~ 1,
                             T ~ fadeAlpha(date, date_lim_conf, as.Date("2020-04-15"))),
           alpha = 0.2 * alpha) 
  
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

geodata <- read_csv("data/ch/geodata.csv")

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
  filter(var != "residential") %>% 
  mutate(var = factor(var, labels = var_labels[var_labels != "Residential"])) %>% 
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
         shape =  guide_legend(title = "Lag between changes")) +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_blank())

ggsave(preg, filename = "COVID-pomp/results/figs/mobility_reg_coef.png", width = 6.5, height = 4.5)

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
         shape =  guide_legend(title = "Lag between changes")) +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_blank())

ggsave(pcorr, filename = "COVID-pomp/results/figs/mobility_corr_coef.png", width = 6.5, height = 4.5)


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
  arrange(rollmean) %>% 
  slice(1) %>% 
  select(ShortName, var, rollmean) %>% 
  mutate(rollmean = -rollmean) %>% 
  spread(var, rollmean) %>% 
  magrittr::set_colnames(c("ShortName", "gp", "pa", "re", "rr", "ts", "wp")) %>% 
  inner_join(r0_change) %>% 
  mutate_at(c("median", "q025", "q975"), function(x) 100*x) %>% 
  group_by(ShortName) %>% 
  mutate(sd = getSd(mean, q975, q025)) %>% 
  left_join(select(geodata, ShortName, pop2018)) %>% 
  ungroup() %>% 
  # filter(!(ShortName %in% c("BL", "GE", "CH"))) %>%
  mutate(w = 1/sd^2/sum(1/sd^2))


fit <- lm(median ~ ts
          , data = change_join
          , weights = 1/sd^2) 

summary(fit)

fit <- mgcv::gam(median ~ s(ts)
                 , data = change_join
                 , family = Gamma(link = log)
                 , method = "REML"
                 # , subset = !(ShortName %in% c("GE", "BL", "CH"))
                 # , weights = w
) 

summary(fit)
plot(fit)

p_mc <- change_join %>% 
  gather(var, value, -ShortName, -q025, -q975, -median, -mean,-pop2018, -sd, -w)  %>% 
  mutate(var = factor(var, labels = var_labels),
         alpha = case_when(ShortName == "CH" ~ 0.3, T ~ .5)) %>% 
  ggplot(aes(x = value, y = median, ymin = q025, ymax = q975, alpha = I(alpha))) +
  geom_errorbar(width = 0) +
  geom_label(aes(label = ShortName), size = 2) +
  facet_wrap(~var, scales = "free_x") +
  labs(x = "Percent change in activity", y = "Percent change in R0") +
  theme_bw()

ggsave(p_mc, filename = "COVID-pomp/results/figs/mob_cantons.png", width = 6.5, height = 4.5)


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
