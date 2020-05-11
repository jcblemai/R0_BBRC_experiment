
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(glue)
library(foreach)
library(doSNOW)
library(itertools)
source("COVID-pomp/scripts/graph_utils.R")
source("COVID-pomp/scripts/utils.R")

Sys.setlocale("LC_TIME", "C")

# Function
getSd <- function(mu, q025, q975) {
  sd1 <- (mu - q025)/2
  sd2 <- (q975 - mu)/2
  return(mean(c(sd1, sd2)))
}

fadeAlpha <- function(d, d_start, d_end, alpha_end = 0.1) {
  1 - (1-alpha_end) * as.numeric(difftime(d, d_start, units = "days"))/as.numeric(difftime(d_end, d_start, units = "days"))
}

readCantons <- function(f){
  df <- read_csv(f)
  if(str_detect(f, "NE"))
    df <- mutate(df, date = as.Date(date, "%d.%m.%Y"))
  df
}

getCorr <- function(var) {
  ind <- !is.na( change_join[[var]])
  pho <- cor(change_join$median[ind], change_join[[var]][ind])
  ci <- psychometric::CIr(pho, sum(ind))
  data.frame(cor = pho, cor.low = ci[1], cor.high = ci[2], var = var)
}

yearsToFormat <- function(x) {
  yearsToDate(x) %>% format("%B-%d")
}

resString <- function(df, title) {
  glue("{title}\nCH: {format(df$median[df$ShortName == 'CH'], digits = 2)}",
       " ({format(df$q025[df$ShortName == 'CH'], digits = 2)} - {format(df$q975[df$ShortName == 'CH'], digits = 2)})\n",
       "Range medians: {format(min(df$median), digits = 2)} - {format(max(df$median), digits = 2)}")
}

# Load data --------------------------------------------------------------------

fdata <- list.files(path = "data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2/",
                    pattern = "COVID19_Fallzahlen_Kanton_*", full.names = TRUE)

data <- purrr::map(fdata, readCantons) %>% 
  bind_rows() %>% 
  rename(ShortName = abbreviation_canton_and_fl) %>% 
  select(date, ShortName, ncumul_conf, current_hosp,
         current_icu, ncumul_deceased, ncumul_released)

data <- data %>% arrange(date) %>% group_by(ShortName) %>%
  mutate(cases = c(NA, diff(ncumul_conf)),
         deaths = c(NA, diff(ncumul_deceased)),
         cum_deaths = ncumul_deceased,
         hosp_curr = current_hosp,
         icu_curr = current_icu,
         discharged = c(ncumul_released[1], diff(ncumul_released)),
         delta_hosp = c(hosp_curr[1], diff(hosp_curr)),
         delta_ID = delta_hosp + discharged) %>%
  mutate(hosp_incid = NA)  %>%
  ungroup()

# Plot of R0s -------------------------------------------------------------------------
config <- load_config("pomp_config.yaml")
geodata <- read_csv(config$setup)

canton_dict <- geodata$Name
names(canton_dict) <- geodata$ShortName
canton_dict <- c(canton_dict, CH = "Switzerland")

#  Time windows
tw_left <- as.Date(config$timewindow_R0$left)
tw_right <- as.Date(config$timewindow_R0$right)
cantons_to_analyze <- config$places

# load filtered varibles
filterstats <- readRDS("COVID-pomp/results/filtered_states_all.rds") %>%
  filter(ShortName %in% config$places)

# Limit of date at which values of R0 are considered confident
date_lim_conf <- as.Date("2020-04-10")

# Define NPIS
npis <- data.frame(
  date = as.Date(c("2020-02-28", "2020-03-13", "2020-03-16", "2020-03-20")),
  label = seq(1:4)) %>%
  mutate(time = dateToYears(date))

alpha_npi <- .8
places <- sort(unique(filterstats$ShortName))
places <- c("CH", places[places != "CH"])

cnt_plots <- lapply(places, function(cnt) {
  
  cnt_label <- data.frame(date = as.Date("2020-04-25"),
                          y = ifelse(cnt == "CH", 3.75, 3.5), label = canton_dict[cnt])
  toplot <- filterstats %>% 
    filter(var == "Rt", ShortName %in% cnt) %>% 
    mutate(alpha = case_when(date < date_lim_conf ~ 1,
                             T ~ fadeAlpha(date, date_lim_conf, as.Date("2020-04-25"))),
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
    coord_cartesian(ylim = c(0, 4), xlim = as.Date(c(min(filterstats$date), as.Date("2020-04-25")))) +
    geom_label(data = cnt_label, aes(y = y, label = label), 
               size = ifelse(cnt == "CH", 4, 2.5), hjust = 1, nudge_x = -1) +
    ylab(ifelse(cnt == "CH", "Basic reproduction number", "")) +
    scale_x_date(date_breaks = ifelse(cnt == "CH", "2 weeks", "1 month"), date_labels = "%B-%d") +
    xlab("") +
    theme_bw()  +
    geom_line(aes(y = median), color = "white", lty = 2, size = ifelse(cnt == "CH", .5, .3))
  
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
                 c(1, 1, 1, 2, 3, 4),
                 c(1, 1, 1, 5, 6, 7),
                 c(8, 9, 10, 11, 12, 13)))
pr <- arrangeGrob(grobs = cnt_plots, layout_matrix = lay)
ggsave(plot = pr, filename = "COVID-pomp/results/figs/all_R0.png", width = 11, height = 4.5)


# Date crossing 1 --------------------------------------------------------------
r0_reduction <- readRDS("COVID-pomp/results/R0_reduction.rds")
date_start = as.Date("2020-03-05")
pdate <-  r0_reduction %>% 
  filter(var == "t1") %>% 
  inner_join(r0_reduction %>% 
               filter(var == "t1frac") %>% 
               ungroup() %>% 
               select(ShortName, mean) %>% 
               rename(frac = mean)) %>% 
  ggplot(aes(x = ShortName)) +
  geom_point(aes(y = median), size  = 2.5) + 
  geom_errorbar(aes(ymin = q25 , ymax = q75), size = 1.5, width  = 0) +
  geom_errorbar(aes(ymin = q025 , ymax = q975), width  = 0) + 
  coord_flip() +
  scale_y_continuous(breaks = seq(1, 21, by = 3), labels = format(date_start  + seq(1, 21, by = 3), "%B %d")) +
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
r0_time <- r0_reduction %>% 
  ungroup() %>% 
  filter(var == "t1") %>% 
  select(ShortName, mean, median, q025, q975)
r0_right <- r0_reduction %>% 
  ungroup() %>% 
  filter(var == "r0_right") %>% 
  select(ShortName, mean, median, q025, q975)
p_change <-r0_change  %>%
  inner_join(r0_left,
             by = "ShortName", suffix = c(".change", ".r0")) %>% 
  ggplot(aes(x = median.r0, y = median.change)) +
  geom_point()

# Associations with mobility -----------------------------------------------------------------

var_labels <- c("Grocery & Pharmacy", "Parks",
                "Residential", "Retail & Recreation",
                "Transit Stations", "Workplace")

gdata_CH <- read_csv("data/google_mobility_CH.csv") %>% 
  mutate(var = factor(var, labels = var_labels))

r_knot <- filter(filterstats, var  == "Rt") %>% 
  select(date, time, ShortName, median) %>% 
  rename(R0 = median) %>% 
  group_by(ShortName) %>% 
  mutate(date = seq.Date(min(date), by = "1 days", length.out = length(date))) %>% 
  arrange(date) %>% 
  mutate(baseline = mean(R0[!is.na(R0)][1:5]),
         relative = (R0/baseline - 1)*100)

mob_plots <- lapply(places, function(cnt) {
  
  cnt_label <- data.frame(date = as.Date(ifelse(cnt == "CH", "2020-03-03", "2020-04-18")), 
                          y =  ifelse(cnt == "CH", -75, 110), 
                          label = canton_dict[cnt])
  toplot <- gdata_CH %>% 
    filter(ShortName %in% cnt)
  
  if(cnt == "CH") {
    guidec <- guide_legend(title = "Activity type", nrow = 3)
  }  else {
    guidec <-"none" 
  }
  
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
    geom_label(data = cnt_label, aes(y = y, label = label), 
               size = ifelse(cnt == "CH", 4, 2.5), hjust = 1, nudge_x = -1) +
    ylab(ifelse(cnt == "CH", "Percent change", "")) +
    guides(color = guidec) +
    xlab("") +
    ggthemes::scale_color_few()+
    scale_x_date(date_breaks = ifelse(cnt == "CH", "2 weeks", "1 month"), date_labels = "%B-%d") 
  
  if (cnt == "CH") {
    p <- p + 
      geom_point(data = npis, aes(x = date, y = 135), shape = 21,
                 color = "black", fill = "white", size = 7, alpha =1) +
      geom_text(data = npis, aes(x = date, y = 135, label = label)) +
      theme(
        legend.position = c(0.99, .975),
        legend.background = element_rect(color = c("#424242"), size = .2),
        legend.key.height = unit(.355, "cm"),
        legend.key.width = unit(.25, "cm"),
        legend.spacing.x = unit(0.075, "cm"),
        legend.text.align = 0, 
        legend.justification = c(1, 1),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9))
  }
  
  p
})

lay_mob <- lay
prmob <- arrangeGrob(grobs = mob_plots, layout_matrix = lay_mob)

# FIGURE 2
ggsave(plot = prmob, filename = "COVID-pomp/results/figs/all_mobility.png", width = 11, height = 4.5)


crosscorrs <- read_csv("COVID-pomp/results/mobility_cross_correlations.csv")  %>% 
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
  labs(x = "Regression coefficient", y = "") +
  guides(color = guide_legend(title = "Lag between changes"),
         shape =  guide_legend(title = "Lag between changes"))

ggsave(preg, filename = "COVID-pomp/results/figs/mobility_reg_coef.png", width = 8.5, height = 4.5)

corr_title <- "Lag of maximum\ncross-correlation"
pcorr <- ggplot(crosscorrs, aes(x = corr, y = ShortName, fill = lag), color = "black") + 
  geom_vline(aes(xintercept = 0), lty =2, size = .4) +
  geom_errorbarh(aes(xmin = corr.low, xmax = corr.high), height = 0) +
  geom_point(pch = 21, alpha = 1, size = 2.5) +
  facet_wrap(~var, nrow = 2) +
  theme_bw() +
  scale_fill_gradient2() +
  scale_color_gradient2() +
  labs(x = "Cross-correlation", y = "")  +
  guides(fill = guide_legend(title = corr_title)) 

ggsave(pcorr, filename = "COVID-pomp/results/figs/mobility_corr_coef.png", width = 8.5, height = 4.5)

change_join <- gdata_CH %>% 
  group_by(ShortName, var) %>% 
  filter((var != "Residential" & rollmean < 0) |
           (var == "Residential" & rollmean > 0)) %>% 
  arrange(-abs(rollmean)) %>% 
  slice(1) %>% 
  select(ShortName, var, rollmean) %>% 
  spread(var, rollmean) %>% 
  magrittr::set_colnames(c("ShortName", "gp", "pa", "re", "rr", "ts", "wp")) %>% 
  inner_join(r0_change) %>% 
  mutate_at(c("median", "q025", "q975"), function(x) 100*x) %>% 
  group_by(ShortName) %>% 
  mutate(sd = getSd(mean, q975, q025)) %>% 
  left_join(select(geodata, ShortName, pop2018)) %>% 
  ungroup() %>% 
  filter(ShortName != "CH") %>% 
  filter(!(ShortName %in% c("BL", "GE", "GR"))) %>%
  mutate(w = 1/sd^2/sum(1/sd^2))

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


# Changepoint models ------------------------------------------------------------------
mobility_change_points <- read_csv("COVID-pomp/results/mob_changepoints.csv")

p_prob_change <- mobility_change_points %>% 
  filter(cp == "cp_1") %>% 
  group_by(ShortName) %>% 
  summarise(p = mean(pvar_after_R0, na.rm = T)) %>%
  arrange(p) %>% 
  ggplot(aes(x = p, y = ShortName)) +
  geom_bar(stat = "identity") +
  labs(x = "Probability changepoint in R0\npreceeds changepoint in mobility", y = "") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10))

# Supplementary figure
ggsave(p_prob_change, filename = "COVID-pomp/results/figs/prob_change_preceed.png", width = 6, height = 6)

# FIGURE 4
p_mob_combined <- arrangeGrob(p_prob_change, pcorr, widths = c(1, 3))
ggsave(p_mob_combined, filename = "COVID-pomp/results/figs/pcomp_mob.png", width = 10, height = 4)

cp_labels <- c(
  cp_1 = "Start decrease",
  cp_2 = "Stabilize"
)

p_change <- mobility_change_points %>% 
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

ggsave(p_change, filename = "COVID-pomp/results/figs/mob_change_dates.png", width = 10, height = 4)

prob_labels <- c(
  prob_1 = as.character(npis$date[1]),
  prob_2 =  as.character(npis$date[2]),
  prob_3 =  as.character(npis$date[3]),
  prob_4 =  as.character(npis$date[4])
)

p_prob <- mobility_change_points %>% 
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

mobility_change_points %>%
  filter(var %in% c("Retail & Recreation", "Transit Stations", "Workplace"), cp =="cp_1") %>% 
  group_by(ShortName) %>% 
  summarise_at(vars(contains("prob")), mean)

# Incidence proportion ---------------------------------------------------------------------

incid_prop <- filterstats %>% 
  filter(var == "tot_I") %>% 
  group_by(ShortName) %>% 
  filter(time >= dateToYears(as.Date("2020-04-24"))) %>% 
  arrange(time) %>%
  slice(1) %>%
  left_join(select(geodata, ShortName, pop2018))

incid_prop$pop2018[is.na(incid_prop$pop2018)] <- sum(geodata$pop2018)

incid_prop_stat <- incid_prop %>% 
  mutate(median_prev = median/pop2018, 
         q025_prev = q025/pop2018,
         q975_prev = q975/pop2018) %>% 
  select(ShortName, contains("prev"))

ch_sf <- sf::st_read("data/ch/shp/ch.shp") %>% 
  left_join(incid_prop %>% 
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
  guides(fill = guide_colorbar(title = "Incidence proportion [%]", title.position = "top"))

ggsave(p_prev, filename = "COVID-pomp/results/figs/incid_prop.png", width = 6, height = 5)


# Results for text -------------------------------------------------------------

cat(resString(r0_left, "R0 left"))
cat(resString(r0_right, "R0 left"))
cat(resString(r0_change, "R0 change"))
cat(resString(r0_time %>% mutate_at(vars(median, contains("q")), function(x) date_start + x), 
              "incid_propalence"))

cat(resString(select(incid_prop_stat, ShortName, contains("prev")) %>% 
                set_colnames(c("ShortName", "median", "q025", "q975")), 
              "incid_propalence"))
