
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(magrittr)



# load data --------------------------------------------------------------------
date_left <- as.Date("2020-02-16")
date_right <- as.Date("2020-03-29")

data <- read_csv("data/google_mobility/retail_recreation.csv", col_names = F) %>% 
  set_colnames(c("time", "change")) %>% 
  mutate(date = date_left + round(time * difftime(date_right, date_left, units = "days")))

load("results/filtered_covidVD_noalpha_noW_newLLh.rda")

filter_stats <- filter_dists %>% 
  group_by(time, parset, var) %>% 
  summarise(mean = mean(value, na.rm = T),
            q025 = quantile(value, 0.025, na.rm = T),
            q975 = quantile(value, 0.975, na.rm = T),
            q25 = quantile(value, 0.25, na.rm = T),
            q75 = quantile(value, 0.75, na.rm = T)) %>%
  ungroup() %>% 
  mutate(date = yearsToDate(time),
         parset = factor(parset))

R0_baseline <- with(filter_stats, mean(mean[var == "Rt" & date <= "2020-03-01"], na.rm = T))

p <- filter_stats %>% 
  filter(var == "Rt") %>% 
  mutate_at(vars(mean, contains("q")), function(x) x = x/R0_baseline - 1) %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = parset), alpha = .2) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = parset), alpha = .2) +
  geom_line(data = data, aes(x = date, y = change)) +
  theme_bw()  +
  labs(y = "% change from baseline", x = "")

ggsave(p, filename = "results/retail_R0_VD.png", width = 6, height = 4)
