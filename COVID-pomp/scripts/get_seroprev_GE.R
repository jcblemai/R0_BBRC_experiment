library(tidyverse)
library(magrittr)
source("COVID-pomp/scripts/utils.R")

tot_I <- readRDS("COVID-pomp/results/filtered_states_all.rds") %>% 
  filter(var == "tot_I", ShortName == "GE")

pop <- read_csv("data/ch/geodata.csv")
pop_GE <- pop$pop2018[pop$ShortName == "GE"]
week_lim <- as.Date(c("2020-04-20", "2020-04-24")) - 15

tot_I_week <- filter(tot_I, time >= dateToYears(week_lim[1]), time <= dateToYears(week_lim[2])) %>% 
  summarise_at(vars(median, q025, q975), function(x) mean(x/pop_GE)) 

sero_data <- tribble(
  ~date, ~median, ~q025, ~q975,
  "2020-04-06", 2.9, 0.5, 6.5,
  "2020-04-16", 6.7, 3.1, 11,
  "2020-04-21", 9.4,  4.7, 15.2
) %>% mutate(date = as.Date(date),
             time = dateToYears(date))

# sero_data <- tribble(
#   ~date, ~median, ~q025, ~q975,
#   "2020-04-08", 3.5, 1.6, 5.4,
#   "2020-04-15", 5.5, 3.3, 7.7,
# ) %>% mutate(date = as.Date(date),
#              time = dateToYears(date))

tot_I2 <- tot_I %>% 
  mutate(dates = "ori") %>% 
  rbind(tot_I %>% mutate(date = date + 15,
                         dates = "shifted")) %>% 
  mutate_at(vars("mean", "median", contains("q")), function(x) x/pop_GE*100)

filter(tot_I2, date >= min(sero_data$date)) %>% 
  ggplot(aes(x = date, y = median, ymin = q025, ymax = q975)) +
  geom_ribbon(alpha = .1, aes(fill = dates)) +
  geom_ribbon(alpha = .3, aes(ymin = q25, ymax = q75, fill = dates)) +
  geom_point(data = sero_data) +
  geom_errorbar(data = sero_data, width = 0) +
  theme_bw()

cat(glue::glue("Seroprevalence in GE for week: {week_lim[1]}-{week_lim[2]}\n{format(tot_I_week$median, digits = 2)}",
               " (95% QR: {format(tot_I_week$q025, digits = 2)} - {format(tot_I_week$q975, digits = 2)})"))
