library(tidyverse)
source("COVID-pomp/scripts/utils.R")
# BAG data
data <- readxl::read_excel("data/ch/2020-04-17_11-50-22_swiss_reporting_ssph[4].xlsx")

ts2h <- data %>% 
  filter(hospitalisation == 1) %>% 
  filter(manifestation_dt > "2020-02-01") %>% 
  # filter(manifestation_dt > "2020-03-15") %>% 
  mutate(time_symp_to_hosp = difftime(hospdatin, manifestation_dt, units = "days")) %>% 
  filter(time_symp_to_hosp >= 0) %>% 
  .$time_symp_to_hosp %>% 
  as.numeric()

hist(ts2h, 30)

fitErland(ts2h)
fitdistrplus::fitdist(ts2h[ts2h>0] , "lnorm")

