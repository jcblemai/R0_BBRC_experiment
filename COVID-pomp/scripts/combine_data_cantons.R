library(tidyverse)

data_ch <- list.files("data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2/", full.names = T) %>% 
  purrr::map(read_csv) %>% 
  bind_rows() %>% 
  filter(abbreviation_canton_and_fl != "FL") %>% 
  select(-abbreviation_canton_and_fl, -time, -source) %>% 
  group_by(date) %>% 
  summarise_all(sum, na.rm = T)


write_csv(data_ch, path = "data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2//COVID19_Fallzahlen_Kanton_CH_total.csv")
