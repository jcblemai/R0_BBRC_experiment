library(tidyverse)

data_ch <- list.files("data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2/", full.names = T) %>% 
  purrr::map(read_csv) %>% 
  bind_rows() %>% 
  filter(abbreviation_canton_and_fl != "FL") %>% 
  select(-abbreviation_canton_and_fl, -time, -source) %>% 
  group_by(date) %>% 
  summarise_all(sum, na.rm = T) %>% 
  arrange(date) %>% 
  mutate(deaths = c(ncumul_deceased[1], diff(ncumul_deceased))) %>% 
  filter(date < "2020-04-10")

write_csv(data_ch, path = "data/ch/cases/covid_19/fallzahlen_kanton_total_csv_v2/COVID19_Fallzahlen_Kanton_CH_total.csv")
write_csv(data_ch, path = "~/Desktop/openzh_CH_aggregate.csv")

data_alth <- read_csv("data/ch/deaths_CH_althaus.csv")

data_oz2 <- read_csv("data/ch/cases/covid_19/COVID19_Fallzahlen_CH_total_v2.csv")

data_oz2_deaths <- select(data_oz2, date, abbreviation_canton_and_fl, ncumul_deceased) %>% 
  magrittr::set_colnames(c("date","canton","cum_deaths"))%>% 
  mutate(source = "oz")

data_oz2_deaths <- data_oz2_deaths %>% 
  rbind(data_oz2_deaths %>% group_by(date, source) %>% 
          summarise(cum_deaths = sum(cum_deaths, na.rm =T)) %>%
          mutate(canton = "CH") %>% 
          ungroup %>% 
          select(one_of(colnames(data_oz2_deaths))))
        
data_db <- read_csv("data/ch/covid19_fatalities_switzerland_openzh.csv") %>% 
  gather(canton, deaths, -Date) %>% 
  rename(date = Date, cum_deaths = deaths) %>% 
  mutate(source = "db")

data_comp <- inner_join(data_oz2_deaths, data_db, by = c("date","canton"), suffix = c(".oz", ".db"))
data_comp <- rbind(data_oz2_deaths, data_db)

# Compare cantons
ggplot(data_comp, aes(x = cum_deaths.oz, y = cum_deaths.db)) +
  geom_point(size = .5) +
  facet_wrap(~canton, scales = "free")

# Compare all

# ------------------------------------------------------------------------------


data_oz2 <- data_oz2 %>% 
  group_by(date) %>% 
  summarise(cum_deaths = sum(ncumul_deceased, na.rm =  T)) %>% 
  mutate(deaths = cum_deaths - lag(cum_deaths, 1)) %>% 
  filter(date < "2020-04-10") %>% 
  mutate(source = "oz")



deathsOZ2 <- data2 %>% 
  filter(canton == "CH") %>% 
  group_by(date) %>% 
  summarise(cum_deaths = sum(cum_deaths, na.rm = T)) %>% 
  mutate(source = "oz2",
         deaths = c(NA, diff(cum_deaths)))

ggplot(rbind(data_oz2, deathsOZ2), aes(x= date, y = cum_deaths, color = source)) +
  geom_point()
  

ggplot(data_ch, aes(x = date, y = deaths)) +
  geom_point() +
  geom_point(data = data_alth, col = "red") +
  geom_point(data = data_oz2, col = "blue")
  
