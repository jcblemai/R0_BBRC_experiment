library(tidyverse)
library(timetk)
# R0 Data
filterstats <- readRDS("COVID-pomp/results/filtered_states_all.rds")
r_knot <- filter(filterstats, var  == "Rt") %>% 
  select(date, time, ShortName, median) %>% 
  rename(R0 = median) %>% 
  group_by(ShortName) %>% 
  mutate(date = seq.Date(min(date), by = "1 days", length.out = length(date))) %>% 
  arrange(date) %>% 
  mutate(baseline = mean(R0[!is.na(R0)][1:5]),
         relative = (R0/baseline - 1)*100)
rm(filterstats)

# Mobility data 
gdata <- read_csv("data/ch/Global_Mobility_Report.csv") %>% 
  filter(country_region == "Switzerland") %>% 
  replace_na(list(sub_region_1 = "CH"))


cantons_dict <- tribble(
  ~name, ~ShortName,
  "CH", "CH",
  "Fribourg", "FR",
  # "Aargau", "AG",
  "Geneva", "GE",
  "Canton of Bern", "BE",
  "Vaud", "VD",
  "Valais", "VS",
  "Neuchâtel", "NE",
  "Zurich", "ZH", 
  "Jura" , "JU",
  "Ticino", "TI",
  "Basel-Landschaft", "BL",
  # "Lucerne", "LU",
  "Basel City", "BS",
  "Uri", "UR",
  "Grisons", "GR")

fillGaps <- function(dates, values) {
  
  if(sum(is.na(values)) > 0.5 *length(values)) {
    return(rep(NA, length(values)))
  } else {
    zoo_obj <- zoo::zoo(values, dates)
    zoo_obj <- zoo::na.approx(zoo_obj, rule=2)
    return(as.vector(zoo_obj))
  }
}
gdata_CH <- gdata %>% 
  select(-contains("not_enough")) %>% 
  inner_join(cantons_dict, by = c("sub_region_1" = "name")) %>% 
  select(-contains("region")) %>% 
  gather(var, value, -ShortName, -date) %>% 
  mutate(var = str_replace_all(var, "_percent_change_from_baseline", ""))

gdata_CH <- gdata_CH %>%
  group_by(ShortName, var) %>% 
  arrange(date) %>% 
  mutate(filled = fillGaps(date, value)) %>% 
  mutate(rollmean = roll_apply_vec(filled,  ~ mean(.), .period = 7))

write_csv(gdata_CH, path = "data/google_mobility_CH.csv")

getCorr <- function(x, y) {
  crosscorr <- ccf(x, y, lag.max = 7, plot = F, na.action = na.omit)
  crosscov <- ccf(x, y, lag.max = 7, plot = F, na.action = na.omit, type = "covariance")
  lim <-  qnorm((1 + 0.95)/2)/sqrt(crosscorr$n.used)
  lastLag <- which.max(abs(crosscorr$acf) >= lim) %>% .[length(.)]
  ilag <-  which.max(abs(crosscorr$acf))[1]
  sdx <- sd(x)
  sdy <- sd(y)
  coef <- coef(lm(r ~ value, data = data.frame(r = x, value = y)))[2]
  max.cor <- crosscorr$acf[ilag]
  max.coef <- max.cor * sdx/sdy
  max.cov <- crosscov$acf[ilag]
  pval <- 2 * (1 - pnorm(abs(max.cor), mean = 0, sd = 1/sqrt(crosscorr$n.used)))
  
  data.frame(corr = max.cor,
             corr.low = psychometric::CIr(max.cor, crosscorr$n.used)[1],
             corr.high = psychometric::CIr(max.cor, crosscorr$n.used)[2],
             cov = max.cov,
             max.coef = max.coef,
             lag = crosscorr$lag[ilag],
             pval = pval) %>% 
    mutate(max.coef.low = corr.low * sdx/sdy,
           max.coef.high = corr.high * sdx/sdy)
}

joined <- left_join(gdata_CH %>% mutate(julian = lubridate::yday(date)), 
                    r_knot %>% ungroup %>%  mutate(julian = lubridate::yday(date)), 
                    by = c("ShortName", "julian")) %>% 
  group_by(ShortName, var) %>% 
  mutate(relative = case_when(is.na(relative) ~ 0, T ~ relative),
         rollmean = case_when(is.na(rollmean) ~ 0, T ~ rollmean)) %>% 
  select(relative, rollmean, date.x) %>% 
  drop_na()

crosscorrs <- joined %>% 
  filter(date.x >= "2020-03-05", date.x <= "2020-04-05") %>% 
  filter(var != "residential") %>% 
  group_by(ShortName, var) %>% 
  group_map(~getCorr(.x$relative, .x$rollmean) %>% 
              mutate(ShortName = .y[[1]], var = .y[[2]])) %>% 
  bind_rows() %>% 
  mutate(lag.type = case_when(lag == 0 ~ "simultaneous",
                              lag < 0 ~ "lagging",
                              lag > 0 ~ "preceding"))

write_csv(crosscorrs, path= "COVID-pomp/results/mobility_cross_correlations.csv")

ggplot(crosscorrs, aes(x = lag)) + 
  geom_histogram() +
  facet_wrap(~var)


ggplot(crosscorrs, aes(x = lag, y = ShortName, fill = corr, color = corr)) + 
  geom_bar(stat= "identity") + 
  geom_point() +
  facet_wrap(~var) +
  scale_fill_gradient2() +
  scale_color_gradient2()

ggplot(crosscorrs, aes(x = ShortName, y = corr, color = lag.type)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = corr.low, ymax = corr.high), width = 0) +
  geom_hline(aes(yintercept = 0), lty = 3) +
  coord_flip() +
  facet_wrap(~var) +
  theme_bw() +
  ggthemes::scale_color_few()

p <- ggplot(gdata_CH, aes(x = date)) +
  geom_line(aes(y = filled, color = var), alpha = .3) +
  geom_line(aes(y = rollmean, color = var)) +
  geom_line(data = r_knot, aes(y = relative)) +
  facet_wrap(~ShortName) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  theme_bw() +
  coord_cartesian(ylim = c(-100, 150)) +
  ggthemes::scale_color_few()


ggplot(crosscorrs, aes(x = max.coef, y = ShortName)) + 
  geom_point() +
  geom_errorbarh(aes(xmin = max.coef.low, xmax = max.coef.high), height = 0) +
  geom_vline(aes(xintercept = 1), lty =3, col ="red") +
  geom_vline(aes(xintercept = 0), lty =2, col ="red") +
  facet_wrap(~var) +
  theme_bw()