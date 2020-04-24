library(dplyr)
library(readr)
library(stringr)
if (F) {         # or whatever values you use to test.
  ti_str <- '2020-01-31'
  tf_str <- '2020-08-31'
  foldername <- '../data/ch/'
  setupname <- 'Current'
}

# Function
getSd <- function(mu, q025, q975) {
  sd1 <- (mu - q025)/2
  sd2 <- (q975 - mu)/2
  return(mean(c(sd1, sd2)))
}

# Percent of baseline for TestIsolate
rfrac_testisolate <- 0.7
print(getwd())
# Observed reductions
r0_reduction <- read_csv(paste0(foldername, "R0_reduction.csv"))

# Range of reductions: R0_after = reduction * R0_before
change_R0_range <- with(filter(r0_reduction, var == "r0change"), c(mean(q025), mean(q975)))
change_R0_mean <- with(filter(r0_reduction, var == "r0change"), mean(mean))
change_R0_sd <- getSd(change_R0_mean, change_R0_range[1], change_R0_range[2])

# Assume reduction on the 20th of March
#start_descent <- as.Date("2020-03-15")
start_descent <- as.Date("2020-03-14")
end_descent <- as.Date("2020-03-28")
start_ascent <- as.Date("2020-05-01")
end_ascent <- as.Date("2020-05-03")

places <- read_csv(paste0(foldername,'geodata.csv'))
dates <- seq.Date(as.Date(ti_str), as.Date(tf_str), 1)
baseline_dates <- which(dates <= start_descent)
descent_dates <- which(dates > start_descent & dates <= end_descent)
low_dates <- which(dates > end_descent & dates <= start_ascent)
ascent_dates <- which(dates > start_ascent & dates <= end_ascent)
release_dates <- which(dates > end_ascent)

# Random draw of baseline R0
baseline_R0s <- runif(dim(places)[1], 2, 3)
names(baseline_R0s) <-  places$ShortName
NPI <- as.data.frame(matrix(0, dim(places)[1],length(dates)))
colnames(NPI) <- as.Date(dates)
rownames(NPI) <- places$ShortName

# Set Baseline cantons either with inferred R0s when available or with mean
for (cnt in places$ShortName) {
  baseline_R0 <- baseline_R0s[cnt]
  NPI[cnt, baseline_dates] <- baseline_R0
  rval_r <- baseline_R0 * truncnorm::rtruncnorm(n = 1, a = 0, b = 2, mean = change_R0_mean, sd = change_R0_sd)
  rval_r <- baseline_R0 / baseline_R0 * runif(1, .4, .9) #hack !
  NPI[cnt, low_dates] <- rval_r
  
  # HERER IS THE DIFFERENCE BETWEEN SCENARIOS
  if (grepl("Current", setupname)) {
    release_rval <- rval_r
  } else if (grepl("Stopped", setupname)) {
    release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = baseline_R0, sd = 0.07)
  }  else if (grepl("R0-0dot9", setupname)) {
      release_rval <- truncnorm::rtruncnorm(n = 1, a = 0.8, b = 1, mean = 0.9, sd = 0.07)
  }  else if (grepl("R0-1dot2", setupname)) {
    release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = 1.2, sd = 0.07)
  } else if (grepl("R0-1dot5", setupname)) {
    release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = 1.5, sd = 0.07)
  } else if (grepl("TestIsolate", setupname)) {
    release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = baseline_R0 * rfrac_testisolate, sd = 0.07)
  }
  # Release value
  NPI[cnt, release_dates] <- release_rval
  # Transitions
  NPI[cnt, descent_dates] <- baseline_R0 + seq_along(descent_dates) * (rval_r-baseline_R0)/length(descent_dates)
  NPI[cnt, ascent_dates] <- rval_r + seq_along(ascent_dates) * (release_rval-rval_r)/length(ascent_dates)
}

# Plot to check
# as_tibble(t(NPI)) %>%
#   mutate(date = dates) %>%
#   gather(canton, value, -date) %>%
#   ggplot(aes(x = date, y = value, col = canton)) +
#   geom_line()
