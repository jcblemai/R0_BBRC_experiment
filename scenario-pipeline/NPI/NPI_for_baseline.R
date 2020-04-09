library(dplyr)
library(readr)
library(stringr)
if (F) {         # or whatever values you use to test.
  ti_str <- '2020-01-31'
  tf_str <- '2020-08-31'
  foldername <- 'data/ch/'
  setupname <- 'TestIsolate'
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
r0_reduction <- read_csv(paste0(foldername, "r0_reduction.csv"))

# Compute the mean and bounds of the initial R0 from which to sample the simulations from
left_R0_range <- with(filter(r0_reduction, var == "r0_left"), c(mean(q025), mean(q975)))
left_R0_mean <- with(filter(r0_reduction, var == "r0_left"), c(mean(mean)))
left_R0_sd <- getSd(left_R0_mean, left_R0_range[1], left_R0_range[2])
# Range of reductions: R0_after = reduction * R0_before
right_R0_range <- with(filter(r0_reduction, var == "r0_right"), c(mean(q025), mean(q975)))
right_R0_mean <- with(filter(r0_reduction, var == "r0_right"), mean(mean))
right_R0_sd <- getSd(right_R0_mean, right_R0_range[1], right_R0_range[2])

# Assume reduction on the 20th of March
start_descent <- as.Date("2020-03-15")
end_descent <- as.Date("2020-03-22")
start_ascent <- as.Date("2020-05-01")
end_ascent <- as.Date("2020-05-03")

places <- read_csv(paste0(foldername,'geodata.csv'))
dates <- seq.Date(as.Date(ti_str), as.Date(tf_str), 1)
baseline_dates <- which(dates <= start_descent)
descent_dates <- which(dates > start_descent & dates <= end_descent)
low_dates <- which(dates > end_descent & dates <= start_ascent)
ascent_dates <- which(dates > start_ascent & dates <= end_ascent)
release_dates <- which(dates > end_ascent)

NPI <- as.data.frame(matrix(0, dim(places)[1],length(dates)))
colnames(NPI) <- as.Date(dates)
rownames(NPI) <- places$ShortName

# Set Baseline cantons either with inferred R0s when available or with mean
for (cnt in places$ShortName) {
  if (cnt %in% r0_reduction$ShortName) {
    # Basline
    ind_l <- r0_reduction$var == "r0_left" & r0_reduction$ShortName == cnt
    mul <- r0_reduction$mean[ind_l]
    q025l <- r0_reduction$q025[ind_l]
    q975l <- r0_reduction$q975[ind_l]
    sdl <- getSd(mul, q025l,  q975l) 
    # R0 value at baseline
    rval_l <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = mul, sd = sdl)
    NPI[cnt, baseline_dates] <- rval_l
    
    # Low dates (during measures)
    ind_r <- r0_reduction$var == "r0_right" & r0_reduction$ShortName == cnt
    mur <- r0_reduction$mean[ind_r]
    q025r <- r0_reduction$q025[ind_r]
    q975r <- r0_reduction$q975[ind_r]
    sdr <- getSd(mur, q025r,  q975r) 
    # R0 value during measures
    rval_r <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = mur, sd = sdr)
    NPI[cnt, low_dates] <- rval_r
    
    # HERER IS THE DIFFERENCE BETWEEN SCENARIOS
    # release_rval is the value after releasing measures
    if (grepl("Current", setupname)) {
      # In Current keep R0 value same to during measures
      release_rval <- rval_r
    } else if (grepl("Stopped", setupname)) {
      # In Stopped go back to value similar to before measures
      release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = mul, sd = sdl)
    } else if (grepl("TestIsolate", setupname)) {
      # In contact tracing ?
      release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = mul * rfrac_testisolate, sd = sdl)
    }
  } else {
    rval_l <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = left_R0_mean, sd = left_R0_sd)
    NPI[cnt, baseline_dates] <- rval_l
    rval_r <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = right_R0_mean, sd = right_R0_sd)
    NPI[cnt, low_dates] <- rval_r
    
    # HERER IS THE DIFFERENCE BETWEEN SCENARIOS
    #    # HERER IS THE DIFFERENCE BETWEEN SCENARIOS
  
    if (grepl("Current", setupname)) {
      release_rval <- rval_r
    } else if (grepl("Stopped", setupname)) {
      release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = left_R0_mean, sd = left_R0_sd)
    } else if (grepl("TestIsolate", setupname)) {
      release_rval <- truncnorm::rtruncnorm(n = 1, a = 0, b = 3.5, mean = left_R0_mean * rfrac_testisolate, sd = left_R0_sd)
    }
  }
  # Release value
  NPI[cnt, release_dates] <- release_rval
  # Transitions
  NPI[cnt, descent_dates] <- rval_l + seq_along(descent_dates) * (rval_r-rval_l)/length(descent_dates)
  NPI[cnt, ascent_dates] <- rval_r + seq_along(ascent_dates) * (release_rval-rval_r)/length(ascent_dates)
}

# Plot to check
# as_tibble(t(NPI)) %>%
#   mutate(date = dates) %>%
#   gather(canton, value, -date) %>%
#   ggplot(aes(x = date, y = value, col = canton)) +
#   geom_line()
