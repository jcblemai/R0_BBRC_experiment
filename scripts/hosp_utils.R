# Crude estimates (biased)
computeCFRBiased <- function(cumul_cases, cumul_deaths, n_lag) {
  cumul_deaths/lag(cumul_cases, n_lag)
}

convolveManual <- function(x, y) {
  res <- rep(0, length(x))
  for (i in 1:length(x)) {
    for (j in 1:i) {
      res[i] <- res[i] + x[j] * y[i - j + 1]
    }
  }
  res
}

# Adjusted estimates following Nishiura et al. 2009
computeU <- function(C_t, f_t2d) {
  n_elem_C <- length(C_t)
  convolveManual(C_t, f_t2d)/C_t
}

profileLik <- function(ps, uC, D) {
  liks <- lapply(ps, function(p) {
    dbinom(D, ceiling(uC), p, log = T)
  }) %>% 
    unlist()
  
  rliks <- 2 * (max(liks) - liks)
  chilim <- qchisq(0.95, 1)
  lb <- which(rliks < chilim)[1]
  mle_i <- which.max(liks)[1]
  ub <- which(rliks[mle_i:length(ps)] > chilim)[1] + mle_i - 1
  q025 <- ps[lb]
  q975 <- ps[ub]
  data.frame(mle = ps[mle_i], q025 = q025, q975 = q975)
}

computeCCFR <- function(t2d, C_t, D_t, dates, distr) {
  
  f_t2d_params <- fitdistrplus::fitdist(t2d, distr = distr)$estimate
  if (distr == "lnorm") {
    f_t2d <- dlnorm(seq(1, 40), f_t2d_params[1], f_t2d_params[2])
  } else {
    f_t2d <- dgamma(seq(1, 40), f_t2d_params[1], f_t2d_params[2])
  }
  
  # As defined in Nishiura et al. 2009
  u_t <- computeU(C_t, f_t2d)
  uC_t <- u_t * C_t
  ps <- seq(0, 1, by = 1e-2)
  ccfrs <- do.call(
    rbind, 
    lapply(1:length(u_t), function(i) {profileLik(ps, uC_t[i], D_t[i])})) %>% 
    mutate(date = dates)
  return(ccfrs)
}


fitDist <- function(df, col, distr, zero_replace = 0.5) {
  data <- df[[col]]
  data[data == 0] <- zero_replace
  fitdist(data[data>0], dist = distr)
}

generateDenstiy <- function(fit) {
  obs_range <- seq(0.1, max(fit$data) * 1.1, by = .1)
  
  if (fit$distname == "gamma") {
    dens <- dgamma(obs_range, fit$estimate[1], fit$estimate[2]) 
  } else if (fit$distname == "lnorm") {
    dens <- dlnorm(obs_range, fit$estimate[1], fit$estimate[2]) 
  }  else {
    dens <- dexp(obs_range, fit$estimate[1]) 
  }
  return(data.frame(x = obs_range, y = dens))
}

