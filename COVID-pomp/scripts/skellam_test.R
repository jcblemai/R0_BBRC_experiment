
library(pomp)
library(skellam)

# Define C code for the skellam distribution
dskellam.C <- "double dskellam(int k, double mu1, double mu2, int give_log) {
  double ll, lik;
  double s12 = sqrt(mu1 * mu2);
  ll = -mu1 - mu2 + 2 * s12 + 0.5 * k * (log(mu1) - log(mu2)) + log(bessel_i(2 * s12, abs(k), 2));

  if (give_log == 0) {
    lik = exp(ll);
  } else {
    lik = ll;
  }
  return lik;
};"

rskellam.C <- "int rskellam(double mu1, double mu2) {
  int k;
  k = rpois(mu1) - rpois(mu2);
  return k;
};"


mu1 <- 10
mu2 <- 1
mus = c(mu1 = mu1, mu2 = mu2)
nobs <- 1000
data <- data.frame(time = seq(1:nobs),
                   obs = rskellam(nobs, mu1, mu2))
check_fit <- skellam.mle(data$obs)

skellam_test <- pomp(data,
                     times = "time",
                     t0 = 1,
                     rinit = Csnippet("x1 = mu1; x2 = mu2;"),
                     rprocess = discrete_time(Csnippet("x1 = mu1; x2 = mu2;")),
                     dmeasure = Csnippet("lik = dskellam(obs, x1, x2, give_log);"),
                     rmeasure = Csnippet("obs = rskellam(mu1, mu2);"),
                     partrans =  parameter_trans(log = c("mu1", "mu2")),
                     params = c("mu1" = 10, "mu2" = 1),
                     paramnames = c("mu1", "mu2"),
                     statenames = c("x1", "x2"),
                     obsnames = c("obs"),
                     globals = paste(dskellam.C, rskellam.C, sep = "\n")
                     )

library(ggplot2)
simulate(skellam_test, 10, format = "data.frame", include.data = F) %>% 
  ggplot(aes(x = obs)) + geom_histogram(position = "dodge") + geom_histogram(data = data, aes(x = obs), fill = "red", position = "dodge")

# fit <- mif2(skellam_test,
#      params = c("mu1" = 1, "mu2" = 1),
#      Np = 200,
#      Nmif = 50,
#      cooling.type = "geometric",
#      cooling.fraction.50 = 0.1,
#      rw.sd = rw.sd(mu1 = 0.02, mu2 = 0.02),
#      verbose = F)

coef(fit)
logLik(pfilter(skellam_test, params = mus, 1000))
# check_fit
sum(dskellam(data$obs, mu1, mu2, log = T))
