
# Define C code for the skellam distribution
dskellam.C <- "double dskellam(int k, double mu1, double mu2, int give_log) {
  double ll, lik;
  
  if (mu1 == 0) {
    lik = dpois(-k, mu2, give_log);
  } else if (mu2 == 0) {
    lik = dpois(k, mu1, give_log);
  } else {
  double s12 = sqrt(mu1 * mu2);
  ll = -mu1 - mu2 + 2 * s12 + 0.5 * k * (log(mu1) - log(mu2)) + log(bessel_i(2 * s12, abs(k), 2));
  if (give_log == 0) {
    lik = exp(ll);
  } else {
    lik = ll;
  }
  }
  
  return lik;
};"

rskellam.C <- "int rskellam(double mu1, double mu2) {
  int k;
  k = rpois(mu1) - rpois(mu2);
  return k;
};"
