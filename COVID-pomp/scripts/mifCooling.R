mifCooling <- function(n_mif, n_data, alpha) {
  cooling <- matrix(rep(1, n_mif * n_data), ncol = n_data)
  for (i in 1:n_mif) {
    for(j in 1:n_data) {
      cooling[i, j] <- alpha^((j - 1 + (i-1) * n_data)/(50 * n_data))
    }
  }
  return(cooling)
}

alphaDiff <- function(alpha, alpha_end, n_mif, n_data) {
  (alpha_end - mifCooling(n_mif, n_data, alpha)[n_mif, 1])^2
}

findAlpha <- function(n_mif, n_data, alpha_end = .1) {
  optimise(alphaDiff, interval = c(0, 1), alpha_end = alpha_end, n_mif = n_mif, n_data = n_data)$minimum
}

