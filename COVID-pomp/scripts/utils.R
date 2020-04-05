erlangLL <- function(lambda, x, k) {
  -1 * sum(dgamma(x, shape = k , rate = lambda, log = T))
}

fitErland <- function(x, ks = 1:10, zero_replace = .5) {
  x[x == 0] <- zero_replace
  x <- x[x>0]
  bind_rows(
    lapply(ks, 
           function(k) 
             optimize(erlangLL, x = x, k = k, interval = c(0, 1)) %>% 
             as_tibble() %>% 
             mutate(k = k, 
                    ll = -1 * objective,
                    lambda = minimum,
                    mean = k/lambda,
                    variance = k/(lambda^2)) %>%
             select(k, ll, lambda, mean, variance)
           
    )
  ) %>% 
    arrange(desc(ll)) %>% 
    slice(1)
}



# function to draw random vriables from a Drichillet distribution using independent gammas
rdirichlet <- function (n, alpha) {
  if (!is.null(salpha <- dim(alpha))) {
    x <- matrix(rgamma(prod(salpha) * n, rep(as.vector(as.matrix(alpha)), each = n)), ncol = salpha[2], byrow = F)
    sm <- x %*% rep(1, salpha[2])
  } else {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
  }
  return(x/as.vector(sm))
}
