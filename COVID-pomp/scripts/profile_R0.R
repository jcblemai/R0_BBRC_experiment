# Title: POMP model of COVID19 in CH, cantonal level
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(doSNOW)
library(pomp)
library(magrittr)
library(foreach)
library(itertools)
library(lubridate)
library(parallel)
library(glue)
select <- dplyr::select

option_list = list(
  optparse::make_option(c("-c", "--config"), action="store", default='pomp_config.yaml', type='character', help="path to the config file"),
  optparse::make_option(c("-a", "--asindex"), action="store", default=0, type='numeric', help="whether to use the index of a slurm array"),
  optparse::make_option(c("-b", "--basepath"), action="store", default="COVID-pomp/", type='character', help="base path"),
  optparse::make_option(c("-j", "--jobs"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-o", "--cores"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-r", "--run_level"), action="store", default=1, type='numeric', help="run level for MIF"),
  optparse::make_option(c("-p", "--place"), action="store", default='CH', type='character', help="name of place to be run, a Canton abbrv. in CH"),
  optparse::make_option(c("-l", "--likelihood"), action="store", default='d-deltah', type='character', help="likelihood to be used for filtering")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- yaml::read_yaml(opt$config)

source(glue("{opt$b}scripts/skellam.R"))
source(glue("{opt$b}scripts/mifCooling.R"))
source(glue("{opt$b}scripts/utils.R"))

if (opt$a == 1 & Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
  array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  opt$place <- config$places[array_id]
}

# Setup ------------------------------------------------------------------------
# Read Rscript arguments
canton <- opt$place

# Which likelihood components to use?
# deltah: balances of inputs and outputs from hospitals
# c: cases
# d: total deaths
lik_components <- str_split(opt$likelihood, "-")[[1]]#c("deltah", "c", "d")
lik_log <- str_c(str_c("ll_", lik_components), collapse = "+")
lik <- str_c(str_c("ll_", lik_components), collapse = "*")
downweight <- opt$downweight
# Test for cases in the likelihood
ll_cases <- "c" %in% lik_components
suffix <- buildSuffix(name = config$name, 
                      place = canton,
                      lik_components = lik_components, 
                      sdfrac = config$sdfrac*100, 
                      params_to_fit = config$parameters_to_fit)

# Level of detail on which to run the computations
run_level <- opt$run_level
sir_Np <- c(1e3, 3e3, 6e3)
sir_Nmif <- c(2, 20, 150)
sir_Ninit_param <- c(2, opt$jobs, opt$jobs)
sir_NpLL <- c(1e3, 1e4, 1e4)
sir_Nreps_global <- c(2, 5, 20)

ll_filename <- glue("{opt$b}results/loglik_exploration_{suffix}.csv")
mif_filename <- glue("{opt$b}results/profiling_mif_sir_sde_{suffix}.rda")
ll_prof_filename <- glue("{opt$b}results/profiling_loglik_{suffix}.csv")

# Initial parameters -----------------------------------------------------------

# Load MLE
liks <- read_csv(ll_filename)

best_params <- liks %>%
  arrange(desc(loglik)) %>% 
  select(-contains("log")) %>% 
  slice(1) %>% as.vector()

load(glue("{opt$b}interm/pomp_{suffix}.rda"))
params <- covid@params
# values of the random walks standard deviations
# values of the random walks standard deviations

# values of the random walks standard deviations
rw.sd_rp <- 0.02  # for the regular (process and measurement model) parameters
rw.sd_rp_prior <- 0.005  # for the more uncertain random walk on std_X
rw.sd_ivp <- 0.2  # for the initial value parameters
rw.sd_ivp_prior <- 0.05  # for the initial value parameters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_rp_prior, rw.sd_ivp, rw.sd_ivp_prior), c("regular", "prior", "ivp", "ivpprior"))

# lower bound for positive parameter values
min_param_val <- 1e-5 
# define the bounds for the parameters to estimate
parameter_bounds <- tribble(
  ~param, ~lower, ~upper,
  # Process noise
  "std_X", 1, 4, #in log-scale
  # Initial conditions
  "I_0", 2e-5, 4e-5
)

if(ll_cases) {
  parameter_bounds <- rbind(parameter_bounds, 
                            tribble(~param, ~lower, ~upper,
                                    # Measurement model
                                    "k", .1, 10,
                                    "epsilon", 0.2, 0.5))}

if (!is.null(config$parameters_to_fit)) {
  # additional params to fit
  other_bounds <- mapply(x = names(config$parameters_to_fit),
                         y = config$parameters_to_fit, 
                         function(x, y) tibble(param = x, lower = y$lower, upper = y$upper),
                         SIMPLIFY = F) %>% 
    bind_rows()
  parameter_bounds <- rbind(parameter_bounds, other_bounds)
}

# if (!ll_cases) {
#   parameter_bounds <- parameter_bounds %>% filter(!(param %in% c("k", "epsilon")))
# }

# convert to matrix for ease
parameter_bounds <- set_rownames(as.matrix(parameter_bounds[, -1]), parameter_bounds[["param"]])

# Sequence of R0s over which to profile
R0_seq <- seq(1.5, 4.5, length.out = opt$jobs)
n_init <- 6   # number of initializations per parameter value

# create random vectors of initial parameters given the bounds
init_params <- sobolDesign(lower = parameter_bounds[, "lower"],
                           upper = parameter_bounds[, "upper"], 
                           nseq = n_init) 

init_params <- bind_rows(lapply(1:length(R0_seq), function(x) init_params))

# convert certain parameters back to log-scale 
# param_logbound_names <- c("std_X")
# init_params[, param_logbound_names] <- exp(init_params[, param_logbound_names])

param_fixed_names <- names(params)[unlist(lapply(names(params), function(x) !(x %in% rownames(parameter_bounds))))]

# bind with the fixed valued parameters
init_params <- cbind(init_params, 
                     matrix(rep(params[param_fixed_names],
                                each = sir_Ninit_param[run_level]),
                            nrow = sir_Ninit_param[run_level]) %>% 
                       set_colnames(param_fixed_names))

# Set R0 to profiling values
init_params$R0_0 <- rep(R0_seq, each = n_init)

if (!is.null(config$parameters_to_fit)) {
  other_rw <- lapply(names(config$parameters_to_fit), 
                     function(x) glue(", {x} = {rw.sd_param['regular']}")) %>% 
    unlist() %>% 
    str_c(sep = " ")
} else {
  other_rw <- ""
}
tvary <- dateToYears(as.Date(config$tvary))
rw_text <- glue("rw.sd( std_X  =ifelse(time>={tvary}, {rw.sd_param['regular']}, {rw.sd_param['regular']/10}) 
                , k  = {ifelse(ll_cases, rw.sd_param['regular'], 0)}
                , epsilon   = {ifelse(ll_cases, rw.sd_param['regular'], 0)}
                 , I_0  = ivp({rw.sd_param['ivp']})
                {other_rw})")

job_rw.sd <- eval(parse(text = rw_text))

# MIF --------------------------------------------------------------------------
# parallel computations
cl <- makeCluster(opt$cores)
registerDoSNOW(cl)

# MIF it! 
t1 <- system.time({
  mf <- foreach(parstart = iter(init_params, by = "row"),
                .inorder = F, 
                .packages = "pomp",
                .errorhandling = "stop") %dopar% {
                  mif2(covid,
                       params = unlist(parstart),
                       Np = sir_Np[run_level],
                       Nmif = sir_Nmif[run_level],
                       cooling.type = "geometric",
                       cooling.fraction.50 = findAlpha(sir_Nmif[run_level], ncol(covid@data), 0.05),
                       rw.sd = job_rw.sd,
                       verbose = F)
                }})

save(t1, mf, file = mif_filename)

cat("----- Done MIF, took", round(t1["elapsed"]/60), "mins \n")

# Log-lik ------------------------------------------------------------------

t2 <- system.time({
  prof_liks <- foreach(mfit = mf,
                       .inorder = T, 
                       .packages = "pomp",
                       .combine = rbind,
                       .errorhandling = "stop"
  ) %dopar% {
    # compute log-likelihood estimate by repeating filtering with the given param vect
    ll <-  logmeanexp(
      replicate(sir_Nreps_global[run_level],
                logLik(
                  pfilter(covid,
                          params = pomp::coef(mfit),
                          Np = sir_NpLL[run_level])
                )
      ), se = TRUE)
    # save to dataframe
    data.frame(loglik = ll[1], loglik_se = ll[2], t(coef(mfit)))
  }
})

write_csv(prof_liks, path = ll_prof_filename, append = file.exists(ll_prof_filename))

cat("----- Done LL, took", round(t2["elapsed"]/60), "mins \n")
stopCluster(cl)


# Compute profile --------------------------------------------------------------
source(glue("{opt$b}scripts/mcap.R"))
prof_liks <- read_csv(ll_prof_filename) %>% 
  arrange(R0_0) %>% 
  group_by(R0_0) %>% 
  arrange(desc(loglik)) %>% 
  filter(loglik > -450) %>%
  slice(1:5)%>%
  filter(R0_0 > 2 , R0_0 < 4.5)

R0_mcap <- mcap(prof_liks$loglik, prof_liks$R0_0, lambda = 1)

p <- ggplot(R0_mcap$fit, aes(x = parameter)) +
  geom_line(aes(y = quadratic), color = "blue") +
  geom_line(aes(y = smoothed), color = "red") +
  geom_point(data = prof_liks, aes(x = R0_0, y = loglik)) +
  geom_vline(aes(xintercept = R0_mcap$ci[1]), lty = 2, size = .3) +
  geom_vline(aes(xintercept = R0_mcap$ci[2]), lty = 2, size = .3) +
  geom_hline(aes(yintercept = R0_mcap$llmax - R0_mcap$delta), lty = 2, size = .3) +
  theme_minimal()
p
ggsave(p, filename = glue("{opt$b}results/figs/profiling_R0_{suffix}.png"), width = 5, height = 4)
paste0(format(R0_mcap$mle, digits = 3), " (",
       paste(format(R0_mcap$ci, digits = 3), collapse = "-"),
       ")")

ggplot(prof_liks, aes(x = R0_0)) +
  geom_point(aes(y = loglik)) + 
  # geom_errorbar(aes(ymin = loglik - loglik_se * 1.96, ymax = loglik + loglik_se * 1.96), width = 0) +
  geom_smooth(aes(y = loglik))

