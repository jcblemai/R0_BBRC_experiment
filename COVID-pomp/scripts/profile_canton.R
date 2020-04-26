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
library(optparse)

select <- dplyr::select

option_list = list(
  make_option(c("-c", "--config"), action="store", default='pomp_config.yaml', type='character', help="path to the config file"),
  make_option(c("-a", "--asindex"), action="store", default=0, type='numeric', help="whether to use the index of a slurm array"),
  make_option(c("-b", "--basepath"), action="store", default="COVID-pomp/", type='character', help="base path"),
  make_option(c("-j", "--jobs"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  make_option(c("-o", "--cores"), action="store", default=detectCores(), type='numeric', help="number of cores used"),
  make_option(c("-r", "--run_level"), action="store", default=1, type='numeric', help="run level for MIF"),
  make_option(c("-p", "--place"), action="store", default='CH', type='character', help="name of place to be run, a Canton abbrv. in CH"),
  make_option(c("-l", "--likelihood"), action="store", default='d-deltah', type='character', help="likelihood to be used for filtering"),
  make_option(c("-x", "--to_profile"), action="store", default='R0_0', type='character', help="Parameter over which to profile"),
  make_option(c("-n", "--n_prof"), action="store", default=4, type='numeric', help="Number of initializations per parameter"),
  make_option(c("-s", "--suffix"), default = "", type = "character", help = "custom suffix to add")
)
opt <- parse_args(OptionParser(option_list=option_list))
config <- yaml::read_yaml(opt$config)

source(glue("{opt$b}scripts/skellam.R"))
source(glue("{opt$b}scripts/mifCooling.R"))
source(glue("{opt$b}scripts/utils.R"))

if (opt$a == 1 & Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
  array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  opt$place <- config$places[array_id]
}

# Level of detail on which to run the computations
run_level <- opt$run_level
sir_Np <- c(1e3, 3e3, 6e3)
sir_Nmif <- c(2, 20, 150)
sir_Ninit_param <- c(2, opt$jobs, opt$jobs)
sir_NpLL <- c(1e3, 1e4, 1e4)
sir_Nreps_global <- c(2, 5, 20)

# Setup ------------------------------------------------------------------------
place <- opt$place

# Parse likelihood components to use
parsed_lik <- parseLikelihood(opt$likelihood)

# Check cases in the likelihood
use_case_incid <- "c" %in% parsed_lik$components

suffix <- buildSuffix(
  name = config$name,
  place = place,
  lik_components = parsed_lik$components,
  params_to_fit = config$parameters_to_fit,
  other =  c(config$sdfrac * 100, opt$suffix)
)

suffix_prof <- str_c(suffix, str_c("_prof-", opt$to_profile, collapse = "-"))
mif_prof_filename <- glue("{opt$b}results/profiling_mif_{suffix_prof}.csv")
ll_prof_filename <- glue("{opt$b}results/profiling_loglik_{suffix_prof}.csv")

# Initial parameters -----------------------------------------------------------
source(glue("{opt$b}scripts/{config$model}"))
load(glue("{opt$b}interm/pomp_{suffix}.rda"))
params <- covid@params

# values of the random walks standard deviations
rw.sd_rp <- 0.02 # for the regular (process and measurement model) parameters
rw.sd_ivp <- 0.2 # for the initial value parameters
rw.sd_param <- set_names(c(rw.sd_rp, rw.sd_ivp), c("regular", "ivp"))

pars_to_profile <- lapply(config$profile_bounds, function(x) seq(x$lower, x$upper, length.out = opt$jobs))
pars_to_profile <- pars_to_profile[opt$to_profile]
init_params <- getInitParams(config = config, 
                             params = params,
                             param_fixed_names = param_fixed_names,
                             use_case_incid = use_case_incid, 
                             npar = opt$n_prof,
                             profile = T,
                             pars_to_profile = pars_to_profile)

# Time before which we constrain R0
tvary <- dateToYears(as.Date(config$tvary))

job_rw.sd <- getRWSD(config = config, rw.sd_param = rw.sd_param, tvary = tvary,
                     profile = T, pars_to_profile = pars_to_profile)

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

save(t1, mf, file = mif_prof_filename)

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

