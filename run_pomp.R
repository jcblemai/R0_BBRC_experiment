source("COVID-pomp/scripts/utils.R")

option_list = list(
  optparse::make_option(c("-c", "--config"), action="store", default='pomp_config.yaml', type='character', help="path to the config file"),
  optparse::make_option(c("-j", "--jobs"), action="store", default=parallel::detectCores(), type='numeric', help="number of cores used"),
  optparse::make_option(c("-r", "--run_level"), action="store", default=1, type='numeric', help="run level for MIF"),
  optparse::make_option(c("-n", "--nfilter"), action="store", default=10, type='numeric', help="Number of filtering iterations"),
  optparse::make_option(c("-w", "--downweight"), action="store", default=0, type='numeric', help="downweight ikelihood to be used for filtering")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
config <- load_config(opt$c)

sh_script <- readLines("run-all_pomp_template.sh")
sh_script <- gsub("%l", paste(lapply(config$likelihoods$to_compute, function(x) paste0("'", x, "'")), collapse = " "), sh_script)
sh_script <- gsub("%c", paste(lapply(config$places, function(x) paste0("'", x, "'")), collapse = " "), sh_script)
sh_script <- gsub("%r", opt$run_level, sh_script)
sh_script <- gsub("%j", opt$jobs, sh_script)
sh_script <- gsub("%n", opt$nfilter, sh_script)
sh_script <- gsub("%w", opt$downweight, sh_script)

write(sh_script, file = "run-all_pomp.sh")

system("nohup bash run-all_pomp.sh &")

# Post-process 
# source("COVID-pomp/scripts/postprocess_results.R")

