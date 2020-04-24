
rc <- list(length=length(scn_dirs))
for (i in 1:length(scn_dirs)) {
  rc[[i]] <- load_hosp_sims_filtered(scn_dirs[i],
                                     num_files = 1000,
                                     name_filter = '',
                                     #post_process = hosp_post_process,
                                     #geoid_len = geoid_len,
                                     #padding_char = padding_char
                                     )
  rc[[i]]$scenario_num <- i
  rc[[i]]$scenario_label <- scenario_labels[[i]]
}

max1 <- rc[[1]] %>% select(icu_curr, geoid, time, scenario_label, sim_num) %>% group_by(geoid, sim_num) %>% mutate(max = max(icu_curr)) %>% ungroup() %>% filter(time == as.Date("2020-01-31")) %>% group_by(geoid) %>% summarize(quants = quantile(max, probs = .5))
max2 <- rc[[2]] %>% select(icu_curr, geoid, time, scenario_label, sim_num) %>% group_by(geoid, sim_num) %>% mutate(max = max(icu_curr)) %>% ungroup() %>% filter(time == as.Date("2020-01-31")) %>% group_by(geoid) %>% summarize(quants = quantile(max, probs = .5))
max3 <- rc[[3]] %>% select(icu_curr, geoid, time, scenario_label, sim_num) %>% group_by(geoid, sim_num) %>% mutate(max = max(icu_curr)) %>% ungroup() %>% filter(time == as.Date("2020-01-31")) %>% group_by(geoid) %>% summarize(quants = quantile(max, probs = .5))