
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(sf)
library(magrittr)
library(foreach)
library(itertools)

# Under what mobility threshold to cut
amount_thresh <- 1e2

mobility_zones <- st_read("data/ch/MNTP/Verkehrszonen_Schweiz_NPVM_2017.shp") %>% 
  as_tibble() %>% 
  select(ID, N_KT) 

mobility_matrix <- read_delim("data/ch/MNTP/DWV_PW_2017_CH.mtx",
                              # n_max = 10,
                              skip = 8, 
                              col_names = c("from", "to", "amount"),
                              trim_ws = T, delim = " ")

mobility_matrix <- inner_join(mobility_matrix, 
                              mobility_zones, by = c("from" = "ID")) %>% 
  inner_join(mobility_zones, by = c("to" = "ID"), suffix = c(".from", ".to"))


mobmat <- foreach(ms = isplit(mobility_matrix, mobility_matrix$N_KT.from),
        .combine = rbind,
        .packages = c("dplyr", "magrittr")) %do% {
          ms$value %>% 
            group_by(N_KT.to) %>% 
            summarise(amount = sum(amount)) %>% 
            mutate(N_KT.from = ms$key[[1]])
        } 

mobmatfinal <- filter(mobmat, amount > amount_thresh) %>% 
  select(N_KT.from, N_KT.to, amount) %>% 
  set_colnames(c("ori", "dest", "amount")) %>% 
  filter(!str_detect(ori, "Enk"), !str_detect(dest, "emk"))

write_csv(mobmatfinal, path = "data/ch/mobility_matrix_adj.csv")

ggplot(mobmatfinal, aes(x = ori, y = dest, fill = amount)) +
  geom_raster() +
  scale_fill_viridis_c()

