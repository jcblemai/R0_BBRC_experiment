
# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(foreach)


flattenNames <- function(str) {
  str <- iconv(str,from="UTF-8",to="ASCII//TRANSLIT") %>%
    str_to_lower() %>% 
    str_replace_all(" ", "_") %>% 
    str_replace_all("\\*|\\(|\\)", "") %>% 
    str_replace_all("^[0-9]_", "")
  return(str)
}

extractHospData <- function(df, cols) {
  df <- set_colnames(df, flattenNames(colnames(df))) %>% 
    mutate(uid = no_interne,
           birthday = as.Date(date_de_naissane, format = "%d.%m.%Y"),
           age = difftime(Sys.Date(), birthday, units = "days"),
           age = round(age/365),
           date_in = as.Date(arrivee_hopital_le_date, format = "%d.%m.%Y"),
           date_out = as.Date(sortie_hopital_date, format = "%d.%m.%Y"),
           icu_in = as.Date(debut_soins_intensifs_date, format = "%d.%m.%Y"),
           icu_out = as.Date(fin_soins_intensifs_date, format = "%d.%m.%Y"),
           deceased = str_trim(decede) == "X") 
  
  # df$uid <- with(df, str_c(birthday, premieres_lettres_nom, premieres_lettres_prenom)) %>% 
  # lapply(digest::digest) %>% 
  # unlist()
  
  return(select(df, one_of(cols)))
}

# Hospitalizatoin data  ------------------------------------------------------------------------

hfiles <- dir("analysis_CH/docs/VD/", pattern = "*SII", full.names = T)
hfiles_dates <- str_extract(hfiles, "(?<=_)[0-9]+(?=_|\\.)") %>% as.Date("%Y%m%d")

sel_cols <- c("uid", "age", "date_in", "date_out", "icu_in", "icu_out", "deceased")

hosp_data <- foreach(hfile = hfiles,
                     hfile_date = hfiles_dates,
                     .combine = rbind) %do% {
                       
                       hdate <- as.character(hfile_date) %>% 
                         str_replace_all("-", "")
                       
                       hosp_in <- readxl::read_excel(hfile, sheet = hdate) %>% 
                         extractHospData(sel_cols)
                       
                       hosp_out <- readxl::read_excel(hfile, sheet = str_c(hdate, " sortie")) %>%
                         extractHospData(cols = sel_cols)
                       
                       rbind(hosp_in, hosp_out) %>% 
                         mutate(report = hdate)
                     } 

hosp_out <- hosp_data %>% 
  filter(!is.na(uid)) %>% 
  group_by(uid) %>% 
  filter(!is.na(date_out))

hosp_current <- hosp_data %>% 
  filter(!(uid %in% hosp_out$uid)) %>% 
  group_by(uid) %>% 
  arrange(desc(report)) %>% 
  slice(1)

hosp_data <- rbind(hosp_out, hosp_current) %>% 
  mutate(outcome = case_when(deceased ~ "deceased",
                             !is.na(date_out) ~ "discharged",
                             T ~ "hospitalized"))

write_csv(hosp_data, path = "analysis_CH/data/VD/hospitalization_data.csv")



