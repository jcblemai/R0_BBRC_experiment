
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
           deceased = str_trim(decede) == "X",
           transfert = `transfert_autre_hopital_nom/localite`
    )
  
  # df$uid <- with(df, str_c(birthday, premieres_lettres_nom, premieres_lettres_prenom)) %>% 
  # lapply(digest::digest) %>% 
  # unlist()
  
  return(select(df, one_of(cols)))
}

# Hospitalizatoin data  ------------------------------------------------------------------------

hfiles <- dir("docs/VD/", pattern = "*SII", full.names = T)
hfiles_dates <- str_extract(hfiles, "(?<=_)[0-9]+(?=_|\\.)") %>% as.Date("%Y%m%d")

sel_cols <- c("uid", "age", "date_in", "date_out", "icu_in", "icu_out", "deceased", "transfert")

hosp_data <- foreach(hfile = hfiles,
                     hfile_date = hfiles_dates,
                     .combine = rbind) %do% {
                       
                       hdate <- as.character(hfile_date) %>% 
                         str_replace_all("-", "")
                       
                       hosp_in <- readxl::read_excel(hfile, sheet = hdate) %>% 
                         extractHospData(sel_cols) %>% 
                         mutate(sheet = "in")
                       hosp_out <- readxl::read_excel(hfile, sheet = str_c(hdate, " sortie")) %>%
                         extractHospData(cols = sel_cols) %>% 
                         mutate(sheet = "out")
                       
                       rbind(hosp_in, hosp_out) %>% 
                         mutate(report = hdate)
                     } 

# Remove duplicates
hosp_data <- distinct(hosp_data) %>% filter(!(icu_in < date_in) | is.na(icu_in))
hosp_data$deceased[is.na(hosp_data$date_out)] <- NA

# Compute total durations for each patient
clean_data <- hosp_data %>%
  filter(!is.na(date_out) | !is.na(icu_in) | !is.na(icu_out)) %>%
  select(-report) %>%
  distinct() %>%
  arrange(uid) %>%
  group_by(uid, icu_in) %>%
  
  # Remove repeated information on ICUS
  arrange(icu_out) %>% slice(1) %>%
  ungroup()

clean_data <- clean_data %>%
  rbind(hosp_data %>%
          filter(is.na(date_out), is.na(icu_in), !(uid %in% clean_data$uid)) %>%
          select(-report) %>%
          distinct())


# number of ICU stays
clean_data %>%
  filter(!is.na(icu_in)) %>%
  group_by(uid) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# number of hosp stays
clean_data %>%
  distinct(uid, date_in) %>%
  group_by(uid) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


# transfers
clean_data %>%
  filter(sheet == "out", !is.na(transfert)) %>%
  distinct(uid, transfert) %>%
  group_by(transfert) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


hosp_data %>%
  group_by(uid) %>%
  summarise(n= n())

hosp_out <- hosp_data %>% 
  filter(!is.na(uid)) %>% 
  group_by(uid) %>% 
  filter(!is.na(date_out))

hosp_current <- hosp_data %>% 
  filter(!(uid %in% hosp_out$uid)) %>% 
  group_by(uid) %>% 
  arrange(desc(report)) %>% 
  slice(1)

hosp_data %>%
  group_by(uid) %>%
  summarise(n= n())

