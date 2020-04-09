flattenNames <- function(str) {
  str <- iconv(str,from="UTF-8",to="ASCII//TRANSLIT") %>%
    str_to_lower() %>% 
    str_replace_all(" ", "_") %>% 
    str_replace_all("\\*|\\(|\\)", "") %>% 
    str_replace_all("^[0-9]_", "")
  return(str)
}
