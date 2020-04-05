
##'Function to plot county final sizes
##'
##'@param country_final_dat sf object with county final sizes
##'
##'@return ggplot object, final sizes with CI (y) by county (x)
##'
make_final_plot <- function(country_final_dat){
  p <- ggplot(country_final_dat,
              aes(x=reorder(county.name,-median), y=median,ymin=pi_low, ymax=pi_high )) +
    geom_pointrange() +
    ylab("Infections") +
    xlab("County") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}