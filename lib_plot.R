library("ggplot2")

x_year <- list(scale_x_date(date_breaks="1 year", date_labels="%Y"), theme(panel.grid.minor.x=element_blank()))
x_season <- scale_x_date(date_breaks="2 months", date_labels="%b")
date_x <- x_season
