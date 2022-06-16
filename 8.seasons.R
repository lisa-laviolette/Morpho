#
# Inspect seasonality among years
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("lubridate")
library("castr")
source("lib_plot.R")


load("7.Rdata")

# reformat and smooth the seasonal component
ds <- ds |> select(-raw, -trend, -remainder, -deseason) |>
  group_by(var) |> mutate(seasonal=slide(seasonal, k=2, n=2, fun=mean, na.rm=T)) |> ungroup() |>
  spread(key="var", val="seasonal")

# create a new fake time coordinate to gather all data over one "virtual" year
ds$ydate <- ds$date
year(ds$ydate) <- 2010
ds$year <- factor(year(ds$date))


# TODO rename things from the get go, it would be easier
select(ds, year, ydate, Richness=MRic, Divergence=MDiv, Evenness=MEve) |>
  gather(key="var", val="val", -ydate, -year) |>
  mutate(var=factor(var, levels=unique(var))) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val, colour=year)) +
  facet_grid(var~., scales="free_y") +
  x_season + theme(axis.title=element_blank())
ggsave("plots/season_indices.pdf", width=6, height=8)

select(ds, year, ydate, NO3=no3, `Chl a`=chla, `Zooplankton`=conc) |>
  gather(key="var", val="val", -ydate, -year) |>
  mutate(var=factor(var, levels=unique(var))) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val, colour=year)) +
  facet_grid(var~., scales="free_y") +
  x_season + theme(axis.title=element_blank())
ggsave("plots/season_env.pdf", width=6, height=8)

