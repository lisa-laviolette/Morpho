#
# Build the table of environmental data from 2009 to 2017
# Integrate environmental data between 0 and 75 meter
#
# (c) 2018 Caroline Cailleton, caroline.cailleton@orange.fr
# last update : 2018-07-01
#

library("tidyverse")
library("lubridate")
# devtools::install_github("jiho/castr")
library("castr")

# read data
env <- read_csv("data/radehydro_ctd.csv", col_types=cols(), guess_max=100000) |>
  mutate(date=as.Date(date_time))

# eliminate data with bad qc
envq <- select(env, date, depth, starts_with("q_"))
envv <- select(env, str_replace(names(envq), "q_", ""))
env <- left_join(
  gather(envv, key=var, val=val, -date, -depth),
  gather(envq, key=var, val=qc, -date, -depth) |> mutate(var=str_replace(var, "q_", ""))
)
env$val[! env$qc %in% c(2, 5, 21, 22)] <- NA
env <- select(env, -qc)

# get pH
ph <- read_csv("data/radehydro_pH.csv", col_types=cols()) |>
  mutate(date=as.Date(sampling_date)) |>
  # keep only relevant variables
  select(date, depth, pH, pCO2) |>
  # reduce replicates
  group_by(date, depth) |> summarise_all(mean, na.rm=T) |> ungroup() |>
  #filter updated PH data to then only bind new data
  filter(date > "2017-07-25")  |>
  # convert to tall format
  gather(key="var", value="val", pH, pCO2) |>
# add it to the rest
env <- bind_rows(env, ph)


# subsample data
env <- env |>
  # keep only some years
  filter(date >= "2000-01-01", date < "2020-01-01") |> #changed to updates data 2018/2019
  # and some variables
  filter(var %in% c("temperature", "salinity", "sigma_theta", "fluorescence", "oxygen_winkler", "no3", "no2", "po4", "sioh4", "poc", "pon", "chla", "pH", "pCO2"))

# change some variable names for simplicity
env <- spread(env, var, val) |>
  rename(density=sigma_theta, oxygen=oxygen_winkler, fluo=fluorescence) |>
  gather(var, val, -date, -depth)

# cleanup some values
ggplot(env) + geom_histogram(aes(x=val), bins=50) + facet_wrap(~var, scales="free")
ggplot(env) + geom_point(aes(x=date, y=val, colour=depth), shape=".") + facet_wrap(~var, scales="free")
env <- spread(env, key=var, value=val)
env$salinity[env$salinity < 36] <- NA
env$density[env$salinity < 36] <- NA
env$density[env$density < 23] <- NA
env <- gather(env, key=var, value=val, -date, -depth)


# integrate all variables over depth
envi <- env |>
  group_by(date, var) |>
  summarise(val=integrate(val, depth, from=0, to=75, fun="mean")) |>
  ungroup()

# plot resulting series
ggplot(envi) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")

# further remove some strange values
envi <- spread(envi, key=var, value=val)
envi$density[envi$density < 25] <- NA
envi <- gather(envi, key=var, value=val, -date)
ggplot(envi) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")

e <- envi
save(e, file="2_woimg.Rdata")
