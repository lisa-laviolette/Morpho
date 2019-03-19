#
# Put time series on a regular time coordinate
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("castr")

load("2.Rdata") # env
load("5.Rdata") # zoo

# define a new regular coordinate
reg_dates <- seq(from=min(e$date)-1, to=max(e$date), by=7)
# Nb: shift by 1 day because we get more matches in env data
sum(reg_dates %in% e$date)


## Regularise environemental data ----

# 1. leaving NAs when data is not there
ew <- spread(e, key="var", value="val") %>%
  # compute closest regular date
  mutate(
    reg_date = closest(from=date, within=reg_dates),
    diff_date = abs(reg_date - date) %>% as.numeric()
  ) %>%
  # remove points too far from a regular date
  filter(diff_date <= 3) %>%
  # assume the closest date on the regular coordinate as the new coord
  select(-date, -diff_date) %>%
  rename(date=reg_date)
# add missing values when the data is too far from a regular date
reg <- tibble(date=reg_dates)
er <- left_join(reg, ew, by="date")

# 2. interpolating through NAs
ef <- e %>%
  group_by(var) %>%
  do({
    tibble(
      date=reg_dates,
      var=.$var[1],
      val=interpolate(.$date, .$val, xout=reg_dates, method="linear")
    )
  }) %>%
  ungroup() %>%
  spread(key="var", val="val")

# plot two versions
ggplot(gather(er, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")
ggplot(gather(ef, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")


## Regularise zooplankton data ----

reg_dates <- reg_dates[reg_dates >= "2009-01-01"]

# 1. leaving NAs when data is not there
dw <- d %>%
  # compute closest regular date
  mutate(
    reg_date = closest(from=date, within=reg_dates),
    diff_date = abs(reg_date - date) %>% as.numeric()
  ) %>%
  # remove points too far from a regular date
  filter(diff_date <= 3) %>%
  # assume the closest date on the regular coordinate as the new coord
  select(-date, -diff_date) %>%
  rename(date=reg_date)
# add missing values when the data is too far from a regular date
reg <- tibble(date=reg_dates)
dr <- left_join(reg, dw, by="date")

# 2. interpolating through NAs
df <- gather(d, key="var", val="val", -date) %>%
  group_by(var) %>%
  do({
    tibble(
      date=reg_dates,
      var=.$var[1],
      val=interpolate(.$date, .$val, xout=reg_dates, method="linear")
    )
  }) %>%
  ungroup() %>%
  spread(key="var", val="val")

# plot again
ggplot(gather(dr, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")
ggplot(gather(df, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")

## Combine and save ----

dr <- left_join(dr, er)
df <- left_join(df, ef)

save(dr, df, er, ef, file="6.Rdata")
