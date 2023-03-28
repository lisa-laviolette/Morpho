#
# Put time series on a regular time coordinate
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("castr")

load("2_incl2020.Rdata") # env
load("5_incl2020.Rdata") # zoo

# define a new regular coordinate
reg_dates <- seq(from=min(e$date)-1, to=max(e$date), by=7)
# Nb: shift by 1 day because we get more matches in env data
sum(reg_dates %in% e$date)


## Regularise environemental data ----

# 1. leaving NAs when data is not there
ew <- spread(e, key="var", value="val") |>
  # compute closest regular date
  mutate(
    reg_date = closest(from=date, within=reg_dates),
    diff_date = abs(reg_date - date) |> as.numeric()
  ) |>
  # remove points too far from a regular date
  filter(diff_date <= 3) |>
  # assume the closest date on the regular coordinate as the new coord
  select(-date, -diff_date) |>
  rename(date=reg_date)
# add missing values when the data is too far from a regular date
reg <- tibble(date=reg_dates)
er <- left_join(reg, ew, by="date")

# 2. interpolating through NAs
ef <- e |>
  group_by(var) |>
  do({
    tibble(
      date=reg_dates,
      var=.$var[1],
      val=interpolate(.$date, .$val, xout=reg_dates, method="linear")
    )
  }) |>
  ungroup() |>
  spread(key="var", val="val")

# plot two versions
ggplot(gather(er, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")
ggplot(gather(ef, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")


## Regularise zooplankton data ----

reg_dates <- reg_dates[reg_dates >= "2009-01-01"]
reg_dates <- reg_dates[reg_dates < "2021-01-01"]

# 1. leaving NAs when data is not there
dw <- d |>
  # compute closest regular date
  mutate(
    reg_date = closest(from=date, within=reg_dates),
    diff_date = abs(reg_date - date) |> as.numeric()
  ) |>
  # remove points too far from a regular date
  filter(diff_date <= 3) |>
  # assume the closest date on the regular coordinate as the new coord
  select(-date, -diff_date) |>
  rename(date=reg_date)

# add missing values when the data is too far from a regular date
reg <- tibble(date=reg_dates)
dr <- left_join(reg, dw, by="date")

# 2. interpolating through NAs
df <- gather(d, key="var", val="val", -date) |>
  group_by(var) |>
  do({
    tibble(
      date=reg_dates,
      var=.$var[1],
      val=interpolate(.$date, .$val, xout=reg_dates, method="linear")
    )
  }) |>
  ungroup() |>
  spread(key="var", val="val")


# 3. check size of NA gabs and fill those that are small (1, 2 weeks); leave others NA (option used in manuscript)
gaps<- df |>
  #mutate(week= rep(1:52, 11),
  #       year= year(gaps$date)) |>
  left_join(dr, by="date") |>
  arrange(date) |>
  mutate(gapID = data.table::rleid(is.na(nb_morphs.y)))

gaps2<- gaps |>  filter(is.na(nb_morphs.y)) |>
  count(gapID, name = "gapsize") |>
  mutate(todo=case_when(gapsize <=2 ~ "fill",
                        gapsize > 2 ~"leaveNA"))
gapsNA<- gaps2$gapID[gaps2$todo=="leaveNA"]

dcomb<- gaps |>
  mutate(conc = ifelse(gapID %in% gapsNA, conc.y, conc.x),
         MDiv = ifelse(gapID %in% gapsNA, MDiv.y, MDiv.x),
         MEve = ifelse(gapID %in% gapsNA, MEve.y, MEve.x),
         MRic = ifelse(gapID %in% gapsNA, MRic.y, MRic.x),
         nb_morphs = ifelse(gapID %in% gapsNA, nb_morphs.y, nb_morphs.x),
         TPielou = ifelse(gapID %in% gapsNA, TPielou.y, TPielou.x),
         TRic = ifelse(gapID %in% gapsNA, TRic.y, TRic.x),
         TShannon = ifelse(gapID %in% gapsNA, TShannon.y, TShannon.x),
         Dim.1_mean = ifelse(gapID %in% gapsNA, Dim.1_mean.y, Dim.1_mean.x),
         Dim.2_mean = ifelse(gapID %in% gapsNA, Dim.2_mean.y, Dim.2_mean.x),
         Dim.3_mean = ifelse(gapID %in% gapsNA, Dim.3_mean.y, Dim.3_mean.x),
         Dim.4_mean = ifelse(gapID %in% gapsNA, Dim.4_mean.y, Dim.4_mean.x),
         Dim.1_sd = ifelse(gapID %in% gapsNA, Dim.1_sd.y, Dim.1_sd.x),
         Dim.2_sd  = ifelse(gapID %in% gapsNA, Dim.2_sd.y, Dim.2_sd.x),
         Dim.3_sd  = ifelse(gapID %in% gapsNA, Dim.3_sd.y, Dim.3_sd.x),
         Dim.4_sd  = ifelse(gapID %in% gapsNA, Dim.4_sd.y, Dim.4_sd.x),
         Dim.1_q25 = ifelse(gapID %in% gapsNA, Dim.1_q25.y, Dim.1_q25.x),
         Dim.2_q25 = ifelse(gapID %in% gapsNA, Dim.2_q25.y, Dim.2_q25.x),
         Dim.3_q25 = ifelse(gapID %in% gapsNA, Dim.3_q25.y, Dim.3_q25.x),
         Dim.4_q25 = ifelse(gapID %in% gapsNA, Dim.4_q25.y, Dim.4_q25.x),
         Dim.1_q75 = ifelse(gapID %in% gapsNA, Dim.1_q75.y, Dim.1_q75.x),
         Dim.2_q75 = ifelse(gapID %in% gapsNA, Dim.2_q75.y, Dim.2_q75.x),
         Dim.3_q75 = ifelse(gapID %in% gapsNA, Dim.3_q75.y, Dim.3_q75.x),
         Dim.4_q75 = ifelse(gapID %in% gapsNA, Dim.4_q75.y, Dim.4_q75.x),
  ) |>
  arrange(date) |>
  select(date, conc, MDiv, MEve, MRic, nb_morphs, TPielou, TRic, TShannon,
         Dim.1_mean, Dim.2_mean, Dim.3_mean, Dim.4_mean, Dim.1_sd, Dim.2_sd, Dim.3_sd, Dim.4_sd,
         Dim.1_q25, Dim.2_q25, Dim.3_q25, Dim.4_q25, Dim.1_q75, Dim.2_q75, Dim.3_q75, Dim.4_q75)

#gaps<- subset(dr, is.na(dr$nb_morphs))
#gaps$year<-year(gaps$date)
#gaps$MonDay<-paste(month(gaps$date), day(gaps$date), sep="-")

# plot again
ggplot(gather(dr, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")
ggplot(gather(df, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")
ggplot(gather(dcomb, var, val, -date)) + geom_path(aes(date, val), na.rm=T) + facet_wrap(~var, scales="free_y")

## Combine and save ----
dcomb<-left_join(dcomb, ef) #or decide which option to take
save(dr, df, dcomb, er, ef, file="6_incl2020.Rdata")
