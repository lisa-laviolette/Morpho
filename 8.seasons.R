#
# Inspect seasonality among years
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("lubridate")
library("castr")
source("lib_plot.R")


load("7_incl2020.Rdata")

# reformat and smooth the seasonal component
ds <- ds |> select(-raw, -trend, -remainder) |> #deseason doesn't exist? removed (", -deseason")
  group_by(var) |> mutate(seasonal=slide(seasonal, k=2, n=2, fun=mean, na.rm=T)) |> ungroup() |>
  #group_by(var, seasonal) |> slice(1) |>
  spread(key="var", val="seasonal")

# create a new fake time coordinate to gather all data over one "virtual" year
ds$ydate <- ds$date
year(ds$ydate) <- 2010
ds$year <- factor(year(ds$date))


# TODO rename things from the get go, it would be easier
select(#ds, year, ydate, Richness=MRic, Divergence=MDiv, Evenness=MEve) |>
  #ds, year, ydate, Richness=TRic, Shannon=TShannon, Pielou=TPielou) |>
  #ds, year, ydate, Temperature=temperature, Salinity=salinity, Density=density, Oxygen=oxygen, NO2=no2, NO3=no3, PO4=po4, SiOH4=sioh4, POC=poc, PON=pon,  `Chl a`=chla, Fluorescence=fluo) |>
  ds, year, ydate, "Zoo. concentration (ind/L)"=conc, "No. morphs"= nb_morphs, "Temperature (°C)"=temperature, "Salinity (psu)"=salinity, "Density (g/cm^3)"=density, "Oxygen (mL/L)"=oxygen, "NO2 (mL/L)"=no2, "NO3 (mL/L)"=no3, "PO4 (mL/L)"=po4, "SiOH4 (mL/L)"=sioh4, "POC (μg/L)"=poc, "PON (μg/L)"=pon,  "Chl a (μg/L)"=chla) |>
  #ds, year, ydate, "No. morphs"=nb_morphs, Concentration=conc) |>
  #ds, year, ydate, Ax1=Dim.1, Ax2=Dim.2, Ax3=Dim.3, Ax4=Dim.4) |>
  filter(year=="2009") |> #with periodic decomposition
  gather(key="var", val="val", -ydate, -year) |>
  mutate(var=factor(var, levels=unique(var))) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val)) + #, colour=year
  #facet_grid(var~., scales="free_y") +
  facet_wrap(.~var, scales="free_y") +
  #x_season +
  scale_x_date(date_labels = "%b") +
  xlab("Month") +
  ylab("") +
  theme_bw() +
  theme(#axis.title=element_blank()
        panel.grid.minor=element_blank())  +
  ggtitle("Seasonality of zooplankton and environmental variables")
ggsave("plots/season_envzoo_incl2020.pdf", width=8, height=8) #8x4 for 2/4 panels, for all env 8x8, axes: 10x4

select(ds, year, ydate, NO3=no3, `Chl a`=chla, PO4=po4, `Zooplankton`=conc) |>
  gather(key="var", val="val", -ydate, -year) |>
  mutate(var=factor(var, levels=unique(var))) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val, colour=year)) +
  facet_grid(var~., scales="free_y") +
  #x_season +
  theme(axis.title=element_blank()) +
  ggtitle("Copepods only")
ggsave("plots/season_env.pdf", width=6, height=8)

select(ds, year, ydate, Ax1=Dim.1, Ax2=Dim.2, Ax3=Dim.3, Ax4=Dim.4) |>
  gather(key="var", val="val", -ydate, -year) |>
  mutate(var=factor(var, levels=unique(var))) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val, colour=year)) +
  facet_wrap(var~., scales="free_y") +
  #x_season +
  theme(axis.title=element_blank()) #+
  ggtitle("No copepods, all morpho space")
ggsave("plots/season_axes_comp_hor.pdf", width=8, height=4)


#calculate intensity of seasonality
#define spring period
ssn<- ds |> select(-date, -period) |>
  gather(key="var", val="val", -year, -ydate) |>
  mutate(season= case_when(month(ydate) %in% c(2:5) ~ "spring",
                           TRUE ~ "rest")) |>
  select(-ydate) |>
  group_by(year, var, season) |>
  summarise(mean=mean(val),
            sd=sd(val),
            meanmin = mean(val) + min(val)) |>
  ungroup()

#calculate intensity
ssn_int<- ssn |>
  select(-sd, -mean) |>
  spread(key=season, val=meanmin) |>
  group_by(year, var) |>
  mutate(intensity=(spring-rest)/rest)

#test for significance (spring vs. rest of the year)
ssn_sig<- ssn |> group_by(var) |>
  summarize_at(vars(mean), ~wilcox.test(. ~ season)$p.value) |>
  mutate(sig = case_when(
    mean < 0.001 ~ "***",
    mean < 0.01  ~ "**",
    mean < 0.05  ~ "*",
    mean < 0.1   ~ ".",
    TRUE      ~ ""))
  #do(  #TRY WITH DO
  #  w<-wilcox.test(mean~season, data=.)) |>
  #summarise(var, wilcox=w$p.value)
  #bind_cols(
  #  glance(w) |> select(p.value)
  #)

write_tsv(ssn_sig, "plots/season_stats.tsv")
