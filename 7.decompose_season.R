#
# Remove seasonality
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("ggforce")
library("lubridate")
library("stlplus")

load("6_incl2020.Rdata")

## Decompose all data on 2009-2017 ----
# get number of points per year
count(er, year(date)) |> count(n)
# -> most at 52

ds <- dcomb |> #decide whether to use filled nas (df), na (dr,"default"), only small gaps filled (dcomb)
  select(!TPielou:TShannon) |> #exclude when doing coponly
  gather(key="var", val="val", -date) |>
  group_by(var) |> do({
    s <- stlplus(.$val, t=.$date, n.p=52, s.window="periodic", t.window=52*5) #s.window=5, t.window=52*5 (default, when giving no t.window: 99); twin now doesn't matter because we're looking at trend + remainder combined
    df <- s$data[,1:4]
    df$date <- s$time
    df
  }) |>
  ungroup()

#to check seasonal sub-series (weekly)
# test<- dcomb |>select(date, MRic)
# s <- stlplus(dr$MRic, n.p=52, s.window=5, t.window=NULL)
# plot(s)
#
# plot_seasonal(s, col = c("darkgray", "black"), lwd = 2, xlab = "Time", ylab = "Centered Seasonal + Remainder")
#
#
# test$MRic[seq(1, 1000, by=50)]

# force order of variables
ds$var <- factor(ds$var, levels=c("temperature", "salinity", "density", "oxygen", "no2", "no3", "po4", "sioh4", "poc", "pon","fluo", "chla",
  "conc", "nb_morphs",
  "MRic", "MDiv", "MEve",
  #"TRic", "TShannon", "TPielou", #remove when doing coponly
  "Dim.1_mean", "Dim.2_mean", "Dim.3_mean", "Dim.4_mean", "Dim.1_sd", "Dim.2_sd", "Dim.3_sd", "Dim.4_sd", "Dim.1_q25", "Dim.2_q25", "Dim.3_q25", "Dim.4_q25", "Dim.1_q75", "Dim.2_q75", "Dim.3_q75", "Dim.4_q75"
))

# plot all decompositions
dst <- ds |>
  mutate(rest=trend+remainder) |>
  select(-trend, -remainder) |>
  gather(key="component", value="val", raw:rest, -date) |>
  mutate(component=factor(component, levels=unique(component))) #|>
  #filter(var %in% c("conc", "nb_morphs", "MRic")) |> droplevels() #added to check decomposition for few vars only

pdf("plots/decomp_periodic_incl2020_ext.pdf", width=8, height=8)
p <- lapply(levels(dst$var), function(v){
  print(
    ggplot(filter(dst, var==v), aes(date, val)) +
      geom_path(na.rm=TRUE) +
      facet_grid(component~., scales="free_y") +
      labs(title=v)
  )
})
dev.off()

## Decompose environmental data on 2000-2017 ----

es <- er |>
  gather(key="var", val="val", -date) |>
  group_by(var) |>do({
    s <- stlplus(.$val, t=.$date, n.p=52, s.window=5, t.window=52*8)
    df <- s$data[,1:4]
    df$date <- s$time
    df
  }) |>
  ungroup()

# force order of variables
es$var <- factor(es$var, levels=c("temperature", "salinity", "density", "oxygen", "no2", "no3", "po4", "sioh4", "poc", "pon", "fluo", "chla"))

# plot all decompositions
est <- es |>
  mutate(rest=trend+remainder) |>
  select(-trend, -remainder) |>
  gather(key="component", value="val", raw:rest, -date) |>
  mutate(component=factor(component, levels=unique(component)))

pdf("plots/decomp_env_periodic_incl2020.pdf", width=8, height=8)
p <- lapply(levels(est$var), function(v){
  print(
    ggplot(filter(est, var==v), aes(date, val)) +
      geom_path(na.rm=TRUE) +
      facet_grid(component~., scales="free_y") +
      labs(title=v)
  )
})
dev.off()


# save to disk
#ds_f<-ds
#ds_comb<-ds
save(ds, es, file="7_incl2020.Rdata")
