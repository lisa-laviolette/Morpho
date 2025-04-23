#
# Remove seasonality
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("ggforce")
library("lubridate")
library("stlplus")

load("6_incl2020.Rdata") # dcomb = zooplancton régularisé / er = environnement régularisé

# Pour chaque variable, décomposer la série temporelle en :
# - tendance (trend)
# - saisonnalité (seasonal)
# - résidu (remainder)
# La combinaison trend + remainder = signal sans saisonnalité
## Decompose all data on 2009-2017 ----
# get number of points per year
count(er, year(date)) |> count(n)
# -> most at 52

ds <- dcomb |>
  select(!TPielou:TShannon) |> #exclude when doing coponly
  gather(key="var", val="val", -date) |>
  group_by(var) |>
  do({
    s <- stlplus(.$val, t=.$date, n.p=52, s.window="periodic", t.window=52*5)
    df <- s$data[,1:4]
    df$date <- s$time
    df
  }) |>
  ungroup()

# Organisation des variables dans un ordre contrôlé pour faciliter les tracés
ds$var <- factor(ds$var, levels=c(
  "temperature", "salinity", "density", "oxygen", "no2", "no3", "po4", "sioh4", "poc", "pon", "chla",
  "conc", "nb_morphs", "MRic", "MDiv", "MEve",
  "TRic", "TShannon", "TPielou", #remove when doing coponly
  "Dim.1_mean", "Dim.2_mean", "Dim.3_mean", "Dim.4_mean", "Dim.1_sd", "Dim.2_sd", "Dim.3_sd", "Dim.4_sd", "Dim.1_q25", "Dim.2_q25", "Dim.3_q25", "Dim.4_q25", "Dim.1_q75", "Dim.2_q75", "Dim.3_q75", "Dim.4_q75"
  ))

# plot all decompositions
# Visualisation par variable des composantes décomposées (1 PDF par variable)

dst <- ds |>
  mutate(rest=trend+remainder) |>
  select(-trend, -remainder) |>
  gather(key="component", value="val", raw:rest, -date) |>
  mutate(component=factor(component, levels=unique(component)))

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


## ÉTAPE 2 : DÉCOMPOSITION DES DONNÉES ENVIRONNEMENTALES (2000–2017) ----

# Même logique de décomposition avec une fenêtre de tendance plus longue (8 ans)

es <- er |>
  gather(key="var", val="val", -date) |>
  group_by(var) |>do({
    s <- stlplus(.$val, t=.$date, n.p=52, s.window="periodic", t.window=52*8)
    df <- s$data[,1:4]
    df$date <- s$time
    df
  }) |>
  ungroup()

# Ordre contrôlé des variables environnementales
es$var <- factor(es$var, levels=c("temperature", "salinity", "density", "oxygen", "no2", "no3", "po4", "sioh4", "poc", "pon", "chla"))

# Visualisation des composantes (PDF)
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
save(ds, es, file="7_incl2020.Rdata")
