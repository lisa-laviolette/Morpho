#
# Remove seasonality
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("ggforce")
library("lubridate")
library("stlplus")

load("6.Rdata")

## Decompose all data on 2009-2017 ----

# get number of points per year
count(er, year(date)) %>% count(n)
# -> most at 52

ds <- dr %>%
  gather(key="var", val="val", -date) %>%
  group_by(var) %>% do({
    s <- stlplus(.$val, t=.$date, n.p=52, s.window=5, t.window=52*5)
    df <- s$data[,1:4]
    df$date <- s$time
    df
  }) %>%
  ungroup()

# force order of variables
ds$var <- factor(ds$var, levels=c(
  "temperature", "salinity", "density", "oxygen",
  "no2", "no3", "po4", "sioh4", "poc", "pon",
  "fluo", "chla",
  "conc", "nb_morphs",
  "MRic", "MDiv", "MEve",
  "TRic", "TShannon", "TPielou"
))

# plot all decompositions
dst <- ds %>% gather(key="component", value="val", raw:remainder) %>%
  mutate(component=factor(component, levels=unique(component)))

pdf("plots/decomposition.pdf", width=8, height=8)
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

es <- er %>%
  gather(key="var", val="val", -date) %>%
  group_by(var) %>% do({
    s <- stlplus(.$val, t=.$date, n.p=52, s.window=5, t.window=52*8)
    df <- s$data[,1:4]
    df$date <- s$time
    df
  }) %>%
  ungroup()

# force order of variables
es$var <- factor(es$var, levels=intersect(levels(ds$var), unique(es$var)))

# save to disk
save(ds, es, file="7.Rdata")
