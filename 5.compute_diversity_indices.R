#
# Compute morphological and taxonomic diversity indexes for all dates
# Plot and compare them
#
# (c) 2018 GNU General Public License v3
#     Jean-Olivier Irisson, irisson@normalesup.org
#     Caroline Cailleton, caroline.cailleton@orange.fr

library("tidyverse")
library("castr")
source("lib_plot.R")

# for functional diversity indexes
source("lib_functional_diversity.R")

# for taxonomic diversity indexes
library("vegan")

load("4_incl2020.Rdata")

## Compute functional diversity indexes ----

# weights matrix: date x morphs
weights <- z |>
  # concentration per date per morph
  group_by(date, morph_nb) |> summarise(conc=sum(conc)) |>
  # convert to wide format and fill with 0s
  spread(morph_nb, conc, fill=0) |>
  as.data.frame()
# convert to matrix with named rows
rownames(weights) <- as.character(weights$date)
weights <- select(weights, -date)
weights <- as.matrix(weights)

# morphological space matrix: morphs x axes
space <- centers |>select(Dim.1:Dim.4) |> as.matrix()
rownames(space) <- centers$morph_nb

# compute index
fd <- multidimFD4(space, weights, verb=TRUE)

# convert result into data.frame
d <- as.data.frame(fd)
d$date <- row.names(d) |> as.Date()
d <- as.tibble(d) |>
  select(date, nb_morphs=Nb_sp, conc=Tot_weight, MRic=FRic, MDiv=FDiv, MEve=FEve)


## Compute taxonomic diversity indexes ----

# compute community matrix: date x taxon
comm <- z |>
  group_by(date, taxon) |> summarise(conc=sum(conc)) |>
  spread(key=taxon, value=conc, fill=0) |>
  as.data.frame()
rownames(comm) <- as.character(comm$date)
comm <- select(comm, -date)
comm <- as.matrix(comm)

# compute taxonomic diversity indexes
d <- mutate(d,
  TRic = specnumber(comm),
  TShannon = diversity(comm, index="shannon"),
  TPielou = TShannon / log(TRic)
  )


## Inspect indexes ----

# plot time series
dt <- gather(d, key="index", val="val", conc:TPielou)
dt$index <- factor(dt$index, levels=unique(dt$index))
ggplot(dt) +
  geom_path(aes(x=date, y=val)) + facet_grid(index~., scales="free_y")

# plot only morphological indices
d |> rename(Richness=MRic, Divergence=MDiv, Evenness=MEve) |>
  gather(key="index", val="val", Richness:Evenness) |>
  mutate(index=factor(index, levels=unique(index))) |>
ggplot() +
  geom_point(aes(x=date, y=val), colour="grey50", size=0.5) +
  geom_path(aes(x=date, y=slide(val, k=3, mean, n=3, na.rm=T)), colour="red", size=0.75) +
  # geom_smooth(aes(x=date, y=val), colour="red", size=0.75, se=F, n=500, span=0.04) +
  theme(axis.title.y=element_blank()) +
  #x_year +
  facet_grid(index~., scales="free_y")
ggsave("plots/morpho_indexes_series_incl2020_coponly.pdf", width=8, height=5)


# compare taxo and morphological indexes
ggplot(d) +
  geom_point(aes(x=TRic/max(TRic), y=MRic), alpha=0.5, shape=16) +
  geom_abline(slope=1, intercept=0)
cor.test(d$TRic, d$MRic, method="spearman")

ggplot(d) +
  geom_point(aes(x=TShannon, y=MDiv), alpha=0.5, shape=16) +
  geom_abline(slope=1, intercept=0)
cor.test(d$TShannon, d$MDiv, method="spearman")

ggplot(d) +
  geom_point(aes(x=TPielou, y=MEve), alpha=0.5, shape=16) +
  geom_abline(slope=1, intercept=0)
cor.test(d$TPielou, d$MEve, method="spearman")


#add also mean coordinates along PCA axes to df
ax<- z |> group_by(date) |> summarise_at(vars(Dim.1:Dim.4), funs(mean, sd, quantile(., 0.25), quantile(., 0.75)))
colnames(ax)<-c("date", "Dim.1_mean", "Dim.2_mean", "Dim.3_mean", "Dim.4_mean", "Dim.1_sd", "Dim.2_sd", "Dim.3_sd", "Dim.4_sd", "Dim.1_q25", "Dim.2_q25", "Dim.3_q25", "Dim.4_q25", "Dim.1_q75", "Dim.2_q75", "Dim.3_q75", "Dim.4_q75")
d<- left_join(d, ax, by="date")

save(d, file="5_incl2020_coponly.Rdata")
