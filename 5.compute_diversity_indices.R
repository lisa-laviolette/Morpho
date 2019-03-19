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
# TODO review this code!

# for taxonomic diversity indexes
library("vegan")


load("4.Rdata")

## Compute functional diversity indexes ----

# weights matrix: date x morphs
weights <- z %>%
  # concentration per date per morph
  group_by(date, morph_nb) %>% summarise(conc=sum(conc)) %>%
  # convert to wide format and fill with 0s
  spread(morph_nb, conc, fill=0) %>%
  as.data.frame()
# convert to matrix with named rows
rownames(weights) <- as.character(weights$date)
weights <- select(weights, -date)
weights <- as.matrix(weights)

# morphological space matrix: morphs x axes
space <- centers %>% select(Dim.1:Dim.4) %>% as.matrix()
rownames(space) <- centers$morph_nb

# compute index
fd <- multidimFD4(space, weights, verb=TRUE)

# convert result into data.frame
d <- as.data.frame(fd)
d$date <- row.names(d) %>% as.Date()
d <- as.tibble(d) %>%
  select(date, nb_morphs=Nb_sp, conc=Tot_weight, MRic=FRic, MDiv=FDiv, MEve=FEve)


## Compute taxonomic diversity indexes ----

# compute community matrix: date x taxon
comm <- z %>%
  group_by(date, taxon) %>% summarise(conc=sum(conc)) %>%
  spread(key=taxon, value=conc, fill=0) %>%
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
d %>% rename(Richness=MRic, Divergence=MDiv, Evenness=MEve) %>%
  gather(key="index", val="val", Richness:Evenness) %>%
  mutate(index=factor(index, levels=unique(index))) %>%
ggplot() +
  geom_point(aes(x=date, y=val), colour="grey50", size=0.5) +
  geom_path(aes(x=date, y=slide(val, k=3, mean, n=3, na.rm=T)), colour="red", size=0.75) +
  # geom_smooth(aes(x=date, y=val), colour="red", size=0.75, se=F, n=500, span=0.04) +
  theme(axis.title.y=element_blank()) +
  x_year +
  facet_grid(index~., scales="free_y")
d %>% rename(Richness=MRic, Divergence=MDiv, Evenness=MEve) %>%
  gather(key="index", val="val", Richness:Evenness) %>%
  mutate(index=factor(index, levels=unique(index))) %>%
ggplot() +
  geom_path(aes(x=date, y=val)) +
  theme(axis.title.y=element_blank()) +
  x_year +
  facet_grid(index~., scales="free_y")
ggsave("plots/morpho_indexes_series.pdf", width=8, height=5)


# compare taxo and morphological indexes
ggplot(d) +
  geom_point(aes(x=TRic/max(TRic), y=MRic), alpha=0.5, shape=16) +
  geom_abline(slope=1, intercept=0)
cor.test(d$TRic, d$MRic, method="spearman")
# S = 6297800, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.5623991

ggplot(d) +
  geom_point(aes(x=TShannon, y=MDiv), alpha=0.5, shape=16) +
  geom_abline(slope=1, intercept=0)
cor.test(d$TShannon, d$MDiv, method="spearman")
# S = 5733300, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.6016251

ggplot(d) +
  geom_point(aes(x=TPielou, y=MEve), alpha=0.5, shape=16) +
  geom_abline(slope=1, intercept=0)
cor.test(d$TPielou, d$MEve, method="spearman")
# S = 5689100, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.6046965


save(d, file="5.Rdata")

# TODO: review this
# #### II. Graphics of each morphological space ####
# tab$mois <- month(tab$date, label=T)
# tab$semaine <- week(tab$date)
# D_index$date <- as.Date(rownames(D_index))
#
# date <- D_index$date
# date <- as.character(date)
# for (i in 1:length(date)){
#   a <- date[i]
#   plotFD2(weight=weight_mat, Serie=D_index$date, Ric=D_index$FRic, Div=D_index$FDiv, Eve=D_index$FEve, date=a, mois=tab$mois, semaine=tab$semaine, coord=centers, Faxes_plot=colnames(centers)[1:5], abond=T)
#   a
# }
#
# # Creating a graph to represent morphological richness and functionnal space for each axis for two date
# plotFD3(weight=weight_mat, Serie=D_index$date, Ric=D_index$FRic, Div=D_index$FDiv, Eve=D_index$FEve, date="2012-02-22", date2="2012-03-21", mois=tab$mois, semaine=tab$semaine, coord=centers, Faxes_plot=colnames(centers)[1:5], abond=T)
#
