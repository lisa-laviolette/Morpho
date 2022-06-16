#
# Groups individuals in PCA space into 200 clusters => morphs
#
# (c) 2018 GNU General Public License v3
#     Jean-Olivier Irisson, irisson@normalesup.org
#     Caroline Cailleton, caroline.cailleton@orange.fr

library("tidyverse")
# devtools::install_github("jiho/wkmeans")
library("wkmeans")

# library("scales")
# require("grid")
# library("png")
# library("dplyr")
# library("ggplot2")

load("3.Rdata")


## Create and inspect morphs ----

# divide into 200 morphs through observation-weighted k-means
# morphs <- wkmeans(select(z, Dim.1:Dim.4), k=200, w=z$conc, iter_max=100, nstart=50, cores=20)
# save(morphs, file="morphs.Rdata")
load("morphs.Rdata")

# add morph number to the full zooplankton data
z$morph_nb <- factor(morphs$cluster)

# get the center of each morph
centers <- as_tibble(morphs$centers)
# compute total concentration per morph
tot_conc <- z |> group_by(morph_nb) |>
  summarise(conc=sum(conc))
# combine the two
centers <- bind_cols(centers, tot_conc)

# plot morphs centers in PCA space
ggplot(centers) +
  geom_point(aes(Dim.1, Dim.2, size=conc), alpha=0.5, shape=16) +
  scale_size(range=c(1,10))
# -> most in the middle but that's OK

# plot morphs in PCA space
ggplot(z) + coord_fixed() +
  geom_point(aes(Dim.1, Dim.2), shape=".") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0)
ggsave("plots/pca.png", width=8, height=6)

ggplot(z) + coord_fixed() +
  geom_point(aes(Dim.1, Dim.2, colour=morph_nb), shape=".") +
  scale_colour_discrete(guide="none") +
  scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0)
ggsave("plots/pca_morphs.png", width=8, height=6)

# Are morphs representative of taxa ?
conc_per_morph <- z |> group_by(morph_nb, taxon) |> summarise(conc=sum(conc))
ggplot(conc_per_morph) + geom_col(aes(x=morph_nb, y=conc, fill=taxon)) + scale_fill_discrete(guide="none")
ggplot(filter(conc_per_morph, taxon != "Copepoda")) + geom_col(aes(x=morph_nb, y=conc, fill=taxon)) + scale_fill_discrete(guide="none")
# -> the morphs do not contain only one species

# # are the morphs present on all dates ?
# tab <- dmesures_sna |> select(objid, date, g_kmeansw) |>
#   group_by(date) |> summarise(nb=n_distinct(g_kmeansw))
# tab$date <- as.Date(tab$date)
#
# ggplot(tab) + geom_path(aes(x=date, y=nb)) + scale_x_date(date_breaks="1 year", labels=date_format("%Y"))
#
# # number of morphs below 100 for all dates
# points <- tab$nb
# supprlignes <- 0
# for (i in 1:nrow(tab)){
#   if (points[i]>100){
#     supprlignes <- c(supprlignes, i)
#   }
# }
# supprlignes <- supprlignes[-1]
# tab2 <- tab[-supprlignes, ]
# tab2
#
# #### 4) Compute concentration for each date ####
# flute <- z |> select(date, conc) |>
#   group_by(date) |> summarise(conc=sum(conc))
# tab$conc <- flute$conc
# ggplot(tab) + geom_path(aes(date, conc))

save(z, centers, file="4.Rdata")
