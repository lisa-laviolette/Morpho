#
# Regroup some taxa to higher level units
# Subsample zooplankton data to compensate for unequal sampling effort
# Transform morphological features to look more normal
#
# (c) 2018 GNU General Public License v3
#     Jean-Olivier Irisson, irisson@normalesup.org
#     Caroline Cailleton, caroline.cailleton@orange.fr

library("tidyverse")
source("lib_plot.R")

## Simplify taxonomy ----

z <- read_csv("data/zoo_incl2020.csv.gz", col_types=cols()) #"..._woimg" ending for current files, in case image extraction changes anything (also for "taxo_grouped file below)

#manually correct wrong date (2017-07-17-> 2018-07-17)
z<- z |> mutate(date = case_when(date(date) == "2017-07-17" ~ `year<-`(date, 2018),
                                TRUE ~ date))

# make a list of available taxa
taxo_base <- z |>
  group_by(lineage, taxon) |>
  summarise(n=n(), conc=sum(conc)) |>
  write_tsv("data/taxo_base_incl2020.tsv")

# read the current state of the taxonomic groupings
taxo_grouped <- read_csv("data/taxo_grouped - Sheet1_woimg.csv", col_types=cols()) |>
  select(lineage, level1) |>
  right_join(taxo_base, by="lineage") |>
  relocate(level1, .after=everything()) |>
  write_tsv("data/taxo_grouped2_incl2020.tsv")

# Now open taxo_grouped.tsv, copy paste it in the Google Sheet and proceed
# to the manual grouping at level1. Once this is finished, export as .csv and:

taxo_grouped <- read_csv("data/taxo_grouped - Sheet1_woimg.csv", col_types=cols()) |>
  select(taxon, level1)

# replace taxon by the grouping
z <- left_join(z, taxo_grouped, by="taxon") |>
  select(-taxon, taxon=level1) |>
  select(date, taxon, lineage, conc, objid, everything()) |>
  # remove extra taxa (which were marked as NA) = mistakes etc.
  filter(!is.na(taxon))
# compute total concentration per taxon
z |> group_by(taxon) |> summarise(conc=sum(conc)) |> arrange(conc)


## Subsample organisms to correct for sampling effort ----

# inspect number of objects on zooscan and concentration per date
effort <- z |> group_by(date) |> summarise(n=n(), conc=sum(conc))

ggplot(effort, aes(date, conc)) + geom_path() + geom_smooth(span=0.05, se=F, n=1000, colour="red") + date_x + scale_y_continuous(trans="sqrt")
# -> seasonal signal in terms of concentration

ggplot(effort) + geom_path(aes(date, n)) + date_x
# -> bias in terms of imaging effort

# remove data before 2009
z <- filter(z, date >= "2009-01-01")

# max number of objects per date
threshold <- 2000

# make random subsampling reproducible
set.seed(1)

zs <- z |> group_by(date) |> do({
  # determine how many objects to keep
  #=`threshold` or the actual number of objects when there are fewer than `threshold`
  n_rows <- nrow(.)
  n <- min(threshold, n_rows)
  # perform the random subsampling
  xs <- sample_n(., size=n, replace=FALSE)
  # compute the subsampling factor and correct concentrations
  sub_prop <- n_rows / n
  xs$conc <- xs$conc * sub_prop
  xs
})

# check that it matches what we had before
effort_s <- zs |> group_by(date) |> summarise(n=n(), conc=sum(conc))
ggplot(mapping=aes(date, conc)) +
  geom_path(data=effort, colour="black") +
  geom_path(data=effort_s, colour="red") +
  scale_y_continuous(trans="sqrt")
# -> the total matches of course, because we corrected

# per group
effort <- z |> group_by(date, taxon) |> summarise(n=n(), conc=sum(conc))
effort_s <- zs |> group_by(date, taxon) |> summarise(n=n(), conc=sum(conc))
ggplot(mapping=aes(date, conc)) +
  geom_path(data=filter(effort, taxon=="Copepoda"), colour="black") +
  geom_path(data=filter(effort_s, taxon=="Copepoda"), colour="red") +
  scale_y_continuous(trans="sqrt")
ggplot(mapping=aes(date, conc)) +
  geom_path(data=filter(effort, taxon=="Annelida"), colour="black") +
  geom_path(data=filter(effort_s, taxon=="Annelida"), colour="red") +
  scale_y_continuous(trans="sqrt")
# -> matches fairly well even for kind of rare taxa


## Transform distribution of image features ----

#' Trim extreme values
#'
#' Replace extreme values by NA
#'
#' @param x a vector
#' @param p proprotion of values to remove
#' @param side on which extreme to remove values (left=low, right=high, both=...both). Can be abbreviated
trim <- function(x, p=0.001, side="right") {
  # check argument (and allow to abbreviate it)
  side <- match.arg(side, choices=c("both", "left", "right"))
  # compute quantiles
  q <- quantile(x, probs=c(p, 1-p), na.rm=T)
  # mask extreme values
  if (side == "left" | side == "both") {
    x[x < q[1]] <- NA
  }
  if (side == "right" | side == "both") {
    x[x > q[2]] <- NA
  }
  return(x)
}

# # plot all histograms
# p <- gather(sample_frac(z, 0.1), key="var", val="val", area:perimmajor) |>
#   ggplot() + geom_histogram(aes(x=val), bins=50) + facet_wrap(~var, scales="free")
# ggsave(p, "plots/histograms_of_features.pdf", width=20, height=10)

# transform features to look more normal
zt <- z |> mutate(
  area = log10(trim(area)),
  mean = trim(mean, side="left"),
  stddev = trim(stddev),
  mode = trim(mode, side="both"),
  min = trim(min),
  max = trim(max, side="left"),
  perim. = log10(trim(perim.)),
  width = log10(trim(width)),
  height = log10(trim(height)),
  major = log10(trim(major)),
  minor = log10(trim(minor)),
  circ. = trim(circ.),
  feret = log10(trim(feret)),
  intden = log10(trim(intden)),
  median = trim(median, side="left"),
  skew = trim(skew, side="both"),
  kurt = trim(kurt),
  `%area` = log1p(trim(`%area`)),
  area_exc = log10(area_exc),
  fractal = trim(fractal, side="both"),
  skelarea = log10(skelarea),
  slope = log10(trim(slope)),
  histcum1 = trim(histcum1, side="left"),
  histcum2 = trim(histcum2, side="left"),
  histcum3 = trim(histcum3, side="left"),
  nb1 = trim(nb1, 0.005),
  nb2 = trim(nb2, 0.005),
  nb3 = trim(nb3, 0.005),
  symetrieh = log10(trim(symetrieh)),
  symetriev = log10(trim(symetriev)),
  symetriehc = trim(symetriehc),
  symetrievc = trim(symetrievc),
  convperim = log10(trim(convperim)),
  convarea = log10(trim(convarea)),
  fcons = trim(fcons),
  thickr = log1p(trim(thickr, side="both")),
  elongation = log10(trim(elongation)),
  range = trim(range, side="both"),
  meanpos = trim(meanpos, side="left"),
  cv = trim(cv),
  sr = trim(sr, side="both"),
  perimferet = trim(perimferet, 0.005),
  perimmajor = log10(trim(perimmajor))
)

# # plot all histograms
# p <- gather(sample_frac(z, 0.2), key="var", val="val", area:perimmajor) |>
#   ggplot() + geom_histogram(aes(x=val), bins=50) + facet_wrap(~var, scales="free")
# ggsave(p, "plots/histograms_of_features_normalised.pdf", width=20, height=10)

# eliminate extreme indviduals
# = more than 5 features are NA
n_na <- select(zt, area:perimmajor) |> apply(1, function(x) {sum(is.na(x))})
zt <- zt[n_na<=5,]

# replace NAs by the mean of the column
sum(!complete.cases(zt)) #count number of ind with missing values

for (col in names(select(zt, area:perimmajor))) {
  zt[[col]][is.na(zt[[col]])] <- mean(zt[[col]], na.rm=TRUE)
}

## Define characteristics of the zooplankton data set ----

z <- zt

nrow(z)
# [1] 587059 #587061(when Miriam re-tried the first time) -> with updated years 2018/2019 (&_woimg) 784104 -> incl 2020: 845812

ncol(select(z, area:perimmajor))
# [1] 45

# save to disk
save(z, file="1_incl2020.Rdata")
