#
# Build a PCA that will represent the basic morphological space
# Illustrate it with plots of the aspect of objects in various sections of the PCA
#
# (c) 2018 Caroline Cailleton, caroline.cailleton@orange.fr
# last update : 2018-07-01
#

library("tidyverse")
library("FactoMineR")
# library("reticulate")
# use_python(system("which python3", intern=TRUE))
# devtools::install_github("jiho/morphr")
library("morphr")

load("1.Rdata")

## Compute PCA and plot result ----

pcaw <- PCA(
  select(z, area:perimmajor),
  row.w=z$conc/sum(z$conc), # NB: weights must be scaled to sum to 1
  graph=FALSE, ncp=5
)

# plot variables
plot(pcaw, choix="var", cex=0.5)
dev.print(pdf, "plots/pca_axes12.pdf", width=8, height=6)
plot(pcaw, choix="var", cex=0.5, axes=3:4)
dev.print(pdf, "plots/pca_axes34.pdf", width=8, height=6)
# -> Interpretation
# Axis 1: circular vs. large
# Axis 2: variable in grey level vs. light
# Axis 3: elongated and symetric vs. complex perimeter
# Axis 4: more difficult but circular vs. complex perimeter and symetric

# check distribution of some variables in taxa
z %>% group_by(taxon) %>% summarise(
  cv=mean(cv),
  mean=mean(mean),
  area=mean(area),
  major=mean(major)
) %>% arrange(mean)


# plot individuals
plot(Dim.2 ~ Dim.1, data=pcaw$ind$coord, pch=".", asp=1)
# -> OK, well spread out

# color per taxon
plot(Dim.2 ~ Dim.1, data=pcaw$ind$coord, pch=".", asp=1, col=factor(z$taxon))
# TODO redo in ggplot2
# -> taxa are somewhat clustered in space but many are still quite spread out


# # histogram of each axis
# ft <- gather(coor_ind, key="var", value="val", Axis1:Axis5)
# ggplot(sample_n(ft, 10000)) + geom_histogram(aes(x=val), bins=50)
# + facet_wrap(~var, scales="free")

# add coordinates in morphological space to zooplankton data
z <- bind_cols(z, as.tibble(pcaw$ind$coord[,1:4]))
# save to disk
save(z, file="3.Rdata")

## Represent space ----

z <- mutate(z,
  path=paste0("/home/jiho/ecotaxa_ml-raw/zooscanPtBWP2/imgs/",objid,".jpg"),
  exists=file.exists(path)
)

zz <- filter(z, exists)
pcaw <- morpho_space(
  zz %>% select(area:perimmajor),
  w=zz$conc/sum(zz$conc)
)
p <- ggmorph_tile(pcaw, imgs=zz$path, steps=13, scale=0.005, dimensions=c(1,2))
ggsave(p, filename="plots/pca_morphs12.png", width=10, height=10)
p <- ggmorph_tile(pcaw, imgs=zz$path, steps=13, scale=0.005, dimensions=c(3,4))
ggsave(p, filename="plots/pca_morphs34.png", width=10, height=10)

## Select objects in specific regions of the PCA space ----
# TODO to re-run, was not redone for SFEcologie2018 on 2018-10-23


# Toolbox functions

# Select elements at the top of x
top <- function(x, p=0.1) {
  x > quantile(x, 1-p)
}
# Select elements at the bottom of x
bottom <- function(x, p=0.1) {
  x < quantile(x, p)
}
# Select elements around a given quantile of x, by default the median
middle <- function(x, p=0.1, p.ref=0.5) {
  between(x, quantile(x, p.ref-p), quantile(x, p.ref+p))
}

#' Rotate a 2D set of points
#' @param x matrix (or object coercible to a matrix) in which the first two columns specify the coordinates of points
#' @param angle in radians to rotate the points by
rotate <- function(x, theta=0) {
  # compute rotation matrix
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), byrow=T, nrow=2)
  # perform rotation by matrix mu
  t(R %*% t(as.matrix(x[,1:2])))
}

#' Subset n points that match several criteria
#' @param X a logical matrix, each column designates whether the object satisfies a given criterion; all criteria must be matched to retain the object
#' @param n how many objects to retain
subset_n <- function(X, n=10) {
  # get indexes of rows where all criteria match
  idx <- which(rowSums(X) == ncol(X))
  # get n of them
  if (length(idx) >= n) {
    idx <- sample(idx, size=n)
  } else {
    stop("Not enough elements")
  }
  return(idx)
}


# get object coords and cos2 in PCA space
x    <- as.tibble(pcaw$ind$coord)
cos2 <- as.tibble(pcaw$ind$cos2)
# add object ids=file names
x$objid <- z$objid

# make this reproducible
set.seed(1)

# in space of axes 1 and 2
obj12 <- list()

# get objects in the middle
im <- apply(x, 2, middle) %>% subset_n()
obj12[[1]] <- x[im,c(1:2,6)]

# get objects at the extremes, rotating the space so as to extract extreme objects in various directions
directions <- seq(0, 2*pi, by=pi/4)
for (i in 2:length(directions)) {
  # rotate points
  theta <- directions[i]
  xx <- x
  xx[,1:2] <- rotate(x[,1:2], theta)
  # get elements at the top extreme of the horizontal axis
  ix <- cbind(top(xx[,1], p=0.05), apply(xx[,-1], 2, middle, p=0.35)) %>% subset_n()
  # ix <- cbind(top(xx[,1], p=0.02), top(rowSums(cos2[,1:2]), p=0.02)) %>% subset_n()
  # NB: selecting through cos2 is not restrictive enough
  obj12[[i]] <- x[ix,c(1:2,6)]
}

# in space of axes 3 and 4
# TODO functionalise this. copying and pasting the code is bad
obj34 <- list()

# get objects in the middle
# im <- apply(x, 2, middle) %>% subset_n()
obj34[[1]] <- x[im,c(3:4,6)]

# get objects at the extremes, rotating the space so as to extract extreme objects in various directions
directions <- seq(0, 2*pi, by=pi/4)
for (i in 2:length(directions)) {
  # rotate points
  theta <- directions[i]
  xx <- x
  xx[,3:4] <- rotate(x[,3:4], theta)
  # get elements at the top extreme of the horizontal axis
  ix <- cbind(top(xx[,3], p=0.05), apply(xx[,-3], 2, middle, p=0.4)) %>% subset_n()
  # ix <- cbind(top(xx[,3], p=0.05), top(rowSums(cos2[,3:4]), p=0.05)) %>% subset_n()
  obj34[[i]] <- x[ix,c(3:4,6)]
}

## Compute morphs and their positions ----

# Compute the path to the image of a given object id
img <- function(id) {
  paste0("/home/jiho/ecotaxa_ml-raw/zooscanPtBWP2/imgs/",id,".jpg")
}

# for both PCA planes
obj <- list(obj12, obj34)
obj <- lapply(obj, function(o) { lapply(o, function(x) {
  # morph the selected objects together
  m <- morph(img(x$objid), adjust_grey=T)
  # compute the morph's coordinates and dimensions
  xy <- colMeans(x[,1:2])
  wh <- dim(m)
  # return both
  list(
    coords=data.frame(x=xy[1], y=xy[2], w=ncol(m), h=nrow(m)),
    morph=m
  )
})})

plot_morphs <- function(morphs, scores=NULL, cos2=NULL, dims=1:2, scale=0.5, ...) {
  # plot scores (=coordinates) of individuals in the target PCA plane
  # plot(scores[,dims], asp=1, pch=".", col=interp_map(rowSums(cos2[,dims])))

  # plot coordinates of morphs
  coords <- map_df(morphs, `[[`, "coords")
  plot(coords[,1:2],
       # be discreet
       pch=".", bty="n",
       # force aspect ratio and expand scales
       asp=1, xlim=range(coords[,1])*1.5, ylim=range(coords[,2])*1.5,
       ...
  )

  # add 0,0 lines
  abline(h=0, v=0, col="grey80")

  # add morph images
  lapply(morphs, function(m) {
    # scale width and height
    m$coords[,3:4] <- m$coords[,3:4] * scale
    # remove white background
    x <- m$morph
    x[x>250] <- NA
    # plot image
    rasterImage(x/255, m$coords$x-m$coords$w/2, m$coords$y-m$coords$h/2, m$coords$x+m$coords$w/2, m$coords$y+m$coords$h/2)
  })
  return(invisible(NULL))
}

# Define the label for a PCA axis
pc_label <- function(i, eig) {
  paste0("PC ", i, " (", round(eig[i,2],1),"%)")
}

# plot both planes side by side (and try to keep the scale of images consistent)
par(mfrow=c(1,2))
plot_morphs(obj[[1]], xlab=pc_label(1, pcaw$eig), ylab=pc_label(2, pcaw$eig), scale=0.035)
plot_morphs(obj[[2]], xlab=pc_label(3, pcaw$eig), ylab=pc_label(4, pcaw$eig), scale=0.025)
dev.print(pdf, "plots/pca_with_morphs.pdf", width=12, height=6)



