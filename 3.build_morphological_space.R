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

load("1_incl2020.Rdata")

#filter only Copepoda for coponly analysis
#z<- z |> filter(taxon=="Copepoda")

## ÉTAPE 1 : CALCUL DE LA PCA PONDÉRÉE ----
pcaw <- PCA(
  select(z, area:perimmajor), # variables morphologiques
  row.w=z$conc/sum(z$conc), # NB: weights must be scaled to sum to 1# pondération selon concentration
  graph=FALSE, ncp=5
)

### ÉTAPE 2 : EXPLORATION DES RÉSULTATS ----

# How many axes to keep?
# Visualisation de la variance expliquée (screeplot) Screeplot: visualization of % of explained variance per axis
library(factoextra)
fviz_screeplot(pcaw, addlabels = TRUE, ylim = c(0, 50))

# Kaiser-Guttman (keep axes with eigenvalues >1)
pcaw$eig

# Affichage des contributions des variables aux axes
plot(pcaw, choix="var", cex=0.5)
dev.print(pdf, "plots/pca_axes12.pdf", width=8, height=6)
plot(pcaw, choix="var", cex=0.5, axes=4:5)
dev.print(pdf, "plots/pca_axes45_anc.pdf", width=8, height=6)
# -> Interpretation
# Axis 1: circular vs. large
# Axis 2: variable in grey level vs. light
# Axis 3: elongated and symetric vs. complex perimeter
# Axis 4: more difficult but circular vs. complex perimeter and symetric

# Visualisation des individus dans le plan PCA (axes 1 et 2)
plot(Dim.2 ~ Dim.1, data=pcaw$ind$coord, pch=".", asp=1)

# color per taxon
plot(Dim.2 ~ Dim.1, data=pcaw$ind$coord, pch=".", asp=1, col=factor(z$taxon))


# coloured plot with selected variables
#to colour arrows by feature-type
cols<-read_csv2("data/traitdescriptors.csv")
cols<- as.data.frame(cols)
library(stringr)
cols$trait2<- str_to_title(cols$trait2)

#to select vars based on contribution to each plane
contribs<-data.frame(pcaw$var$contrib); contribs$feature<- rownames(contribs)
contribs<- contribs %>% mutate(Dim.12 = abs(Dim.1) + abs(Dim.2),
                               Dim.23 = abs(Dim.2) + abs(Dim.3),
                               Dim.34 = abs(Dim.3) + abs(Dim.4))
# Préparation des coordonnées des variables (flèches)
pca_coor<-as.data.frame(pcaw$var$coord)
pca_coor$name2<- rownames(pca_coor)
pca_coor<-left_join(pca_coor, cols[,c(2, 6, 7)], by=c("name2"="name2"))

# Biplot avec flèches colorées selon le type de trait
fviz_pca_var(pcaw, axes=c(1,2), #change axes here
             select.var = list(name =rownames(contribs[order(contribs$Dim.12, decreasing = TRUE),])[1:12]), #change axes/plane here!
             habillage = factor(pca_coor$trait2),
             labelsize=5,
             palette = c("royalblue", "tomato", "goldenrod", "turquoise")
             ) +
  labs(x=paste("PC1 ", round(pcaw$eig[1,2], 1), "%", sep=""), #select axes here!
       y=paste("PC2 ", round(pcaw$eig[2,2], 1), "%", sep=""),
       title=""
       ) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.print(pdf, "plots/pca_colvars_ax12.pdf", width=8, height=6)

# check distribution of some variables in taxa
z |>group_by(taxon) |> summarise(
  cv=mean(cv),
  mean=mean(mean),
  area=mean(area),
  major=mean(major)) |>
  arrange(mean)

## ÉTAPE 4 : AJOUT DES COORDONNÉES PCA AUX DONNÉES ----

# add coordinates in morphological space to zooplankton data
z <- bind_cols(z, as.tibble(pcaw$ind$coord[,1:4]))

# save to disk
save(z, file="3_incl2020.Rdata")

## ÉTAPE 5 : AFFICHAGE DES IMAGES DANS L’ESPACE MORPHOLOGIQUE ----
#visual representation of morphospace
#-> ! requires the actual images !
# Définir les chemins vers les images originales
img_dir <- "~/datasets/pointB_wp2/orig"
z <- z |> mutate(img_path=str_c(img_dir, "/", objid, ".jpg"))
#sum(file.exists(z$img_path))

# Prétraitement des images (rognage + amélioration du contraste)
preprocess <- function(x) {
  x |>
    # remove 31 pixels from the bottom (=the scale bar)
    img_chop(bottom=31) |>
    # change the gamma to see the light objects better
    img_adjust_gamma(gamma=0.7)
}

# To draw a circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
circ <- circleFun(c(0,0),2,npoints = 500)

#to filter most important axes (highest contribution to plane)
contribs<-data.frame(pcaw$var$contrib); contribs$feature<- rownames(contribs)
contribs<- contribs %>% mutate(Dim.12 = abs(Dim.1) + abs(Dim.2),
                               Dim.23 = abs(Dim.2) + abs(Dim.3),
                               Dim.34 = abs(Dim.3) + abs(Dim.4))
var_contrib <- rownames(contribs[order(contribs$Dim.12, decreasing = TRUE),])[1:12]

# Homogénéisation des échelles pour le biplot
#homogenize scaling between individuals & variables for correct biplot
# Change scaling of variables/columns from scaling 1 to 2
eig <- pcaw$eig[,1]
varScores <- as.data.frame(t(t(pcaw$var$coord) / sqrt(eig))) #de-scale
varScores2 <- as.data.frame(t(t(varScores) * sqrt(nrow(varScores) * eig))) #re-scale
pca_coor<-varScores2
pca_coor$name<- rownames(pca_coor)

#to colour features by feature-type
cols<-read_csv2("data/traitdescriptors.csv")
cols<- as.data.frame(cols)
library(stringr)
cols$trait2<- str_to_title(cols$trait2)
pca_coor<-left_join(pca_coor, cols[,c(2, 7)], by=c("name"="name2"))

# Affichage final : morphospace avec images et flèches

ggmorph_tile(pcaw, z$img_path, steps=18, n_imgs=3, fun=preprocess, dimensions=c(1,2)) +
  geom_path(data=circ, aes(x*7,y*7), lty=2, color = "grey", alpha = 0.7) + #adapt scaling of circle to fit the arrows (here: 7 for plane 12; 5 for plane 34)
  geom_hline(yintercept=0, lty=2, color="grey", alpha=0.9) +
  geom_vline(xintercept=0, lty=2, color="grey", alpha=0.9) +
  geom_segment(data=subset(pca_coor, pca_coor$name %in% var_contrib), aes(x=0, xend=Dim.1, y=0, yend=Dim.2, color=factor(trait2)),
               arrow=arrow(length=unit(0.025, "npc"), type = "open"),
               lwd=1.2) +
  geom_text(data=subset(pca_coor, pca_coor$name %in% var_contrib), aes(x=Dim.1*1.05, y=Dim.2*1.05, label=name, color=factor(trait2)),
            lwd=5, hjust=0) +
  scale_color_manual(values=c("royalblue", "tomato", "goldenrod", "turquoise"), name="Type") +
  #ggtitle("(aucune transformation)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="top",
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))
ggsave("plots/morphospace_images_scaling2_ax12_circ_coponly.png", width=10, height=8, dpi=350)
