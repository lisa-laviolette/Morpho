#
# Extracted labels, images and metadata from EcoTaxa
#
# (c) 2018 GNU General Public License v3
#     Jean-Olivier Irisson, irisson@normalesup.org

# NB: Several steps are long; run as job in RStudio.
 OBJECTIF : Télécharger, nettoyer et préparer les métadonnées et images
#            de zooplancton pour des analyses morphologiques

## CHARGEMENT DES PACKAGES ----
suppressMessages(library("tidyverse")) # Pour manipulation de données (ggplot2, dplyr, etc.)
suppressMessages(library("furrr")) # Pour parallélisation avec future + purrr
plan(multisession, workers=20)

# remotes::install_github("jiho/ecotaxar")
library("ecotaxar")
# remotes::install_github("jiho/morphr")
library("morphr") # Fonctions pour analyser l'espace morphologique (PCA, etc.)

# # Répertoire local pour stocker les images(outside of this repo, since they are big)
img_dir <- "~/datasets/pointB_wp2/"

## ÉTAPE 1 : CONNEXION À ECOTAXA ----

message("Download objects metadata from EcoTaxa") # ----

# connect to EcoTaxa
db <- db_connect_ecotaxa()

# list PtB WP2 project
proj_ids <- ids <- c(292, 293, 294, 295, 297, 300, 301, 302, 303, 304, 337, 756, 1608, 2710)

## ÉTAPE 2 : HARMONISATION DES MÉTADONNÉES ----

# # check that all samples are fully validated
# tbl(db, "objects") |>
#   filter(projid %in% proj_ids) |>
#   group_by(sampleid) |> summarise(percent_valid=count(classif_qual=="V")/n()) |>
#   filter(percent_valid < 1)
# # -> OK

# # check zooscan models
# acqs <- tbl(db, "acquisitions") |>
#   filter(projid %in% proj_ids) |> collect() |>
#   group_by(projid) |> do({
#     map <- filter(projs, projid == .$projid[1])$mappingacq
#     map_names(., map) |> select(projid, instrument, hardware, software, imgtype, scan_date)
#   }) |> ungroup() |>
#   left_join(select(projs, projid, title))
# count(acqs, title, software)

# get metadata mappings for each project
projs <- filter(tbl(db, "projects"), projid %in% proj_ids) |> collect()
# check that all mappings are the same
maps <- projs |> select(starts_with("mapping")) |> distinct()
# -> only different in process, keep one
maps <- maps[1,]

# get objects for all projects of interest
zoo <- tbl(db, "objects") |>
  filter(projid %in% !!proj_ids) |>
  # get zooprocess features
  select(projid, sampleid, acquisid, processid, objid, date=objdate, classif_id, classif_qual, n01:n69) |>
  map_names(maps$mappingobj) |>
  # get images
  left_join(
    tbl(db, "images") |>
      select(objid, img_path=file_name),
    by="objid"
  ) |>
  # get info to compute concentration (volume and fractionning rate)
  left_join(
    tbl(db, "samples") |>
      map_names(maps$mappingsample) |>
      select(sampleid, orig_sampleid=orig_id, tot_vol),
    by="sampleid"
  ) |>
  left_join(
    tbl(db, "acquisitions") |>
      map_names(maps$mappingacq) |>
      select(acquisid, sub_part, scan_date),
    by="acquisid"
  ) |>
  # # get info to compute features in real world measures
  # left_join(
  #   tbl(db, "process") |>
  #     map_names(maps$mappingprocess) |>
  #     select(processid, particle_pixel_size_mm, img_resolution),
  #   by="processid"
  # ) |>
  collect()

# reformat some columns
zoo <- zoo |>
  # mutate_at(vars(tot_vol, sub_part, particle_pixel_size_mm, img_resolution), as.numeric)
  mutate_at(vars(tot_vol, sub_part), as.numeric) |>
  mutate(scan_date = scan_date |> na_if("nan") |> parse_date("%Y%m%d"))

# get taxonomic classification
taxo <- tbl(db, "taxonomy") |> collect()
zoo <- mutate(zoo,
  taxon = taxo_name(classif_id, taxo, unique=TRUE),
  lineage = lineage(classif_id, taxo)
)

db_disconnect_ecotaxa(db)


message("Cleanup objects table") # ----

# cleanup useless records
nrow(zoo) #1410757
zoo <- zoo |>
  filter(
    # not validated
    classif_qual == "V",
    # not living
    str_detect(lineage, "^living")
  )
# count(zoo, taxon) |> arrange(taxon) |> print(n=200)

## ÉTAPE 3 : EXTRACTION DES OBJETS ----

zoo <- zoo |> filter(
    # remove some  taxa meaningless for morphological analysis
    # not relevant for plankton studies
    !str_detect(lineage, "Hexapoda"),
    !str_detect(taxon, "seaweed"),
    !str_detect(taxon, "othertocheck"),
    # not sampled quantitatively
    !str_detect(lineage, "Trachymedusae"),
    # incomplete
    !str_detect(taxon, "^part"), !str_detect(taxon, "^trunk"),
    !str_detect(taxon, "^head"), !str_detect(taxon, "^tail"),
    !str_detect(taxon, "^scale"),
    !str_detect(taxon, "^nucleus"), !str_detect(taxon, "^endostyle"),
    # representing several individuals
    !str_detect(taxon, "^colony"), !str_detect(taxon, "^chain"),
    !str_detect(taxon, "^multiple"),
    # not necessarily holoplankton
    !str_detect(taxon, "^egg"),
    # badly imaged
    !str_detect(taxon, "^badfocus")
  )
nrow(zoo)


zoo <- zoo |>
  # useless columns
  select(-classif_id, -classif_qual, -lat_end, -lon_end, -projid, -sampleid, -acquisid, -processid, -orig_sampleid) |>
  # descriptors that are
  select(
    # not meaningful
    -x, -y, -xm, -ym, -bx, -by, -angle, -xstart, -ystart, -tag,
    # or incorrect
    -esd, -centroids, -perimareaexc, -feretareaexc, -circex, -cdexc
  )

# futher remove 0 variance descriptors
sds <- summarise_if(zoo, is.numeric, sd)
sds[which(sds == 0)]
zoo <- select(zoo, -starts_with("comp"))

# compute individual "concentration"
zoo <- mutate(zoo, conc = 1 * sub_part / tot_vol) |>
  select(-sub_part, -tot_vol)

# reorder columns
zoo <- select(zoo, date, taxon, lineage, conc, everything())


message("Download images from EcoTaxa") # ----

# remove objects for which the image seems corrupted
zoo <- filter(zoo, ! objid %in% c(78785686, 34834712))

# define the location of images
orig_img_dir <- str_c(img_dir, "orig")
dir.create(orig_img_dir, showWarnings=FALSE, recursive=TRUE)
cropped_img_dir <- str_c(img_dir, "cropped")
dir.create(cropped_img_dir, showWarnings=FALSE, recursive=TRUE)

## ÉTAPE 4 : FORMATAGE ET TAXONOMIE ----

# Conversion de colonnes en valeurs numériques (volume, fraction)
zoo <- zoo |>
  mutate(
    source_img=str_c("/remote/ecotaxa/vault/", img_path),
    orig_img=str_c(orig_img_dir, "/", objid, ".jpg") |> path.expand(),
    cropped_img=str_c(cropped_img_dir, "/", objid, ".png") |> path.expand()
  )

# remove extra original images
existing_orig_imgs <- list.files(orig_img_dir, full=TRUE)
to_remove <- existing_orig_imgs[!existing_orig_imgs %in% zoo$orig_img]
message("  remove ", length(to_remove), " extra original images")
ok <- unlink(to_remove)

# get (missing) original images (in parallel)
zoo_to_copy <- filter(zoo, ! orig_img %in% existing_orig_imgs)
message("  copy ", nrow(zoo_to_copy), " original images from EcoTaxa")
ok <- future_map2_lgl(
  .x=zoo_to_copy$source_img, .y=zoo_to_copy$orig_img,
  ~file.copy(.x, .y),
  .progress=TRUE
)
all(ok)


# remove extra cropped images
existing_cropped_imgs <- list.files(cropped_img_dir, full=TRUE)
to_remove <- existing_cropped_imgs[!existing_cropped_imgs %in% zoo$cropped_img]
message("  remove ", length(to_remove), " extra cropped images")
ok <- unlink(to_remove)

# crop (missing) images (in parallel)
zoo_to_crop <- filter(zoo, ! cropped_img %in% existing_cropped_imgs)
message("  crop ", nrow(zoo_to_crop), " original images")
ok <- future_walk2(
  .x=zoo_to_crop$orig_img, .y=zoo_to_crop$cropped_img,
  function(x, y) {
    message(x)
    img_read(x) |>
      # remove legend
      img_chop(b=31) |>
      # keep only largest object
      img_extract_largest(threshold=0.1/255, quiet=TRUE) |>
      # write as png
      img_write(y)
  },
  .progress=TRUE
)

message("Save objects to disk") # ----

zoo |>
  select(-source_img, -img_path) |>
  write_csv("data/zoo.csv.gz")
# save(zoo, file="0.Rdata")
