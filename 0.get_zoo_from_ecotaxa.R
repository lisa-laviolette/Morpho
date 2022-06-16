#
# Check data extracted from EcoTaxa and add associated metadata
#
# (c) 2018 GNU General Public License v3
#     Jean-Olivier Irisson, irisson@normalesup.org

library("tidyverse")
# install_github("jiho/ecotaxar")
library("ecotaxar")

# connect to EcoTaxa
db <- db_connect_ecotaxa()
tbl_ecotaxa(db)

# list PtB WP2 project
proj_ids <- ids <- c(292, 293, 294, 295, 297, 300, 301, 302, 303, 304, 337)

# # check that all samples are fully validated
# objects |>
#   filter(projid %in% proj_ids) |>
#   group_by(sampleid) |> summarise(percent_valid=count(classif_qual=="V")/n()) |>
#   filter(percent_valid < 1)
# # -> OK

# # check zooscan models
# acqs <- acquisitions |>
#   filter(projid %in% proj_ids) |> collect() |>
#   group_by(projid) |> do({
#     map <- filter(projs, projid == .$projid[1])$mappingacq
#     map_names(., map) |> select(projid, instrument, hardware, software, imgtype, scan_date)
#   }) |> ungroup() |>
#   left_join(select(projs, projid, title))
# count(acqs, title, software)

# get metadata mappings for each project
projs <- filter(projects, projid %in% proj_ids) |> collect()
# check that all mappings are the same
maps <- projs |> select(starts_with("mapping")) |> distinct()
# -> only different in process
maps <- maps[1,]

# get objects and zooprocess features
obj <- filter(objects, projid %in% proj_ids[1]) |>
  select(objid, date=objdate, classif_id, classif_qual, n01:n69, projid, sampleid, acquisid, processid)
  map_names(maps$mappingobj)

# get metadata: sampled volume, motodata fraction, pixel resolution
zoo <- obj |>
  # get info to compute concentration
  left_join(
    filter(samples, projid %in% proj_ids) |> collect() |>
      map_names(maps$mappingsample) |>
      select(projid, sampleid, orig_sampleid=orig_id, tot_vol)
  ) |>
  left_join(
    filter(acquisitions, projid %in% proj_ids) |> collect() |>
      map_names(maps$mappingacq) |>
      select(projid, acquisid, sub_part)
  ) |>
  # # get info to compute features in real world measures
  # left_join(
  #   filter(process, projid %in% proj_ids) |> collect() |>
  #     map_names(maps$mappingprocess) |>
  #     select(projid, processid, particle_pixel_size_mm, img_resolution)
  # ) |>
  # mutate_at(vars(tot_vol, sub_part, particle_pixel_size_mm, img_resolution), as.numeric)
  mutate_at(vars(tot_vol, sub_part), as.numeric)

# get taxonomic classification
taxo <- extract_taxo(db, zoo$classif_id)
zoo <- mutate(zoo,
  taxon = taxo_name(classif_id, taxo, unique=TRUE),
  lineage = lineage(classif_id, taxo)
)

# cleanup useless records
zoo <- zoo |>
  filter(
    # not validated
    classif_qual == "V",
    # or living
    str_detect(lineage, "^living")
  )

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

# remove some  taxa meaningless for morphological analysis
zoo <- zoo |> filter(
  # not relevant taxomically
  !str_detect(lineage, "Hexapoda"),
  !str_detect(taxon, "othertocheck"), !str_detect(taxon, "seaweed"),
  # not sampled quantitatively
  !str_detect(lineage, "Trachymedusae"),
  # incomplete
  !str_detect(taxon, "part$"), !str_detect(taxon, "multiple$"),
  !str_detect(taxon, "head$"), !str_detect(taxon, "tail$"),
  !str_detect(taxon, "colony$"),
  !str_detect(taxon, "egg$")
)

# compute individual "concentration"
zoo <- mutate(zoo, conc = 1 * sub_part / tot_vol) |>
  select(-sub_part, -tot_vol)

# reorder columns
zoo <- select(zoo, date, taxon, lineage, conc, everything())

write_csv(zoo, "data/zoo.csv.gz")
# save(zoo, file="0.Rdata")
