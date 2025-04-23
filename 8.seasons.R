#
# Inspect seasonality among years
#OBJECTIF : Visualiser la saisonnalité annuelle moyenne et ses variations
#            dans les séries temporelles du zooplancton et de l’environnement
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("lubridate")
library("castr")
source("lib_plot.R")

load("7_incl2020.Rdata") # ds = zooplancton (décomposition STL) / es = environnement

## ÉTAPE 1 : Lisser légèrement la composante saisonnière reformat and smooth the seasonal component
ds <- ds |> select(-raw, -trend, -remainder) |>
  group_by(var) |> mutate(seasonal=slide(seasonal, k=2, n=1, fun=mean, na.rm=T)) |> ungroup() |> # Moyenne glissante (lissage)
  spread(key="var", val="seasonal")  # passage au format large pour faciliter la visualisation

## ÉTAPE 2 : Superposer les années sur un même calendrier fictif ----

# On ramène toutes les années à 2010 pour superposer les cycle
# create a new fake time coordinate to gather all data over one "virtual" year
ds$ydate <- ds$date
year(ds$ydate) <- 2010
ds$year <- factor(year(ds$date))  # conservation de l'année d'origine

## ÉTAPE 3 : Visualiser une année de saisonnalité moyenne (ex: 2009) ----

# Exemple simple avec les variables principales de la communauté

select(
  # Zooplankton Community
  ds, year, ydate, "Zoo. concentration"=conc, TRic, TDiv=TShannon, TEve=TPielou, MRic, MEve, MDiv) |>
  # Environmental variables
  #ds, year, ydate, "Temperature (°C)"=temperature, "Salinity (psu)"=salinity, "Density (g/cm^3)"=density, "Oxygen (mL/L)"=oxygen, "NO2 (mL/L)"=no2, "NO3 (mL/L)"=no3, "PO4 (mL/L)"=po4, "SiOH4 (mL/L)"=sioh4, "POC (μg/L)"=poc, "PON (μg/L)"=pon,  "Chl a (μg/L)"=chla) |>
  # Morphological traits (mean & sd)
  #ds, year, ydate, PC1=Dim.1_mean, PC2=Dim.2_mean, PC3=Dim.3_mean, PC4=Dim.4_mean, PC1.sd=Dim.1_sd, PC2.sd=Dim.2_sd, PC3.sd=Dim.3_sd, PC4.sdy=Dim.4_sd) |>
  filter(year=="2009") |> #with periodic decomposition  # Affichage d’une année donnée (modifiable)
  gather(key="Variable", val="val", -ydate, -year) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val)) +
  facet_wrap(.~Variable, scales="free_y", ncol=1) +
  scale_x_date(date_labels = "%b") +
  xlab("Month") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5))  +
  ggtitle("Zooplankton community")
ggsave("plots/season_smallpanel_env.pdf", width=8, height=6) #8x4 for 2/4 panels, for all env 8x8, axes: 10x4


## ÉTAPE 4 : Afficher plusieurs types de variables dans une même figure ----
#plot with multiple variables per panel (env/ morpho div/ taxo div/ zoo) as for manuscript
library(patchwork)
plot_list <- ds |> select(year, ydate,
                    `Zoo. concentration`=conc, MRichness=MRic, MEveness=MEve, MDivergence=MDiv,
                    Temperature=temperature, "PO4"=po4, PON=pon, "NO3"=no3, `Chl. a`=chla,
                    TRichness=TRic, TDiversity=TShannon, TEveness=TPielou,
                    #Size=Dim.1_mean, Transparency=Dim.2_mean, Circularity=Dim.3_mean, Complexity=Dim.4_mean, #for supp
                    #Size.sd=Dim.1_sd, Transparency.sd=Dim.2_sd, Circularity.sd=Dim.3_sd, Complexity.sd=Dim.4_sd, #for supp
                    #"Salinity"=salinity, "Density"=density, "Oxygen"=oxygen, "NO2"=no2, "POC"=poc, "SiOH4"=sioh4 #for supp
                    ) |>
  mutate_if(is.numeric, scale) |> # Standardisation pour comparaison
  filter(year=="2009") |>
  gather(key="Variable", val="val", -ydate, -year) |>
  mutate(vartype= case_when(Variable %in% c("MRichness", "MEveness", "MDivergence") ~ "Morphological diversity",
                            Variable %in% c("TRichness", "TDiversity", "TEveness") ~ "Taxonomic diversity",
                            Variable %in% c("Zoo. concentration", "Chl. a") ~ "Phyto- and zooplankton",
                            Variable %in% c("Temperature", "PO4", "NO3", "PON")  ~ "Environmental variables",
                            #Variable %in% c("Size", "Transparency", "Circularity", "Complexity") ~ "Morphological traits (mean)", #For supp
                            #Variable %in% c("Size.sd", "Transparency.sd", "Circularity.sd", "Complexity.sd") ~ "Morphological traits (sd)", #For supp
                            #Variable %in% c("Salinity", "Density", "Oxygen", "NO2", "POC", "SiOH4")  ~ "Environmental variables", #For supp
                            TRUE ~ "NA"))  |>
  mutate(across(vartype, factor, levels=c("Morphological diversity","Taxonomic diversity","Phyto- and zooplankton", "Environmental variables"))) |> #to change order of panels (complete community)
  mutate(across(Variable, factor, levels=c("MRichness","MEveness","MDivergence", "TRichness", "TEveness", "TDiversity", "Chl. a", "Zoo. concentration", "Temperature", "NO3", "PO4", "PON"))) |> #to change order vars within panels
  #mutate(across(vartype, factor, levels=c("Morphological diversity","Morphological traits"))) |> #to change order of panels (coponly)
  #mutate(across(Variable, factor, levels=c("Zoo. concentration", "MRichness","MEveness","MDivergence", "Size", "Transparency", "Circularity", "Complexity"))) |> #to change order vars within panels (coponly)
  #mutate(across(vartype, factor, levels=c("Morphological traits (mean)", "Morphological traits (sd)", "Environmental variables"))) |> # to change order of panels (supp)
  #mutate(across(Variable, factor, levels=c("Size", "Transparency", "Circularity", "Complexity","Size.sd", "Transparency.sd", "Circularity.sd", "Complexity.sd", "Salinity", "Density", "Oxygen", "NO2", "POC", "SiOH4"))) |> #to change order vars within panels (supp)
  group_split(vartype)

gg_list = lapply(plot_list, function(x) {
  ggplot(x) +
    geom_path(aes(x=ydate, y=val, linetype=Variable)) +
    facet_wrap(.~vartype, scales="free_y") +
    scale_x_date(date_labels = "%b") +
    xlab("Month") +
    ggtitle("") +
    ylab("") +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.title= element_blank(),
      legend.position="top") +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  })

wrap_plots(gg_list, ncol=2) + plot_annotation(title="Seasonal signal of morphological traits and environment")
ggsave("plots/season_smallpanels_smoothn1_supp_linetype.pdf", width=8, height=5) #8x4 for 2/4 panels, for all env 8x8, axes: 10x4

## ÉTAPE 5 : Variables environnementales complémentaires (pour suppléments) ----
#figure with remaining env variables for Supp:
select(ds, year, ydate, "Salinity (psu)"=salinity, "Density (g/cm^3)"=density, "Oxygen (mL/L)"=oxygen, "NO2 (mL/L)"=no2, "NO3 (mL/L)"=no3, "SiOH4 (mL/L)"=sioh4 ) |>
  mutate_if(is.numeric, scale) |>
  filter(year=="2009") |> #with periodic decomposition
  gather(key="Variable", val="val", -ydate, -year) |>
  ggplot() +
  geom_path(aes(x=ydate, y=val, linetype=Variable)) +
  scale_x_date(date_labels = "%b") +
  xlab("Month") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.title = element_text(hjust = 0.5))  +
  ggtitle("Environmental variables")
ggsave("plots/season_smallpanel_envSUP_linetype.pdf", width=6, height=4) #8x4 for 2/4 panels, for all env 8x8, axes: 10x4
