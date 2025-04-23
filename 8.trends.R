#
# Detect trends in environmental data
## OBJECTIF : Identifier et quantifier les tendances temporelles dans les
#            séries environnementales et de diversité du zooplancton
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("trend") # Pour le test de Mann-Kendall
library("nlme") # Pour la régression linéaire généralisée (GLS)
library("broom") #Pour extraire les statistiques des modèles


load("7_incl2020.Rdata")
## FONCTIONS UTILES ----

# Extraire les stats utiles d’un modèle GLS
glance.gls <- function(m) {
  s <- summary(m)

  # r.squared
  f <- predict(m)
  mss <- sum((f - mean(f))^2)
  rss <- sum(residuals(m)^2)
  rsq <- mss / (mss + rss)

  # residuals
  shap <- shapiro.test(m$residuals)
  a <- pacf(residuals(m, type="normalized"), plot=FALSE)

  tibble(
    r.squared = rsq,

    # sigma = s$sigma,

    statistic = s$tTable[2, "t-value"],
    p.value = s$tTable[2, "p-value"],

    # AIC = s$AIC,
    # BIC = s$BIC,

    intercept = m$coefficients[1],
    slope = m$coefficients[2],

    shapiro.p.value = shap$p.value,
    cor.struct=class(m$modelStruct$corStruct)[1],
    acf1 = a$acf[1],
    acf2 = a$acf[2]
  )
}
# Convertir une p-value en étoiles de significativité
signif_stars <- function(x) {
  case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

## ANALYSE DES TENDANCES ----
# combine the two periods
es$period <- "2000-"
ds$period <- "2009-"

s <- bind_rows(es, ds) # Fusion des deux sources
s$var <- factor(s$var, levels=levels(ds$var))

# Recalcul de la série désaisonnalisée
s <- mutate(s, deseason = trend+ifelse(is.na(remainder), 0, remainder))

# Calcul statistique pour chaque variable et période
# x <- filter(ds, var == "pon")
stats <- s |> group_by(period, var) |>
  do({
    x <- .
    # message(x$var[1], " ", nrow(x))

    # 1. Mann-Kendall trend test
    mkt <- mk.test(x$deseason)

    # 2. GLS regression
    # simple model
    m <- gls(deseason ~ date, data=x) 
      # Détection d’autocorrélation significative

    a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    # if autocorrelation is too strong
    if (abs(a$acf[1]) > 0.2) {
      # add AR1 model on residuals
    # Modèle GLS avec erreur AR1
      m <- gls(deseason ~ date, data=x, cor=corAR1(round(a$acf[1], 1)))
      a <- pacf(residuals(m, type="normalized"), plot=FALSE)
      # # if autocorrelation is still too strong
      # # NB: does not bring a much better fit and is very long
      # if (abs(a$acf[1]) > 0.2) {
      #   # add AR2 model on residuals
      #   cor <- round(a$acf[1:2], 1)
      #   # stabilise errors if the autocorr is too strong
      #   if (sum(cor) > 0.9){ cor <- cor - 0.1 }
      #   m <- gls(deseason ~ date, data=x, cor=corARMA(cor, form=~date, p=2, q=0))
      #   # a <- pacf(residuals(m, type="normalized"), plot=FALSE)
      # }
    }

    # extract diagnostic information for both approaches
    bind_cols(
      glance(mkt) |> select(p.value), # Résultat Mann-Kendall
      glance.gls(m) |> select(r.squared, p.value, intercept, slope, cor.struct, acf=acf1)
     )
  }) |>
  ungroup() |>
  mutate(
    acf = abs(acf),# force valeur absolue de l’autocorr
    signif.mk = signif_stars(p.value...3),
    signif.gls = signif_stars(p.value...5)
  ) |>
  rename(p.value.mk = p.value...3,
         r.squared.mk = r.squared,
         p.value.gls = p.value...5,
         intercept.gls = intercept,
         slope.gls = slope,
         cor.struct.gls = cor.struct,
         acf.gls = acf
         )
# Export des résultats

write_tsv(stats, "plots/stats_periodic_incl2020_ext_coponly.tsv")

## Plots ----

# Sélection de variables à illustrer
vars <- c("Zoo. concentration (ind/L)"="2009-conc", "Temperature (°C)"="2000-temperature", `Morph. Richness`="2009-MRic", "Chl a (μg/L)"="2000-chla", `Morph. Eveness`="2009-MEve", "POC (μg/L)"="2000-poc", `Morph. Divergence`="2009-MDiv", "PON (μg/L)"="2000-pon")
vars <- c("Concentration (ind/L)"="2009-conc", `Morph. Richness`="2009-MRic", `Morph. Eveness`="2009-MEve", `Morph. Divergence`="2009-MDiv") #coponly SUPP
# Filtrage des données à tracer

dp <- s |>
  mutate(id=str_c(period, var)) |>
  filter(id %in% vars) |>
  mutate(id=names(vars)[match(id, vars)]) |>
  mutate(id=factor(id, levels=names(vars))) |>
  filter(date > "2009-01-01")

statsp <- stats |>
  mutate(id=str_c(period, var)) |>
  filter(id %in% vars) |>
  mutate(id=names(vars)[match(id, vars)]) |>
  mutate(id=factor(id, levels=names(vars)))

#when drawing 2nd line for recent years (environmental variables)
#(have to also change period in first var-selection to '2000-')
#vars2 <- c("Zoo. concentration (ind/L)"="2009-conc", "Temperature (°C)"="2009-temperature", `Morph. Richness`="2009-MRic", "Chl a (μg/L)"="2009-chla", `Morph. Eveness`="2009-MEve", "POC (μg/L)"="2009-poc", `Morph. Divergence`="2009-MDiv", "PON (μg/L)"="2009-pon")

#statsp2 <- stats |>
#  mutate(id=str_c(period, var)) |>
#  filter(id %in% vars2) |>
#  mutate(id=names(vars2)[match(id, vars2)]) |>
#  mutate(id=factor(id, levels=names(vars2)))
# Tracé : série + droite de tendance

ggplot() +
  facet_wrap(id~., scales="free_y", ncol=2) +
  geom_path(aes(date, deseason), data=dp, colour="grey20") + # série désaisonnalisée
  geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp, signif.gls %in% c("*", "**", "***")), colour="red", size=0.75, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ggtitle("Long-term trends") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=11),
        panel.grid.minor.y=element_blank())
ggsave("plots/trends_smallpanels.pdf", width=7, height=9)
