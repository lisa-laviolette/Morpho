#
# Detect trends in environmental data
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("trend")
library("nlme")
library("broom")

load("7_incl2020.Rdata")

## Utilities ----

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

signif_stars <- function(x) {
  case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

## Test trends ----


# combine the two periods
es$period <- "2000-"
ds$period <- "2009-"

s <- bind_rows(es, ds)
s$var <- factor(s$var, levels=levels(ds$var))

# compute the deseasonalised part
s <- mutate(s, deseason = trend+ifelse(is.na(remainder), 0, remainder))

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
    a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    # if autocorrelation is too strong
    if (abs(a$acf[1]) > 0.2) {
      # add AR1 model on residuals
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
      glance(mkt) |> select(p.value), #|> set_names(., paste0, ".mk"), #causes error: The size of `nm` (10) must be compatible with the size of `x` (1).
      glance.gls(m) |> select(r.squared, p.value, intercept, slope, cor.struct, acf=acf1) #|> set_names(., paste0, ".gls") #causes error:  The size of `nm` (10) must be compatible with the size of `x` (6).
     )
  }) |>
  ungroup() |>
  mutate(
    #acf.gls = abs(acf.gls),
    #signif.mk = signif_stars(p.value.mk),
    #signif.gls = signif_stars(p.value.gls)
    acf = abs(acf),
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

# View(stats)
#write_tsv(stats, "plots/trend_stats_comb.tsv")
write_tsv(stats, "plots/stats_periodic_incl2020_ext_coponly.tsv")


# TODO check whether copepods are declining more than the reste which could explain the increase of divergence

## Plots ----

vars <- c(Temperature="2000-temperature", NO3="2000-no3", `Fluo (Chl a)`="2000-fluo", `Part. org. Nitrogen`="2000-pon", Zooplankton="2009-conc", `Morph. Divergence`="2009-MDiv")
vars <- c(Temperature="2000-temperature", Salinity="2000-salinity", Density="2000-density", Oxygen="2000-oxygen", NO2="2000-no2", NO3="2000-no3", PO4="2000-po4", SiOH4="2000-sioh4", `Part. org. carbon`="2000-poc", `Part. org. nitrogen`="2000-pon", `Chl. a`="2000-chla")

vars<- c("Zoo. concentration (ind/L)"="2009-conc", "No. morphs"= "2009-nb_morphs", "Temperature (°C)"="2000-temperature", "Salinity (psu)"="2000-salinity", "Oxygen (mL/L)"="2000-oxygen", "PO4 (mL/L)"="2000-po4", "POC (μg/L)"="2000-poc", "PON (μg/L)"="2000-pon",  "Chl a (μg/L)"="2000-chla") #sig 12 years
vars<- c("Zoo. concentration (ind/L)"="2009-conc", "No. morphs"= "2009-nb_morphs", "Temperature (°C)"="2000-temperature", "Salinity (psu)"="2000-salinity", "NO3 (mL/L)"="2000-no3", "PO4 (mL/L)"="2000-po4", "SiOH4 (mL/L)"="2000-sioh4", "POC (μg/L)"="2000-poc", "PON (μg/L)"="2000-pon",  "Chl a (μg/L)"="2000-chla") #sig 21 years

vars <- c(Zooplankton="2009-conc", `nm morphs`="2009-nb_morphs", `Morph. Richness`="2009-MRic", `Morph. Eveness`="2009-MEve", `Morph. Divergence`="2009-MDiv")
vars <- c(`Morph. Richness`="2009-MRic", `Morph. Eveness`="2009-MEve", `Morph. Divergence`="2009-MDiv")
vars <- c(Ax1="2009-Dim.1", Ax2="2009-Dim.2", Ax3="2009-Dim.3", Ax4="2009-Dim.4")

vars <- c("Ax1-mean"="2009-Dim.1_mean", "Ax2-mean"="2009-Dim.2_mean", "Ax3-mean"="2009-Dim.3_mean", "Ax4-mean"="2009-Dim.4_mean", "Ax1-sd"="2009-Dim.1_sd", "Ax2-sd"="2009-Dim.2_sd", "Ax3-sd"="2009-Dim.3_sd", "Ax4-sd"="2009-Dim.4_sd", "Ax1-q25"="2009-Dim.1_q25", "Ax2-q25"="2009-Dim.2_q25", "Ax3-q25"="2009-Dim.3_q25", "Ax4-q25"="2009-Dim.4_q25", "Ax1-q75"="2009-Dim.1_q75", "Ax2-q75"="2009-Dim.2_q75", "Ax3-q75"="2009-Dim.3_q75", "Ax4-q75"="2009-Dim.4_q75")

dp <- s |>
  mutate(id=str_c(period, var)) |>
  filter(id %in% vars) |>
  mutate(id=names(vars)[match(id, vars)]) |>
  mutate(id=factor(id, levels=names(vars)))


statsp <- stats |>
  mutate(id=str_c(period, var)) |>
  filter(id %in% vars) |>
  mutate(id=names(vars)[match(id, vars)]) |>
  mutate(id=factor(id, levels=names(vars)))

#when drawing 2nd line for recent years (environmental variables)
vars2 <- c(Temperature="2009-temperature", Salinity="2009-salinity", Density="2009-density", Oxygen="2009-o2", NO2="2009-no2", NO3="2009-no3", PO4="2009-po4", SiOH4="2009-sioh4", `Part. org. carbon`="2009-poc", `Part. org. nitrogen`="2009-pon", `Chl. a`="2009-chla")

vars2<- c("Zoo. concentration (ind/L)"="2009-conc", "No. morphs"= "2009-nb_morphs", "Temperature (°C)"="2009-temperature", "Salinity (psu)"="2009-salinity", "Oxygen (mL/L)"="2009-oxygen", "PO4 (mL/L)"="2009-po4", "POC (μg/L)"="2009-poc", "PON (μg/L)"="2009-pon",  "Chl a (μg/L)"="2009-chla") #sig 12 years
vars2<- c("Zoo. concentration (ind/L)"="2009-conc", "No. morphs"= "2009-nb_morphs", "Temperature (°C)"="2009-temperature", "Salinity (psu)"="2009-salinity", "NO3 (mL/L)"="2009-no3", "PO4 (mL/L)"="2009-po4", "SiOH4 (mL/L)"="2009-sioh4", "POC (μg/L)"="2009-poc", "PON (μg/L)"="2009-pon",  "Chl a (μg/L)"="2009-chla") #sig 21 years

statsp2 <- stats |>
  mutate(id=str_c(period, var)) |>
  filter(id %in% vars2) |>
  mutate(id=names(vars2)[match(id, vars2)]) |>
  mutate(id=factor(id, levels=names(vars2)))


ggplot() +
  #facet_grid(id~., scales="free_y") +
  facet_wrap(id~., scales="free_y", ncol=2) +
  geom_path(aes(date, deseason), data=dp, colour="grey20") +
  geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp, signif.gls %in% c("*", "**", "***")), colour="red", size=0.75, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) +
  #geom_segment(aes(x=as.Date("2009-01-01", "%Y-%m-%d"), xend=as.Date("2020-12-17", "%Y-%m-%d"), y=intercept.gls+slope.gls, yend=intercept.gls+slope.gls*12), data=statsp2, colour="blue", size=0.75, alpha=0.7) +
  #annotate("segment", x=as.Date("2009-01-01", "%Y-%m-%d"), xend=as.Date("2019-12-17", "%Y-%m-%d"), y=statsp2$intercept.gls+statsp2$slope.gls*2009, yend=statsp2$intercept.gls+statsp2$slope.gls*2020, colour="pink", size=0.75, alpha=0.7) +
  xlab("Date") +
  ggtitle("Trends in mean coordinates along axes of morpho space - coponly") +
  theme(axis.title.y=element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=11))
ggsave("plots/trends_axes_incl2020_ext_coponly.pdf", width=7, height=14)

