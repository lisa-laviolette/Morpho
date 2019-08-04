#
# Detect trends in environmental data
#
# (c) 2018 Jean-Olivier Irisson, GNU General Public License v3

library("tidyverse")
library("trend")
library("nlme")
library("broom")

load("7.Rdata")

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
stats <- s %>% group_by(period, var) %>%
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
      glance(mkt) %>% select(p.value) %>% set_names(., paste0, ".mk"),
      glance(m) %>% select(r.squared, p.value, intercept, slope, cor.struct, acf=acf1) %>% set_names(., paste0, ".gls")
    )
  }) %>%
  ungroup() %>%
  mutate(
    acf.gls = abs(acf.gls),
    signif.mk = signif_stars(p.value.mk),
    signif.gls = signif_stars(p.value.gls)
  )

# View(stats)
write_tsv(stats, "plots/trend_stats.tsv")


# TODO check whether copepods are declining more than the reste which could explain the increase of divergence

## Plots ----

vars <- c(Temperature="2000-temperature", NO3="2000-no3", `Fluo (Chl a)`="2000-fluo", `Part. org. Nitrogen`="2000-pon", Zooplankton="2009-conc", `Morph. Divergence`="2009-MDiv")

dp <- s %>%
  mutate(id=str_c(period, var)) %>%
  filter(id %in% vars) %>%
  mutate(id=names(vars)[match(id, vars)]) %>%
  mutate(id=factor(id, levels=names(vars)))


statsp <- stats %>%
  mutate(id=str_c(period, var)) %>%
  filter(id %in% vars) %>%
  mutate(id=names(vars)[match(id, vars)]) %>%
  mutate(id=factor(id, levels=names(vars)))


ggplot() +
  facet_grid(id~., scales="free_y") +
  geom_path(aes(date, deseason), data=dp, colour="grey20") +
  geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=statsp, colour="red", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank())
ggsave("plots/trends.pdf", width=7, height=10)

ggplot() +
  facet_grid(id~., scales="free_y") +
  geom_path(aes(date, deseason), data=filter(dp, id %in% c("Zooplankton", "Morph. Divergence")), colour="grey20") +
  # geom_smooth(aes(date, deseason), data=filter(dp, id %in% c("Zooplankton", "Morph. Divergence")), colour="blue", se=FALSE, span=0.5) +
  geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=filter(statsp, id %in% c("Zooplankton", "Morph. Divergence")), colour="red", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank())
ggsave("plots/trends_zoo.pdf", width=7, height=4)




