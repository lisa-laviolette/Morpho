## Explore time series of image features ----

summary(zs)

# compute descriptive statistics per date, for all features
desc <- zs %>%
  select(date, area:perimmajor) %>%
  # reshape the data to be able to quickly compute everything per feature
  gather(key=feature, value=value, -date) %>%
  # compute summary statistic
  group_by(date, feature) %>%
  summarise(
    mean=mean(value,na.rm=TRUE),
    median=median(value,na.rm=TRUE),
    min=min(value,na.rm=TRUE),
    max=max(value,na.rm=TRUE),
    q25=quantile(value, p=0.25, na.rm=TRUE),
    q75=quantile(value, p=0.75, na.rm=TRUE)
  ) %>%
  ungroup()

# plot a time series of the variability of a feature
plot_feature <- function(f) {
  ggplot(filter(desc, feature==f), aes(x=date)) +
    # full range
    # geom_ribbon(aes(ymin=min, ymax=max), alpha=0.2) +
    # quantiles
    geom_ribbon(aes(ymin=q25, ymax=q75), alpha=0.4) +
    geom_line(aes(y=median)) +
    # reference line=median of the medians
    geom_hline(aes(yintercept=median(median)), colour="red") +
    # nice looking plot
    date_x + labs(y=f)
}

# inspect size related features through time
plot_feature("area") + scale_y_log10()
plot_feature("perim.") + scale_y_log10()
plot_feature("major") + scale_y_log10()

# inspect grey level related features through time
plot_feature("mean")
plot_feature("stddev")

# -> no obvious shift through time. Maybe more variability in grey level < 2008 but that may very well be related to the fact that there are fewer objects


