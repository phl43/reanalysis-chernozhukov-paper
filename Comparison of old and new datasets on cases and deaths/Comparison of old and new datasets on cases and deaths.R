JHU <- FALSE
policy_smoothing <- FALSE
L.c <- 14
L.d <- 21

dataset_version <- "Original dataset"

source(paste(rootdir,"cases_and_policies/R/dataprep - updated.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/regprep - updated.R", sep="/"), local = TRUE)

df_old <- df %>%
  dplyr::filter(date >= as.Date("2020-03-07") & date <= as.Date("2020-06-03")) %>%
  dplyr::mutate(dataset_version = "Original") %>%
  dplyr::select(state, date, cases, deaths, dataset_version)

dataset_version <- "Up to date epidemic data but data on policies as in paper"

source(paste(rootdir,"cases_and_policies/R/dataprep - updated.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/regprep - updated.R", sep="/"), local = TRUE)

df_new <- df %>%
  dplyr::filter(date >= as.Date("2020-03-07") & date <= as.Date("2020-06-03")) %>%
  dplyr::mutate(dataset_version = "Revised") %>%
  dplyr::select(state, date, cases, deaths, dataset_version)

df_comparison <- dplyr::bind_rows(df_old, df_new)

dir.create("Comparison of old and new datasets on cases and deaths/Figures")

ggplot(df_comparison %>% dplyr::filter(state == "Wisconsin"), aes(x = date, y = cases, group = dataset_version, color = dataset_version)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Cumulative number of cases in Wisconsin before and after the data were revised") +
  xlab("Date") +
  ylab("Cumulative number of cases") +
  scale_color_discrete(name = "Version of dataset") +
  scale_x_date(
    labels = scales::date_format("%m/%d/%Y"),
    date_breaks = "7 day"
  ) +
  labs(caption = "Source: New York Times - Chart by Philippe Lemoine (@phl43)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  ggsave(
    "Comparison of old and new datasets on cases and deaths/Figures/Cumulative number of cases in Wisconsin before and after the data were revised.png",
    width = 12,
    height = 6
  )
