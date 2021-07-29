dataset_version <- "Original dataset"
JHU <- FALSE
policy_smoothing <- TRUE
source(paste(rootdir,"cases_and_policies/R/dataprep - updated.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/regprep - updated.R", sep="/"), local = TRUE)
dataset_version <- ifelse(!policy_smoothing, paste0(dataset_version, " - No policy smoothing"), dataset_version)
dataset_version <- ifelse(JHU, paste0(dataset_version, " - JHU"), dataset_version)

policies_timing_original_dataset <- df %>%
  tibble::as_tibble() %>%
  dplyr::filter(date >= as.Date("2020-03-07") & date <= as.Date("2020-06-03")) %>%
  dplyr::rename(
    `Mandating face masks for public-facing employees` = pmaskbus,
    `Closing K-12 schools` = pk12,
    `Stay-at-home order` = pshelter,
    `Closing non-essential businesses` = pindex
  ) %>%
  tidyr::pivot_longer(cols = c("Mandating face masks for public-facing employees", "Closing K-12 schools", "Stay-at-home order", "Closing non-essential businesses"), names_to = "policy_name", values_to = "policy_value") %>%
  dplyr::group_by(policy_name, date) %>%
  dplyr::summarize(
    proportion = mean(policy_value)
  ) %>%
  dplyr::mutate(date = as.Date(date))

dataset_version <- "Up to date epidemic data and data on policies as in latest CUSP database"
JHU <- FALSE
policy_smoothing <- TRUE
source(paste(rootdir,"cases_and_policies/R/dataprep - updated.R", sep="/"), local = TRUE)
source(paste(rootdir,"cases_and_policies/rmd/regprep - updated.R", sep="/"), local = TRUE)
dataset_version <- ifelse(!policy_smoothing, paste0(dataset_version, " - No policy smoothing"), dataset_version)
dataset_version <- ifelse(JHU, paste0(dataset_version, " - JHU"), dataset_version)

policies_timing_uptodate_dataset <- df %>%
  tibble::as_tibble() %>%
  dplyr::filter(date >= as.Date("2020-03-07") & date <= as.Date("2020-06-03")) %>%
  dplyr::rename(
    `Mandating face masks for public-facing employees` = pmaskbus,
    `Closing K-12 schools` = pk12,
    `Stay-at-home order` = pshelter,
    `Closing non-essential businesses` = pindex
  ) %>%
  tidyr::pivot_longer(cols = c("Mandating face masks for public-facing employees", "Closing K-12 schools", "Stay-at-home order", "Closing non-essential businesses"), names_to = "policy_name", values_to = "policy_value") %>%
  dplyr::group_by(policy_name, date) %>%
  dplyr::summarize(
    proportion = mean(policy_value)
  ) %>%
  dplyr::mutate(date = as.Date(date))

policies_timing_original_dataset_plot <- ggplot(policies_timing_original_dataset, aes(x = date, y = proportion)) +
  geom_line(size = 1, color = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  ggtitle(paste0("Proportion of states where each policy was in effect (original dataset)")) +
  xlab("Date") +
  ylab("Proportion") +
  facet_wrap(~ policy_name, ncol = 2) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

policies_timing_uptodate_dataset_plot <- ggplot(policies_timing_uptodate_dataset, aes(x = date, y = proportion)) +
  geom_line(size = 1, color = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  ggtitle(paste0("Proportion of states where each policy was in effect (up to date dataset)")) +
  xlab("Date") +
  ylab("Proportion") +
  facet_wrap(~ policy_name, ncol = 2) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

plot_grid(
  policies_timing_original_dataset_plot,
  policies_timing_uptodate_dataset_plot,
  labels = c("", ""),
  ncol = 1,
  rel_heights = c(1, 1)
) +
  ggsave(
    "Comparison of the timing of policies in different versions of the dataset/Comparison of the timing of policies in different versions of the dataset.png",
    width = 24,
    height = 24,
    limitsize = FALSE
  )
