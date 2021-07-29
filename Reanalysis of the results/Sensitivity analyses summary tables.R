library(tidyverse)
library(gt)

dataset_list <- c(
  "Original dataset",
  "Up to date epidemic data but data on policies as in paper",
  "Up to date epidemic data and data on policies as in latest CUSP database"
  )

JHU_list <- c(FALSE, TRUE)

policy_smoothing_list <- c(TRUE, FALSE)

summary_tables_without_pmaskall <- tibble()
summary_tables_with_pmaskall <- tibble()

for (dataset in dataset_list) {
  for (JHU in JHU_list) {
    for (policy_smoothing in policy_smoothing_list) {
      dataset_version <- dataset
      dataset_version <- ifelse(!policy_smoothing, paste0(dataset_version, " - No policy smoothing"), dataset_version)
      dataset_version <- ifelse(JHU, paste0(dataset_version, " - JHU"), dataset_version)
      
      data_policy_version <- ifelse(
        dataset == "Up to date epidemic data and data on policies as in latest CUSP database",
        "Revised",
        "Same as paper"
      )
      
      if (dataset == "Original dataset" & JHU == FALSE) {
        data_epidemic_version <- "Original (NYT)"
      } else if (dataset == "Original dataset" & JHU == TRUE) {
        data_epidemic_version <- "Original (JHU)"
      } else if (dataset != "Original dataset" & JHU == FALSE) {
        data_epidemic_version <- "Revised (NYT)"
      } else {
        data_epidemic_version <- "Revised (JHU)"
      }
      
      results <- readRDS(paste0("Reanalysis of the results/Results - ", dataset_version, ".rds"))
      
      sensitivity_analysis_pmaskbus <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Mandating face masks for public-facing employees" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
            )
      
      # those are specifications I consider equally plausible for pk12 (excluding those which included pmaskall)
      sensitivity_analysis_pk12 <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Closing K-12 schools" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
          )
      
      # those are specifications I consider equally plausible for pshelter (excluding those which included pmaskall)
      sensitivity_analysis_pshelter <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Stay-at-home order" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
          )
      
      # those are specifications I consider equally plausible for pindex (excluding those which included pmaskall)
      sensitivity_analysis_pindex <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Closing non-essential businesses" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
          )
      
      sensitivity_analysis_pmaskall_without_pmaskbus <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              !stringr::str_detect(policies_description, "mask mandate for public-facing employees") &
              stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Mandating face masks for everyone in public spaces" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
          )
      
      sensitivity_analysis_pmaskall_with_pmaskbus <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              stringr::str_detect(policies_description, "mask mandate for public-facing employees") &
              stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Mandating face masks for everyone in public spaces" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
          )
      
      sensitivity_analysis_pmaskbus_with_pmaskall <- results %>%
          dplyr::filter(
            period_description == "2020-03-07 to 2020-06-03" &
              stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
              policy == "Mandating face masks for public-facing employees" & type == "Total" &
              (
                (lag >= 7 & lag <= 16 & outcome_description == "weekly growth rate in cases") |
                  (lag >= 18 & lag <= 28 & outcome_description == "weekly growth rate in deaths")
              )
          ) %>%
          dplyr::summarize(
            mean_effect = mean(effect),
            smallest = min(effect),
            largest = max(effect),
            percent_significant_negative = sum(upper < 0) / dplyr::n(),
            percent_significant_positive = sum(lower > 0) / dplyr::n()
          )
      
      summary_table_without_pmaskall <- tribble(
        ~`Policy`,
        ~`Data on cases/deaths`,
        ~`Data on policy`,
        ~`Policy smoothing`,
        ~`Mean effect (range of point estimates)`,
        ~`Percent significant (negative)`,
        "Mandating face masks for employees",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pmaskbus$mean_effect, 2), " (", round(sensitivity_analysis_pmaskbus$smallest, 2), ", ", round(sensitivity_analysis_pmaskbus$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pmaskbus$percent_significant_negative * 100), "%"),
        "Closing K-12 schools",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pk12$mean_effect, 2), " (", round(sensitivity_analysis_pk12$smallest, 2), ", ", round(sensitivity_analysis_pk12$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pk12$percent_significant_negative * 100), "%"),
        "Stay-at-home order",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pshelter$mean_effect, 2), " (", round(sensitivity_analysis_pshelter$smallest, 2), ", ", round(sensitivity_analysis_pshelter$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pshelter$percent_significant_negative * 100), "%"),
        "Closing non-essential businesses",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pindex$mean_effect, 2), " (", round(sensitivity_analysis_pindex$smallest, 2), ", ", round(sensitivity_analysis_pindex$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pindex$percent_significant_negative * 100), "%")
      )
      
      summary_table_with_pmaskall <- tribble(
        ~`Policy`,
        ~`Mandating face masks for employees in model`,
        ~`Data on cases/deaths`,
        ~`Data on policy`,
        ~`Policy smoothing`,
        ~`Mean effect (range of point estimates)`,
        ~`Percent significant (negative)`,
        ~`Percent significant (positive)`,
        "Mandating face masks for everyone",
        "No",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pmaskall_without_pmaskbus$mean_effect, 2), " (", round(sensitivity_analysis_pmaskall_without_pmaskbus$smallest, 2), ", ", round(sensitivity_analysis_pmaskall_without_pmaskbus$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pmaskall_without_pmaskbus$percent_significant_negative * 100), "%"),
        paste0(round(sensitivity_analysis_pmaskall_without_pmaskbus$percent_significant_positive * 100), "%"),
        "Mandating face masks for everyone",
        "Yes",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pmaskall_with_pmaskbus$mean_effect, 2), " (", round(sensitivity_analysis_pmaskall_with_pmaskbus$smallest, 2), ", ", round(sensitivity_analysis_pmaskall_with_pmaskbus$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pmaskall_with_pmaskbus$percent_significant_negative * 100), "%"),
        paste0(round(sensitivity_analysis_pmaskall_with_pmaskbus$percent_significant_positive * 100), "%"),
        "Mandating face masks for employees",
        "Yes",
        data_epidemic_version,
        ifelse(dataset == "Up to date epidemic data and data on policies as in latest CUSP database", "Revised", "Same as paper"),
        ifelse(policy_smoothing, "Yes", "No"), 
        paste0(round(sensitivity_analysis_pmaskbus_with_pmaskall$mean_effect, 2), " (", round(sensitivity_analysis_pmaskbus_with_pmaskall$smallest, 2), ", ", round(sensitivity_analysis_pmaskbus_with_pmaskall$largest, 2), ")"),
        paste0(round(sensitivity_analysis_pmaskbus_with_pmaskall$percent_significant_negative * 100), "%"),
        paste0(round(sensitivity_analysis_pmaskbus_with_pmaskall$percent_significant_positive * 100), "%")
      )
      
      summary_tables_without_pmaskall <- bind_rows(summary_tables_without_pmaskall, summary_table_without_pmaskall)
      summary_tables_with_pmaskall <- bind_rows(summary_tables_with_pmaskall, summary_table_with_pmaskall)
    }
  }
}

gt(
  bind_rows(summary_tables_without_pmaskall),
  groupname_col = "Policy"
  ) %>%
  gtsave("Reanalysis of the results/Figures/Summary table of sensitivity analysis without universal mask mandates.png")

gt(
  bind_rows(summary_tables_with_pmaskall),
  groupname_col = "Policy"
  ) %>%
  gtsave("Reanalysis of the results/Figures/Summary table of sensitivity analysis with universal mask mandates.png")
