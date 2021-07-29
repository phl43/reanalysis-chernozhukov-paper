create_results_summary <- function(policy_names, di, bdi, ci) {
  se <- di * NA
  for (i in 1:length(se)) {
    se[i] <- sd(sapply(bdi, function(b) b[i]))
  }
  
  summary <- tibble::tibble()
  
  for (i in 1:length(policy_names)) {
    summary <- summary %>%
      dplyr::bind_rows(
        dplyr::tribble(
          ~policy, ~type, ~effect, ~se, ~lower, ~upper,
          policy_names[i], "Direct", di[i,1], se[i,1], di[i,1] - qnorm(1 - (1 - ci) / 2) * se[i,1], di[i,1] + qnorm(1 - (1 - ci) / 2) * se[i,1],
          policy_names[i], "Indirect", di[i,2], se[i,2], di[i,2] - qnorm(1 - (1 - ci) / 2) * se[i,2], di[i,2] + qnorm(1 - (1 - ci) / 2) * se[i,2],
          policy_names[i], "Total", di[i,3], se[i,3], di[i,3] - qnorm(1 - (1 - ci) / 2) * se[i,3], di[i,3] + qnorm(1 - (1 - ci) / 2) * se[i,3]
        )
      )
  }
  
  summary
}

perform_analysis <- function (df, policy_names, start_date, end_date, yvar, polvars, bvars, infovars, tvars, statevars, lag) {
  xlist <- list("")
  interactions <- list("month", statevars)
  ilist <- list(interactions)
  iv <- list("0")
  
  sdf <- subset(df, df$date >= start_date & df$date <= end_date)
  
  regs <- mainregressions(sdf, yvar, polvars, bvars, infovars, tvars, xlist, ilist, iv, lag)
  
  summary_results <- tibble::tribble(
    ~policy, ~type, ~effect, ~se, ~lower, ~upper
  )
  
  pib <- lapply(regs$pib[[1]][1:4], function(x) x$reg)
  pbiy <- regs$pbiy[[1]]$reg
  piy <- regs$piy[[1]]$reg
  
  S <- 999
  bs <- bootstrap_felm(sdf, c(pib,list(pbiy,piy)), S=S)
  
  diall <- dieff_table(lapply(regs$pib[[1]][1:4], function(x) x$reg$coef[,1]),
                       regs$pbiy[[1]]$reg$coef[,1],
                       regs$piy[[1]]$reg$coef[,1],
                       policies=c(polvars, infovars[[1]]), nsum=length(polvars))
  
  bdiall <- lapply(1:S, function(i)
    dieff_table(lapply(1:4, function(j) bs[[j]][,i]),
                bs[[5]][,i], bs[[6]][,i],
                policies=c(polvars, infovars[[1]]), nsum=length(polvars))
  )
  
  create_results_summary(policy_names, diall, bdiall, 0.95)
}

# seed random number generator
set.seed(21)

# seed for reproducibility with furrr
# (see https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/ and https://davisvaughan.github.io/furrr/reference/furrr_options.html#reproducible-random-number-generation-rng-)
rng_seed <- 127

yvars_list <- list(
  list(
    varnames = "dlogdc",
    description = "weekly growth rate in cases"
  ),
  list(
    varnames = "dlogdd",
    description = "weekly growth rate in deaths"
  )
)

polvars_list <- list(
  list(
    varnames = c("pmaskbus", "pk12", "pshelter", "pindex"),
    description = "mask mandate for public-facing employees, closure of schools, stay-at-home order, business closure policies"
  ),
  list(
    varnames = c("pmaskall", "pk12", "pshelter", "pindex"),
    description = "mask mandate for everyone in public spaces, closure of schools, stay-at-home order, business closure policies"
  ),
  list(
    varnames = c("pmaskall", "pmaskbus", "pk12", "pshelter", "pindex"),
    description = "mask mandate for public-facing employees, mask mandate for everyone in public spaces, closure of schools, stay-at-home order, business closure policies"
  ),
  list(
    varnames = c("pmaskbus","pk12","pshelter","pmovie","prestaurant","pnonessential"),
    description = "mask mandate for public-facing employees, closure of schools, stay-at-home order, closure of movie theaters, closure of restaurants, closure of other non-essential businesses"
  ),
  list(
    varnames = c("pmaskall","pk12","pshelter","pmovie","prestaurant","pnonessential"),
    description = "mask mandate for everyone in public spaces, closure of schools, stay-at-home order, closure of movie theaters, closure of restaurants, closure of other non-essential businesses"
  ),
  list(
    varnames = c("pmaskall", "pmaskbus","pk12","pshelter","pmovie","prestaurant","pnonessential"),
    description = "mask mandate for public-facing employees, mask mandate for everyone in public spaces, closure of schools, stay-at-home order, closure of movie theaters, closure of restaurants, closure of other non-essential businesses"
  )
)

infovars_list <- list(
  list(
    varnames = c("dlogdc", "logdc"),
    description = "number of cases during the past week in state, weekly cases growth in state"
  ),
  list(
    varnames = c("dlogdc", "logdc", "dlogdc.national", "logdc.national"),
    description = "number of cases during the past week in state, weekly cases growth in state, number of cases during the past week nationally, weekly cases growth nationally"
  ),
  list(
    varnames = c("dlogdd", "logdd"),
    description = "number of deaths during the past week in state, weekly deaths growth in state"
  ),
  list(
    varnames = c("dlogdd", "logdd", "dlogdd.national", "logdd.national"),
    description = "number of deaths during the past week in state, weekly deaths growth in state, number of deaths during the past week nationally, weekly deaths growth nationally"
  ),
  list(
    varnames = c("dlogdc", "logdc", "dlogdd", "logdd"),
    description = "number of cases during the past week in state, weekly cases growth in state, number of deaths during the past week in state, weekly deaths growth in state"
  ),
  list(
    varnames = c("dlogdc", "logdc", "dlogdc.national", "logdc.national", "dlogdd", "logdd", "dlogdd.national", "logdd.national"),
    description = "number of cases during the past week in state, weekly cases growth in state, number of cases during the past week nationally, weekly cases growth nationally, number of deaths during the past week in state, weekly deaths growth in state, number of deaths during the past week nationally, weekly deaths growth nationally"
  )
)

bvars <- c(
  "workplaces",
  "retail",
  "grocery",
  "transit"
)

tvars <- "dlogtests"

statevars <- c(
  "log(Population.2018)",
  "log(Square.Miles)",
  "Percent.Unemployed..2018.",
  "Percent.living.under.the.federal.poverty.line..2018.",
  "Percent.at.risk.for.serious.illness.due.to.COVID",
  "party"
)

if (grepl("Up to date epidemic data and data on policies as in latest CUSP database", dataset_version)) {
  period_list <- list(
    c(as.Date("2020-03-07"), as.Date("2020-06-03")),
    c(as.Date("2020-03-07"), as.Date("2020-12-31"))
  )
} else {
  period_list <- list(
    c(as.Date("2020-03-07"), as.Date("2020-06-03"))
  )
}

policy_names_list <- list(
  c(
    "Mandating face masks for public-facing employees",
    "Closing K-12 schools",
    "Stay-at-home order",
    "Closing non-essential businesses",
    "Combined effect of policy"
  ),
  c(
    "Mandating face masks for everyone in public spaces",
    "Closing K-12 schools",
    "Stay-at-home order",
    "Closing non-essential businesses",
    "Combined effect of policy"
  ),
  c(
    "Mandating face masks for public-facing employees",
    "Mandating face masks for everyone in public spaces",
    "Closing K-12 schools",
    "Stay-at-home order",
    "Closing non-essential businesses",
    "Combined effect of policy"
  ),
  c(
    "Mandating face masks for public-facing employees",
    "Closing K-12 schools",
    "Stay-at-home order",
    "Closing movie theaters",
    "Closing restaurants",
    "Closing other non-essential businesses",
    "Combined effect of policy"
  ),
  c(
    "Mandating face masks for everyone in public spaces",
    "Closing K-12 schools",
    "Stay-at-home order",
    "Closing movie theaters",
    "Closing restaurants",
    "Closing other non-essential businesses",
    "Combined effect of policy"
  ),
  c(
    "Mandating face masks for public-facing employees",
    "Mandating face masks for everyone in public spaces",
    "Closing K-12 schools",
    "Stay-at-home order",
    "Closing movie theaters",
    "Closing restaurants",
    "Closing other non-essential businesses",
    "Combined effect of policy"
  )
)

cases_lags <- 7:16
deaths_lags <- 18:28

nb_specifications <- length(polvars_list) * length(infovars_list) * length(period_list) * length(cases_lags) + 
  length(polvars_list) * length(infovars_list) * length(period_list) * length(deaths_lags)

specifications <- list(
  outcome_description = list(),
  policies_description = list(),
  information_description = list(),
  period_description = list(),
  policy_names = list(),
  start_date = list(),
  end_date = list(),
  yvar = list(),
  polvars = list(),
  bvars = rep(list(bvars), nb_specifications),
  infovars = list(),
  tvars = list(),
  statevars = rep(list(statevars), nb_specifications),
  lag = list()
)

for (i in 1:length(yvars_list)) {
  for (j in 1:length(polvars_list)) {
    for (k in 1:length(infovars_list)) {
      for (l in 1:length(period_list)) {
        if (i == 1) {
          lags <- cases_lags
        } else {
          lags <- deaths_lags
        }
        for (m in lags) {
          specifications$outcome_description <- append(specifications$outcome_description, yvars_list[[i]]$description)
          specifications$policies_description <- append(specifications$policies_description, polvars_list[[j]]$description)
          specifications$information_description <- append(specifications$information_description, infovars_list[[k]]$description)
          specifications$period_description <- append(
            specifications$period_description,
            paste0(
              period_list[[l]][1],
              " to ",
              period_list[[l]][2]
            )
          )
          specifications$policy_names <- append(specifications$policy_names, list(policy_names_list[[j]]))
          specifications$start_date <- append(specifications$start_date, list(period_list[[l]][1]))
          specifications$end_date <- append(specifications$end_date, list(period_list[[l]][2]))
          specifications$yvar <- append(specifications$yvar, yvars_list[[i]]$varname)
          specifications$polvars <- append(specifications$polvars, list(polvars_list[[j]]$varnames))
          # infovars needs to be a list of lists even if there is only ever one vector in the list because mainregressions
          # was written to run the regressions using 2 different sets of information variables
          specifications$infovars <- append(specifications$infovars, list(list(infovars_list[[k]]$varnames)))
          specifications$tvars <- append(specifications$tvars, ifelse(yvars_list[[i]]$varname == "dlogdc", list(tvars), list(NULL)))
          specifications$lag <- append(specifications$lag, list(m))
        }
      }
    }
  }
}

# results <- list()
# for (i in 1:nb_specifications) {
#   results[[i]] <- perform_analysis(
#     df,
#     specifications$policy_names[[i]],
#     specifications$start_date[[i]],
#     specifications$end_date[[i]],
#     specifications$yvar[[i]],
#     specifications$polvars[[i]],
#     specifications$bvars[[i]],
#     specifications$infovars[[i]],
#     specifications$tvars[[i]],
#     specifications$statevars[[i]],
#     specifications$lag[[i]]
#   ) %>%
#     tibble::add_column(
#       outcome_description = specifications$outcome_description[[i]],
#       policies_description = specifications$policies_description[[i]],
#       information_description = specifications$information_description[[i]],
#       period_description = specifications$period_description[[i]],
#       lag = specifications$lag[[i]]
#     )
# }
# results <- dplyr::bind_rows(results)

# on my computer, this generates an error when I try to run it for the first time, but it works when I try again
# (note: this is probably because of a bug in furrr, but I don't have time to look into this, so if that still doesn't
# work for you when you try again, just use the code I commented out above and it will work though it will be slower)
progressr::with_progress({
  p <- progressr::progressor(steps = nb_specifications)
  
  results <- furrr::future_pmap_dfr(
    specifications,
    function(
      outcome_description,
      policies_description,
      information_description,
      period_description,
      policy_names,
      start_date,
      end_date,
      yvar,
      polvars,
      bvars,
      infovars,
      tvars,
      statevars,
      lag
    ) {
      p()
      perform_analysis(
        df,
        policy_names,
        start_date,
        end_date,
        yvar,
        polvars,
        bvars,
        infovars,
        tvars,
        statevars,
        lag
      ) %>%
        tibble::add_column(
          outcome_description = outcome_description,
          policies_description = policies_description,
          information_description = information_description,
          period_description = period_description,
          lag = lag
        )
    },
    .options = furrr::furrr_options(seed = rng_seed)
  )
})

saveRDS(
  results,
  paste0(
    "Reanalysis of the results/Results - ",
    dataset_version,
    ".rds"
  )
)

dir.create("Reanalysis of the results/Figures")
figures_directory <- paste("Reanalysis of the results/Figures", dataset_version, sep = "/")
dir.create(figures_directory)
dir.create(
  paste(
    figures_directory,
    "Effects of non-pharmaceutical interventions with lag of 14 days for cases and 21 days for deaths",
    sep = "/"
  )
)

for (i in 1:length(period_list)) {
  period <- paste(
    period_list[[i]][1],
    "to",
    period_list[[i]][2]
  )
  dir.create(
    paste(
      figures_directory,
      "Effects of non-pharmaceutical interventions with lag of 14 days for cases and 21 days for deaths",
      period,
      sep = "/"
    )
  )
  for (j in 1:length(polvars_list)) {
    dir.create(
      paste(
        figures_directory,
        "Effects of non-pharmaceutical interventions with lag of 14 days for cases and 21 days for deaths",
        period,
        paste(polvars_list[[j]]$varnames, collapse = ", "),
        sep = "/"
      )
    )
    for (k in 1:length(infovars_list)) {
      dir.create(
        paste(
          figures_directory,
          "Effects of non-pharmaceutical interventions with lag of 14 days for cases and 21 days for deaths",
          period,
          paste(polvars_list[[j]]$varnames, collapse = ", "),
          paste(infovars_list[[k]]$varnames, collapse = ", "),
          sep = "/"
        )
      )
      
      cases_effects <- results %>%
        dplyr::filter(
          outcome_description == "weekly growth rate in cases" &
            policies_description == polvars_list[[j]]$description &
            information_description == infovars_list[[k]]$description &
            period_description == period &
            lag == 14
        ) %>%
        dplyr::mutate(
          policy = forcats::fct_relevel(policy, policy_names_list[[j]])
        )
      deaths_effects <- results %>%
        dplyr::filter(
          outcome_description == "weekly growth rate in deaths" &
            policies_description == polvars_list[[j]]$description &
            information_description == infovars_list[[k]]$description &
            period_description == period &
            lag == 21
        ) %>%
        dplyr::mutate(
          policy = forcats::fct_relevel(policy, policy_names_list[[j]])
        )
      
      model_description <- paste0(
        "policies = ",
        polvars_list[[j]]$description,
        "; information = ",
        infovars_list[[k]]$description,
        "; period = ",
        period
      )
      model_description <- stringr::str_wrap(model_description, 100)
      
      cases_effects_plot <- ggplot(cases_effects, aes(x = policy, y = effect, group = type, color = type)) +
        geom_pointrange(
          aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.75)
        ) +
        theme_bw() +
        ggtitle(
          "Weekly cases growth (lag of 14 days)"
        ) +
        xlab("Policy") +
        ylab("Effect") +
        scale_color_discrete(name = "Type") +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      deaths_effects_plot <- ggplot(deaths_effects, aes(x = policy, y = effect, group = type, color = type)) +
        geom_pointrange(
          aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.75)
        ) +
        theme_bw() +
        ggtitle(
          "Weekly deaths growth (lag of 21 days)"
        ) +
        xlab("Policy") +
        ylab("Effect") +
        scale_color_discrete(name = "Type") +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      plot_title <- ggdraw() + 
        draw_label(
          "Effects of non-pharmaceutical interventions",
          fontface = "bold"
        )
      
      plot_subtitle <- ggdraw() + 
        draw_label(
          model_description,
          size = 10
        )
      
      path <- paste(
        figures_directory,
        "Effects of non-pharmaceutical interventions with lag of 14 days for cases and 21 days for deaths",
        period,
        paste(polvars_list[[j]]$varnames, collapse = ", "),
        paste(infovars_list[[k]]$varnames, collapse = ", "),
        "Results.png",
        sep = "/"
      )
      
      plot_grid(
        plot_title,
        plot_subtitle,
        cases_effects_plot,
        deaths_effects_plot,
        labels = c("", "", "", ""),
        ncol = 1,
        rel_heights = c(0.5, 1.5, 10, 10)
      ) +
        ggsave(
          path,
          width = length(polvars_list[[j]]$varnames) * 4,
          height = 12 + 3,
          limitsize = FALSE
        )
    }
  }
}

# those are the specifications for which the results are reported in the paper for pmaskbus
(paper_results_pmaskbus <- results %>%
  dplyr::filter(
    period_description == "2020-03-07 to 2020-06-03" &
    !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
    stringr::str_detect(policies_description, "business closure policies") &
    policy == "Mandating face masks for public-facing employees" & type == "Total" & (
      (lag == 14 & outcome_description == "weekly growth rate in cases" & information_description %in% c("number of cases during the past week in state, weekly cases growth in state", "number of cases during the past week in state, weekly cases growth in state, number of cases during the past week nationally, weekly cases growth nationally")) |
      (lag == 21 & outcome_description == "weekly growth rate in deaths" & information_description %in% c("number of deaths during the past week in state, weekly deaths growth in state", "number of deaths during the past week in state, weekly deaths growth in state, number of deaths during the past week nationally, weekly deaths growth nationally"))
    )
  )
)

# those are the specifications for which the results are reported in the paper for pk12
(paper_results_pk12 <- results %>%
  dplyr::filter(
    period_description == "2020-03-07 to 2020-06-03" &
    !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
    stringr::str_detect(policies_description, "business closure policies") &
    policy == "Closing K-12 schools" & type == "Total" & (
      (lag == 14 & outcome_description == "weekly growth rate in cases" & information_description %in% c("number of cases during the past week in state, weekly cases growth in state", "number of cases during the past week in state, weekly cases growth in state, number of cases during the past week nationally, weekly cases growth nationally")) |
      (lag == 21 & outcome_description == "weekly growth rate in deaths" & information_description %in% c("number of deaths during the past week in state, weekly deaths growth in state", "number of deaths during the past week in state, weekly deaths growth in state, number of deaths during the past week nationally, weekly deaths growth nationally"))
    )
  )
)

# those are the specifications for which the results are reported in the paper for pshelter
(paper_results_pshelter <- results %>%
  dplyr::filter(
    period_description == "2020-03-07 to 2020-06-03" &
    !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
    stringr::str_detect(policies_description, "business closure policies") &
    policy == "Stay-at-home order" & type == "Total" & (
      (lag == 14 & outcome_description == "weekly growth rate in cases" & information_description %in% c("number of cases during the past week in state, weekly cases growth in state", "number of cases during the past week in state, weekly cases growth in state, number of cases during the past week nationally, weekly cases growth nationally")) |
      (lag == 21 & outcome_description == "weekly growth rate in deaths" & information_description %in% c("number of deaths during the past week in state, weekly deaths growth in state", "number of deaths during the past week in state, weekly deaths growth in state, number of deaths during the past week nationally, weekly deaths growth nationally"))
    )
  )
)

# those are the specifications for which the results are reported in the paper for pindex
(paper_results_pindex <- results %>%
  dplyr::filter(
    period_description == "2020-03-07 to 2020-06-03" &
    !stringr::str_detect(policies_description, "mask mandate for everyone in public spaces") &
    stringr::str_detect(policies_description, "business closure policies") &
    policy == "Closing non-essential businesses" & type == "Total" & (
      (lag == 14 & outcome_description == "weekly growth rate in cases" & information_description %in% c("number of cases during the past week in state, weekly cases growth in state", "number of cases during the past week in state, weekly cases growth in state, number of cases during the past week nationally, weekly cases growth nationally")) |
      (lag == 21 & outcome_description == "weekly growth rate in deaths" & information_description %in% c("number of deaths during the past week in state, weekly deaths growth in state", "number of deaths during the past week in state, weekly deaths growth in state, number of deaths during the past week nationally, weekly deaths growth nationally"))
    )
  )
)

# those are specifications I consider equally plausible for pmaskbus (excluding those which included pmaskall)
(sensitivity_analysis_pmaskbus <- results %>%
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
    percent_significant = sum(upper < 0) / dplyr::n())
  )

# those are specifications I consider equally plausible for pk12 (excluding those which included pmaskall)
(sensitivity_analysis_pk12 <- results %>%
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
    percent_significant = sum(upper < 0) / dplyr::n())
)

# those are specifications I consider equally plausible for pshelter (excluding those which included pmaskall)
(sensitivity_analysis_pshelter <- results %>%
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
    percent_significant = sum(upper < 0) / dplyr::n())
  )

# those are specifications I consider equally plausible for pindex (excluding those which included pmaskall)
(sensitivity_analysis_pindex <- results %>%
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
    percent_significant = sum(upper < 0) / dplyr::n())
)

(sensitivity_analysis_pmaskall_without_pmaskbus <- results %>%
    dplyr::filter(
      period_description == "2020-03-07 to 2020-12-31" &
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
      percent_significant = sum(upper < 0) / dplyr::n())
)

(sensitivity_analysis_pmaskall_with_pmaskbus <- results %>%
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
      percent_significant = sum(upper < 0) / dplyr::n())
)

(sensitivity_analysis_pmaskbus_with_pmaskall <- results %>%
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
      percent_significant = sum(upper < 0) / dplyr::n())
)

dir.create(
  paste(
    figures_directory,
    "Prevalence of mask mandates",
    sep = "/"
  )
)

for (i in 1:length(period_list)) {
  period <- paste(
    period_list[[i]][1],
    "to",
    period_list[[i]][2]
  )
  
  dir.create(
    paste(
      figures_directory,
      "Prevalence of mask mandates",
      period,
      sep = "/"
    )
  )
  
  mask_mandates_timing <- df %>%
    tibble::as_tibble() %>%
    dplyr::filter(date >= as.Date(period_list[[i]][1]) & date <= as.Date(period_list[[i]][2])) %>%
    dplyr::group_by(date) %>%
    dplyr::summarize(
      pmaskbus = mean(pmaskbus),
      pmaskall = mean(pmaskall)
    )
  
  employees_mask_mandate_plot <- ggplot(mask_mandates_timing, aes(x = date, y = pmaskbus)) +
    geom_line(size = 1, color = "steelblue") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw() +
    ggtitle("Proportion of states in which masks were mandated for employees of public-facing businesses during the first wave in the US") +
    xlab("Date") +
    ylab("Proportion of states") +
    scale_x_date(
      labels = scales::date_format("%m/%d"),
      date_breaks = "7 day"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  universal_mask_mandate_plot <- ggplot(mask_mandates_timing, aes(x = date, y = pmaskall)) +
    geom_line(size = 1, color = "steelblue") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw() +
    ggtitle("Proportion of states in which masks were mandated for everyone in public spaces during the first wave in the US") +
    xlab("Date") +
    ylab("Proportion of states") +
    scale_x_date(
      labels = scales::date_format("%m/%d"),
      date_breaks = "7 day"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  plot_grid(
    employees_mask_mandate_plot,
    universal_mask_mandate_plot,
    labels = c("", ""),
    ncol = 1,
    rel_heights = c(1, 1)
  ) +
    ggsave(
      paste(
        figures_directory,
        "Prevalence of mask mandates",
        period,
        "Proportion of states in which masks were mandated during the first wave in the US.png",
        sep = "/"
      ),
      width = 12,
      height = 12,
      limitsize = FALSE
    )
}
