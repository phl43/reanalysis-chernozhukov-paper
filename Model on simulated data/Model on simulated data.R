# in order for the code below to work, you need to install this package I created to simulate epidemics,
# which can be found in the "Model based on simulated data" folder
library(episim)
library(lfe)
library(plm)
library(abind)
library(data.table)
library(reshape2)
library(Rcpp)
#library(fixest)
# library(furrr)
# library(progressr)
library(zoo)
#library(broom)
library(magrittr)
library(ggplot2)
library(cowplot)

# this is the R version of the simulate_epidemic function I distributed in the episim package
# note: it doesn't have exactly the same signature as the C++ function because I can't pass user-defined C++
# functions with Rcpp so I have to assume a particular distribution for the generation time, incubation
# period and infection to death delay (taking as argument their parameters), but I can be more flexible in R
# and take user-defined functions as arguments
# simulate_epidemic <- function(
#   id,
#   simulation_length,
#   population,
#   seed_rate,
#   seed_length,
#   R0,
#   ifr,
#   gt,
#   ip,
#   i2d,
#   prob_detection_cases,
#   prob_detection_deaths,
#   policy_names,
#   policy_matrix,
#   policy_effects,
#   initial_state = NULL
# ) {
#   R <- R0 + policy_matrix %*% policy_effects
#   effective_R <- R
# 
#   # initialize the simulation with the last snapshot in initial state if one was passed
#   if (is.null(initial_state)) {
#     start <- 1
#     seed <- rpois(seed_length, seed_rate)
#     infections <- c(seed, rep(0, simulation_length - seed_length))
#     cases <- rep(0, simulation_length)
#     deaths <- rep(0, simulation_length)
#     recorded_deaths <- rep(0, simulation_length)
#     snapshots <- array(
#       rep(0, simulation_length * simulation_length * 4),
#       dim = c(simulation_length, 4, simulation_length)
#     )
#   } else {
#     start <- dim(initial_state)[3] + 1
#     infections <- initial_state[,1,start - 1]
#     cases <- initial_state[,2,start - 1]
#     deaths <- initial_state[,3,start - 1]
#     recorded_deaths <- initial_state[,4,start - 1]
#     if (start - 1 < seed_length) {
#       infections <- infections + c(rep(0, start - 1), rpois(seed_length - start - 1, seed_rate))
#     }
#     snapshots <- abind(
#       initial_state,
#       array(
#         rep(0, (simulation_length - start + 1) * simulation_length * 4),
#         dim = c(simulation_length, 4, simulation_length - start + 1)
#       ),
#       along = 3
#     )
#   }
# 
#   # simulate infections
#   for (i in start:simulation_length) {
#     effective_R[i] <- R[i] * (1 - sum(infections[1:i]) / population)
#     nb_infections <- infections[i]
#     if (nb_infections > 0 & effective_R[i] > 0) {
#       # see https://wellcomeopenresearch.org/articles/5-67 for using a dispersion paramater of 0.1 and
#       # https://www.nature.com/articles/nature04153 for the logic behind modeling secondary transmission
#       # with a negative binominal distribution
#       nb_secondary_infections <- rnbinom(nb_infections, mu = effective_R[i], size = 0.1)
#       for (j in 1:nb_infections) {
#         if (nb_secondary_infections[j] > 0) {
#           generation_times <- round(gt(nb_secondary_infections[j]))
#           secondary_infections_times <- i + generation_times
# 
#           for (k in 1:nb_secondary_infections[j]) {
#             if (secondary_infections_times[k] <= simulation_length) {
#               infections[secondary_infections_times[k]] <- infections[secondary_infections_times[k]] + 1
#             }
#           }
#         }
#       }
#     }
# 
#     # take a snapshot of the current state of the simulation for infections
#     snapshots[,1,i] <- infections
#   }
# 
#   # simulate cases detection
#   for (i in start:simulation_length) {
#     if (infections[i] > 0) {
#       nb_detected_infections <- rbinom(1, infections[i], prob_detection_cases)
#       if (nb_detected_infections > 0) {
#         incubation_periods <- round(ip(nb_detected_infections))
#         symptoms_onset_times <- i + incubation_periods
#         for (j in 1:nb_detected_infections) {
#           # increment the number of cases of the relevant day
#           if (symptoms_onset_times[j] <= simulation_length) {
#             cases[symptoms_onset_times[j]] <- cases[symptoms_onset_times[j]] + 1
#           }
#         }
#       }
#     }
# 
#     # take a snapshot of the current state of the simulation for cases
#     snapshots[,2,i] <- cases
#   }
# 
#   # simulate deaths
#   for (i in start:simulation_length) {
#     if (infections[i] >= 1) {
#       nb_fatal <- rbinom(1, infections[i], ifr)
#       if (nb_fatal > 0) {
#         death_lags <- round(i2d(nb_fatal))
#         for (j in 1:nb_fatal) {
#           death_time <- i + death_lags[j]
#           if (death_time <= simulation_length) {
#             deaths[death_time] <- deaths[death_time] + 1
#           }
#         }
#       }
#     }
# 
#     # take a snapshot of the current state of the simulation for deaths
#     snapshots[,3,i] <- deaths
#   }
# 
#   # simulate deaths detection
#   for (i in start:simulation_length) {
#     if (deaths[i] > 0) {
#       recorded_deaths[i] <- rbinom(1, deaths[i], prob_detection_deaths)
#     }
# 
#     # take a snapshot of the current state of the simulation for recorded deaths
#     snapshots[,4,i] <- recorded_deaths
#   }
# 
#   simulation <- tibble::tibble(
#     id = id,
#     t = 1:simulation_length,
#     R0 = R,
#     R = effective_R,
#     infections = infections,
#     cases = cases,
#     deaths = deaths,
#     recorded_deaths = recorded_deaths
#   )
# 
#   for (i in 1:length(policy_names)) {
#     simulation <- simulation %>%
#       tibble::add_column(
#         "{policy_names[i]}" := policy_matrix[,i]
#       )
#   }
# 
#   list(
#     end_result = simulation,
#     snapshots = snapshots
#   )
# }

create_policy_matrix <- function(policies, simulation_length, nb_simulations) {
  timing_policies <- list()
  policy_matrix <- rep(list(matrix(0, nrow = simulation_length, ncol = length(policies))), nb_simulations)
  
  for (i in 1:nb_simulations) {
    # draw start and end dates for each policy
    timing_policies[[i]] <- tibble::tribble(~policy, ~start, ~end)
    for (j in 1:length(policies)) {
      if (rbinom(1, 1, policies[[j]]$prob)) {
        start <- round(policies[[j]]$start(1))
        start <- ifelse(start > simulation_length, simulation_length, start)
        end <- start + round(policies[[j]]$duration(1))
        
        timing_policies[[i]] <- timing_policies[[i]] %>%
          dplyr::add_row(
            policy = policies[[j]]$long_name,
            start = start,
            end = end
          )
      } else {
        timing_policies[[i]] <- timing_policies[[i]] %>%
          dplyr::add_row(
            policy = policies[[j]]$long_name,
            start = 0,
            end = 0
          )
      }
    }
    
    # fill the matrix based on the previously drawn dates
    for (j in 1:simulation_length) {
      for (k in 1:length(policies)) {
        if (j >= timing_policies[[i]]$start[k] & j <= timing_policies[[i]]$end[k]) {
          policy_matrix[[i]][j, k] = 1
        } else {
          policy_matrix[[i]][j, k] = 0
        }
      }
    }
    
    # for (j in 1:length(policies)) {
    #   policy_matrix[[i]][, j] <- rollmean(policy_matrix[[i]][,j], difflag, fill = rep("extend", difflag - 1), align = "right")
    # }
  }
  
  list(
    timing_policies = timing_policies,
    policy_matrix = policy_matrix
  )
}

# since as explained above the signature of the C++ version of simulate_epidemic is not exactly the same as that
# of the R version, I also need to modify simulate_data if I want to use the C++ version of simulate_epidemic, so
# I commented out the code of the function that uses the R version of simulate_epidemic
# simulate_data <- function (
#   nb_simulations,
#   simulation_length,
#   population,
#   seed_rate,
#   seed_length,
#   R0,
#   gt,
#   ip,
#   ifr,
#   i2d,
#   prob_detection_cases,
#   prob_detection_deaths,
#   policies,
#   initial_state = NULL
# ) {
#   initial_conditions <- tibble::tibble(
#     id = 1:nb_simulations,
#     R0 = R0(nb_simulations),
#     ifr = ifr(nb_simulations)
#   )
# 
#   args <- list(
#     id = as.list(1:nb_simulations),
#     R0 = as.list(initial_conditions$R0),
#     ifr = as.list(initial_conditions$ifr)
#   )
# 
#   policy_names = unlist(purrr::map(policies, ~ .x$short_name))
#   policy_effects = unlist(purrr::map(policies, ~ .x$effect))
# 
#   timing_matrix_policies <- create_policy_matrix(policies, simulation_length, nb_simulations)
#   timing_policies <- timing_matrix_policies$timing_policies
#   args[["policy_matrix"]] <- timing_matrix_policies$policy_matrix
# 
#   progressr::with_progress({
#     p <- progressr::progressor(steps = nb_simulations)
# 
#     simulations <- furrr::future_pmap(
#       args,
#       function(
#         id,
#         R0,
#         ifr,
#         policy_matrix
#       ) {
#         p()
#         simulate_epidemic(
#           id,
#           simulation_length,
#           population,
#           seed_rate,
#           seed_length,
#           R0,
#           ifr,
#           gt,
#           ip,
#           i2d,
#           prob_detection_cases,
#           prob_detection_deaths,
#           policy_names,
#           policy_matrix,
#           policy_effects,
#           initial_state
#         )
#       },
#       .options = furrr::furrr_options(seed = rng_seed)
#     )
#   })
# 
#   end_results <- purrr::map(simulations, ~ .x$end_result) %>%
#     dplyr::bind_rows()
# 
#   snapshots <- purrr::map(simulations, ~ .x$snapshots)
# 
#   list(
#     end_results = end_results,
#     snapshots = snapshots,
#     timing_policies = timing_policies,
#     initial_conditions = initial_conditions
#   )
# }

simulate_data <- function (
  nb_simulations,
  simulation_length,
  population,
  seed_rate,
  seed_length,
  R0,
  ifr,
  shape_gt,
  rate_gt,
  shape_ip,
  rate_ip,
  meanlog_i2d,
  sdlog_i2d,
  prob_detection_cases,
  prob_detection_deaths,
  policies,
  initial_state = NULL
) {
  initial_conditions <- tibble::tibble(
    id = 1:nb_simulations,
    R0 = R0(nb_simulations),
    ifr = ifr(nb_simulations)
  )

  args <- list(
    id = as.list(1:nb_simulations),
    R0 = as.list(initial_conditions$R0),
    ifr = as.list(initial_conditions$ifr)
  )

  policy_names = unlist(purrr::map(policies, ~ .x$short_name))
  policy_effects = unlist(purrr::map(policies, ~ .x$effect))

  timing_matrix_policies <- create_policy_matrix(policies, simulation_length, nb_simulations)
  timing_policies <- timing_matrix_policies$timing_policies
  args[["policy_matrix"]] <- timing_matrix_policies$policy_matrix

  progressr::with_progress({
    p <- progressr::progressor(steps = nb_simulations)

    simulations <- furrr::future_pmap(
      args,
      function(
        id,
        R0,
        ifr,
        policy_matrix
      ) {
        p()
        simulate_epidemic(
          id,
          simulation_length,
          population,
          seed_rate,
          seed_length,
          R0,
          ifr,
          shape_gt,
          rate_gt,
          shape_ip,
          rate_ip,
          meanlog_i2d,
          sdlog_i2d,
          prob_detection_cases,
          prob_detection_deaths,
          policy_names,
          policy_matrix,
          policy_effects,
          initial_state
        )
      },
      .options = furrr::furrr_options(seed = rng_seed)
    )
  })

  # .x$end_result is a list because that's what the C++ verson of simulate_epidemic returns, but it turns out that
  # bind_rows automatically converts it to a tibble, so I don't need to explicitly do the conversion first
  end_results <- purrr::map(simulations, ~ .x$end_result) %>%
    dplyr::bind_rows()

  snapshots <- purrr::map(simulations, ~ .x$snapshots)

  list(
    end_results = end_results,
    snapshots = snapshots,
    timing_policies = timing_policies,
    initial_conditions = initial_conditions
  )
}

estimate_model <- function (data, yvar, pvars, pnames, lag, id_fixed_effect = FALSE) {
  rhs <- paste0(
    paste(
      sprintf("lag(%s, %d)", pvars, lag),
      collapse = " + "
    ),
    ifelse(id_fixed_effect, " | id | 0 | id", " | 0 | 0 | id")
  )
  
  formula <- as.formula(
    paste(yvar, rhs, sep = " ~ ")
  )
  
  fit <- felm(formula, data, keepX = TRUE)
  
  # rhs <- paste(
  #   sprintf("l(%s, %d)", pvars, lag),
  #   collapse = " + "
  # )
  # 
  # formula <- as.formula(
  #   paste(yvar, rhs, sep = " ~ ")
  # )
  # 
  # fit <- fixest::femlm(
  #   formula,
  #   data = data_for_regression,
  #   panel.id = ~ id + t,
  #   cluster = c("id"),
  #   family = "gaussian"
  # )
  
  fit_summary <- broom::tidy(fit, conf.int = TRUE) %>%
    dplyr::filter(
      term != "(Intercept)"
    ) %>%
    dplyr::mutate(
      term = pnames,
      lag = lag
    ) %>%
    dplyr::rename(
      policy = term,
      lower = conf.low,
      upper = conf.high
    )
  
  list(
    fit = fit,
    fit_summary = fit_summary
  )
}

get_degenerate_random_variate_generator <- function (constant) {
  x <- constant
  
  function(n) rep(x, n)
}

create_simulation_function <- function(fit, counterfactual_data, dvar, var, lag, difflag = 7) {
  # use update.felm.formula to get a features matrix based on the counterfactual_data and get the indices of the
  # observations used to fit the model (i. e. for each epidemic, every observation except the lag first observations)
  counterfactual_data$idxcfs <- 1:nrow(counterfactual_data)
  mid <- felm(update.felm.formula(fit$formula, idxcfs ~ .), data = counterfactual_data, keepX = TRUE)
  idx <- mid$response
  
  # create a vector that, for each observation that was used to fit the model, contains the id of the epidemic to
  # which that observation belongs
  ids <- counterfactual_data[idx, "id"]
  
  function(id) {
    # randomly draw coefficients from the asymptotic distribution of the estimators
    V <- vcov(fit)
    e <- t(rmvnorm(1, mean = rep(0, nrow(V)), sigma = V))
    beta <- fit$coefficients + e
    
    # compute the residuals with actual policies based on the coefficients I have just drawn and keep
    # only those for the epidemic whose id was passed as argument to the function
    residuals <- (fit$response - fit$X %*% beta)[ids == id]
    
    # compute the effect of the counterfactual policies based on the coefficients I have just drawn and
    # keep only those for the epidemic whose id was passed as argument to the function
    counterfactual_policy_effect <- (mid$X %*% beta)[ids == id]
    
    # get the index of the first observation of the dependent variable that was used to fit the model for the
    # epidemic whose id was passed as argument to the function (which should be the (lag + 1)th observation)
    t0 <- counterfactual_data$t[idx[counterfactual_data[idx,"id"] == id]][1]
    index_t0 <- which(as.vector(counterfactual_data$t) == t0 & as.vector(counterfactual_data[,"id"]) == id)
    
    # pick the largest of lag and difflag as the length of the initial segment of values used to recursively
    # compute the counterfactual outcomes of interest logdc/logdd
    maxlags <- max(lag, difflag)
    
    # get the first lag observations of logdc/logdd and dlogdc/dlogdd in reverse order for the epidemic
    # whose id was passed as argument to the function
    # Y0 <- as.vector(sapply(1:maxlags, function(k) lag(counterfactual_data[,var], k)[index_t0]))
    # dY0 <- as.vector(sapply(1:maxlags, function(k) lag(counterfactual_data[,dvar], k)[index_t0]))
    Y0 <- rev(counterfactual_data[(index_t0 - lag):(index_t0 - 1),var])
    dY0 <- rev(counterfactual_data[(index_t0 - lag):(index_t0 - 1),dvar])
    
    # filter with zero everywhere except at difflag where it's set to 1 because, by definition of dlogdc
    # and given the model I'm using, logdc(t + lag) = X %*% beta + residual + logdc(t + lag - difflag)
    # (cf. appendix B.4 of the paper and blog post for more details) and same thing for logdd and dlogdd
    filter <- rep(0, maxlags)
    filter[difflag] <- 1
    
    # iteratively compute the counterfactual time series of dlogdc/dlogdd logdc/logdd
    dY1 <- as.vector(counterfactual_policy_effect + residuals)
    Y1 <- as.vector(stats::filter(dY1, filter, method = "recursive", init = Y0))
    
    tibble::tibble(
      t = 1:(length(Y0) + length(Y1)),
      "{var}" := c(rev(Y0), Y1),
      "{dvar}" := c(rev(dY0), dY1)
    )
  }
}

perform_counterfactual_simulation <- function(actual_data, counterfactual_data, yvar, pvars, pnames, lag, difflag = 7, nb_simulations) {
  fit <- estimate_model(actual_data, yvar, pvars, pnames, lag)$fit
  ids <- unique(actual_data$id)
  
  simulate_counterfactual <- create_simulation_function(
    fit,
    counterfactual_data,
    yvar,
    stringr::str_sub(yvar, 2, stringr::str_length(yvar)),
    lag,
    difflag
  )
  
  progressr::with_progress({
    p <- progressr::progressor(steps = length(ids))
    
    counterfactual_simulations <- list()
    for (id in ids) {
      counterfactual_simulations[[id]] <- furrr::future_map_dfr(
        1:nb_simulations,
        function(
          id_simulation
        ) {
          simulate_counterfactual(id) %>%
            dplyr::mutate(
              id = id,
              id_simulation = id_simulation
            )
        },
        .options = furrr::furrr_options(seed = rng_seed)
      )
      
      p()
    }
  })
  
  dplyr::bind_rows(counterfactual_simulations)
}

# seed random number generator
set.seed(21)

# seed for reproducibility with furrr
# (see https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/ and https://davisvaughan.github.io/furrr/reference/furrr_options.html#reproducible-random-number-generation-rng-)
rng_seed <- 127

# length of the difference operator used, in case I want to see how much better the model performs when I set it to 1
difflag <- 7

simulation_length <- 120
population <- 10e6
seed_rate <- 0.5
seed_length <- 30
prob_detection_cases <- 0.1
prob_detection_deaths <- 0.9

# mean_R0 <- 2.5
# sd_R0 <- 0.1
# R0 <- purrr::partial(rnorm, mean = mean_R0, sd = sd_R0)
R0 <- function (n) { rep(2.5, n) }

# I'm using the meta-analytic estimate of the generation time distribution in https://www.medrxiv.org/content/10.1101/2020.11.17.20231548v2
mean_gt <- 4.8
sd_gt <- 1.7
gt <- purrr::partial(rgamma, shape = (mean_gt / sd_gt)^2, rate = mean_gt / sd_gt^2)
#gt <- function (n) { rep(5, n) }

# I'm using the Gamma distributed estimate in https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.5.2000062
mean_ip <- 6.5
sd_ip <- 2.6
ip <- purrr::partial(rgamma, shape = (mean_ip / sd_ip)^2, rate = mean_ip / sd_ip^2)
#ip <- function (n) { rep(7, n)}

# mean_ifr <- 0.005
# sd_ifr <- 0.001
# ifr <- purrr::partial(rnorm, mean = mean_ifr, sd = sd_ifr)
ifr <- function (n) { rep(0.005, n) }

# I'm using the estimate of the infection to death distribution in https://onlinelibrary.wiley.com/doi/10.1111/biom.13462
# (note: the paper doesn't give the value of the parameters, but you can find them by running the code)
mean_i2d <- 25.35593
sd_i2d <- 11.08799
i2d <- purrr::partial(rlnorm, meanlog = 2 * log(mean_i2d) - log(sd_i2d^2 + mean_i2d^2) / 2, sdlog = - 2 * log(mean_i2d) + log(sd_i2d^2 + mean_i2d^2))

# mean_masks_start <- 60
# sd_masks_start <- 25
# mean_masks_duration <- 365
# sd_masks_duration <- 50
# masks_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_masks_start) - log(sd_masks_start^2 + mean_masks_start^2) / 2, sdlog = - 2 * log(mean_masks_start) + log(sd_masks_start^2 + mean_masks_start^2))
# masks_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_masks_duration) - log(sd_masks_duration^2 + mean_masks_duration^2) / 2, sdlog = - 2 * log(mean_masks_duration) + log(sd_masks_duration^2 + mean_masks_duration^2))
# 
# mean_schools_start <- 45
# sd_schools_start <- 25
# mean_schools_duration <- 200
# sd_schools_duration <- 15
# schools_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_schools_start) - log(sd_schools_start^2 + mean_schools_start^2) / 2, sdlog = - 2 * log(mean_schools_start) + log(sd_schools_start^2 + mean_schools_start^2))
# schools_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_schools_duration) - log(sd_schools_duration^2 + mean_schools_duration^2) / 2, sdlog = - 2 * log(mean_schools_duration) + log(sd_schools_duration^2 + mean_schools_duration^2))
# 
# mean_shelter_start <- 55
# sd_shelter_start <- 25
# mean_shelter_duration <- 50
# sd_shelter_duration <- 15
# shelter_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_shelter_start) - log(sd_shelter_start^2 + mean_shelter_start^2) / 2, sdlog = - 2 * log(mean_shelter_start) + log(sd_shelter_start^2 + mean_shelter_start^2))
# shelter_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_shelter_duration) - log(sd_shelter_duration^2 + mean_shelter_duration^2) / 2, sdlog = - 2 * log(mean_shelter_duration) + log(sd_shelter_duration^2 + mean_shelter_duration^2))
# 
# mean_business_start <- 50
# sd_business_start <- 25
# mean_business_duration <- 60
# sd_business_duration <- 25
# business_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_business_start) - log(sd_business_start^2 + mean_business_start^2) / 2, sdlog = - 2 * log(mean_business_start) + log(sd_business_start^2 + mean_business_start^2))
# business_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_business_duration) - log(sd_business_duration^2 + mean_business_duration^2) / 2, sdlog = - 2 * log(mean_business_duration) + log(sd_business_duration^2 + mean_business_duration^2))

mean_masks_start <- 45
sd_masks_start <- 30
mean_masks_duration <- 60
sd_masks_duration <- 30
masks_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_masks_start) - log(sd_masks_start^2 + mean_masks_start^2) / 2, sdlog = - 2 * log(mean_masks_start) + log(sd_masks_start^2 + mean_masks_start^2))
masks_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_masks_duration) - log(sd_masks_duration^2 + mean_masks_duration^2) / 2, sdlog = - 2 * log(mean_masks_duration) + log(sd_masks_duration^2 + mean_masks_duration^2))

mean_schools_start <- 45
sd_schools_start <- 30
mean_schools_duration <- 60
sd_schools_duration <- 30
schools_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_schools_start) - log(sd_schools_start^2 + mean_schools_start^2) / 2, sdlog = - 2 * log(mean_schools_start) + log(sd_schools_start^2 + mean_schools_start^2))
schools_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_schools_duration) - log(sd_schools_duration^2 + mean_schools_duration^2) / 2, sdlog = - 2 * log(mean_schools_duration) + log(sd_schools_duration^2 + mean_schools_duration^2))

mean_shelter_start <- 45
sd_shelter_start <- 30
mean_shelter_duration <- 60
sd_shelter_duration <- 30
shelter_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_shelter_start) - log(sd_shelter_start^2 + mean_shelter_start^2) / 2, sdlog = - 2 * log(mean_shelter_start) + log(sd_shelter_start^2 + mean_shelter_start^2))
shelter_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_shelter_duration) - log(sd_shelter_duration^2 + mean_shelter_duration^2) / 2, sdlog = - 2 * log(mean_shelter_duration) + log(sd_shelter_duration^2 + mean_shelter_duration^2))

mean_business_start <- 45
sd_business_start <- 30
mean_business_duration <- 60
sd_business_duration <- 30
business_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_business_start) - log(sd_business_start^2 + mean_business_start^2) / 2, sdlog = - 2 * log(mean_business_start) + log(sd_business_start^2 + mean_business_start^2))
business_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_business_duration) - log(sd_business_duration^2 + mean_business_duration^2) / 2, sdlog = - 2 * log(mean_business_duration) + log(sd_business_duration^2 + mean_business_duration^2))

policies <- list(
  list(
    long_name = "Mandating masks",
    short_name = "masks",
    prob = 1,
    start = masks_start,
    duration = masks_duration,
    effect = - 1.6 / 3
  ),
  list(
    long_name = "Schools closure",
    short_name = "schools",
    prob = 1,
    start = schools_start,
    duration = schools_duration,
    effect = - 1.6 / 3
  ),
  list(
    long_name = "Stay-at-home order",
    short_name = "stay",
    prob = 1,
    start = shelter_start,
    duration = shelter_duration,
    effect = - 1.6 / 3
  ),
  list(
    long_name = "Non-essential businesses closure",
    short_name = "business",
    prob = 1,
    start = business_start,
    duration = business_duration,
    effect = - 1.6 / 3
  )
)

# this is the code for calling the version of simulate_data that calls the R version of simulate_epidemic,
# so I commented it out because I'm using the C++ version of simulate_epidemic, for which I had to write
# a slightly different version of simulate_data, as explained above
# simulation_results <- simulate_data(
#   nb_simulations = 500,
#   simulation_length = 120,
#   population = population,
#   seed_rate = seed_rate,
#   seed_length = seed_length,
#   R0 = R0,
#   ifr = ifr,
#   gt = gt,
#   ip = ip,
#   i2d = i2d,
#   prob_detection_cases = prob_detection_cases,
#   prob_detection_deaths = prob_detection_deaths,
#   policies = policies
# )

simulation_results <- simulate_data(
  nb_simulations = 500,
  simulation_length = 120,
  population = population,
  seed_rate = seed_rate,
  seed_length = seed_length,
  R0 = R0,
  ifr = ifr,
  shape_gt = (mean_gt / sd_gt)^2,
  rate_gt = mean_gt / sd_gt^2,
  shape_ip = (mean_ip / sd_ip)^2,
  rate_ip = mean_ip / sd_ip^2,
  meanlog_i2d = 2 * log(mean_i2d) - log(sd_i2d^2 + mean_i2d^2) / 2,
  sdlog_i2d = - 2 * log(mean_i2d) + log(sd_i2d^2 + mean_i2d^2),
  prob_detection_cases = prob_detection_cases,
  prob_detection_deaths = prob_detection_deaths,
  policies = policies
)

# check for which simulated epidemics S(t)/N is approximately equal to 1 because otherwise (9) in the paper doesn't hold
attack_rates <- simulation_results$end_results %>%
  dplyr::group_by(id) %>%
  dplyr::summarize(attack_rate = sum(infections) / population) %>%
  dplyr::filter(attack_rate >= 0.005 & attack_rate <= 0.1)

# randomly choose 50 epidemics among those for which S(t)/N is approximately equal to 1 that will be used to
# fit the model, because otherwise (9) in the paper doesn't hold and we expect the results to be biased
sample <- sample(attack_rates$id, 50)

data_for_regression <- simulation_results$end_results %>%
  dplyr::filter(id %in% sample) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(
    dc = rollsum(cases, difflag, fill = rep(0, difflag - 1), align = "right"),
    logdc = ifelse(dc > 0, log(dc), -1),
    dlogdc = logdc - dplyr::lag(logdc, difflag, default = -1),
    dd = rollsum(recorded_deaths, difflag, fill = rep(0, difflag - 1), align = "right"),
    logdd = ifelse(dd > 0, log(dd), -1),
    dlogdd = logdd - dplyr::lag(logdd, difflag, default = -1),
    masks = rollmean(masks, difflag, fill = rep("extend", difflag - 1), align = "right"),
    schools = rollmean(schools, difflag, fill = rep("extend", difflag - 1), align = "right"),
    stay = rollmean(stay, difflag, fill = rep("extend", difflag - 1), align = "right"),
    business = rollmean(business, difflag, fill = rep("extend", difflag - 1), align = "right")
  ) %>%
  dplyr::ungroup() %>%
  pdata.frame(index = c("id", "t"))

policies_timing_actual <- df %>%
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

policies_timing_simulation <- data_for_regression %>%
  dplyr::rename(
    `Mandating face masks` = masks,
    `Closing K-12 schools` = schools,
    `Stay-at-home order` = stay,
    `Closing non-essential businesses` = business
  ) %>%
  tidyr::pivot_longer(cols = c("Mandating face masks", "Closing K-12 schools", "Stay-at-home order", "Closing non-essential businesses"), names_to = "policy_name", values_to = "policy_value") %>%
  dplyr::group_by(policy_name, t) %>%
  dplyr::summarize(
    proportion = mean(policy_value)
  ) %>%
  dplyr::mutate(t = as.integer(t))

policies_timing_actual_plot <- ggplot(policies_timing_actual, aes(x = date, y = proportion)) +
  geom_line(size = 1, color = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  ggtitle(paste0("Proportion of states in the real world where each policy was in effect")) +
  xlab("Date") +
  ylab("Proportion") +
  facet_wrap(~ policy_name, ncol = 2) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

policies_timing_simulation_plot <- ggplot(policies_timing_simulation, aes(x = t, y = proportion)) +
  geom_line(size = 1, color = "steelblue") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  ggtitle(paste0("Proportion of states in the simulation where each policy was in effect")) +
  xlab("Time") +
  ylab("Proportion") +
  facet_wrap(~ policy_name, ncol = 2) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

plot_grid(
  policies_timing_actual_plot,
  policies_timing_simulation_plot,
  labels = c("", ""),
  ncol = 1,
  rel_heights = c(1, 1)
) +
  ggsave(
    "Model on simulated data/Figures/Proportion of states where each policy was in effect (reality and simulation).png",
    width = 24,
    height = 24,
    limitsize = FALSE
  )

data_for_correlation_table_actual <- df[, c("date", "pmaskbus", "pk12", "pshelter", "pindex")] %>%
  tibble::as_tibble() %>%
  dplyr::filter(date >= as.Date("2020-03-07") & date <= as.Date("2020-06-03")) %>%
  dplyr::select(-date) %>%
  dplyr::rename(
    `Mandating face masks for public-facing employees` = pmaskbus,
    `Closing K-12 schools` = pk12,
    `Stay-at-home order` = pshelter,
    `Closing non-essential businesses` = pindex
  )

data_for_correlation_table_simulation <- data_for_regression[, c("masks", "schools", "stay", "business")] %>%
  dplyr::rename(
    `Mandating face masks` = masks,
    `Closing K-12 schools` = schools,
    `Stay-at-home order` = stay,
    `Closing non-essential businesses` = business
  )

modelsummary::datasummary_correlation(
  data_for_correlation_table_actual,
  title = "Correlation table for policies in the real world",
  output = "gt"
) %>%
  gt::gtsave(
    "Model on simulated data/Figures/Correlation table for policies in the real world.png"
  )

modelsummary::datasummary_correlation(
  data_for_correlation_table_simulation,
  title = "Correlation table for policies in the simulation",
  output = "gt"
  ) %>%
  gt::gtsave(
    "Model on simulated data/Figures/Correlation table for policies in the simulation.png"
  )

model_fit_cases <- estimate_model(
  data_for_regression,
  "dlogdc",
  unlist(purrr::map(policies, ~ .x$short_name)),
  unlist(purrr::map(policies, ~ .x$long_name)),
  round(mean_gt + mean_ip)
)$fit_summary

model_fit_deaths <- estimate_model(
  data_for_regression,
  "dlogdd",
  unlist(purrr::map(policies, ~ .x$short_name)),
  unlist(purrr::map(policies, ~ .x$long_name)),
  round(mean_i2d)
)$fit_summary

cases_effects_plot <- ggplot(model_fit_cases, aes(x = policy, y = estimate)) +
  geom_pointrange(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.75),
    color = "steelblue"
  ) +
  theme_bw() +
  ggtitle(
    paste0("Weekly cases growth (lag of ", round(mean_gt + mean_ip), " days)")
  ) +
  xlab("Policy") +
  ylab("Effect") +
  scale_color_discrete(name = "Type") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

deaths_effects_plot <- ggplot(model_fit_deaths, aes(x = policy, y = estimate)) +
  geom_pointrange(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.75),
    color = "steelblue"
  ) +
  theme_bw() +
  ggtitle(
    paste0("Weekly cases growth (lag of ", round(mean_i2d), " days)")
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
    "Effects of non-pharmaceutical interventions according to a model similar to that used in Chernozhukov et al. (2021)\nbut fitted on simulated data",
    fontface = "bold"
  )

dir.create("Model on simulated data/Figures")

plot_grid(
  plot_title,
  cases_effects_plot,
  deaths_effects_plot,
  labels = c("", "", ""),
  ncol = 1,
  rel_heights = c(1.5, 10, 10)
) +
  ggsave(
    "Model on simulated data/Figures/Model fit on simulated data.png",
    width = 18,
    height = 12,
    limitsize = FALSE
  )

nb_counterfactual_simulations <- 10

# compute "real" counterfactual, obtained by running the model that generated the data
# in the first place, but removing each policy in turn
# real_counterfactual_simulation_results <- tibble::tibble()
# for (i in 1:length(sample)) {
#   timing_policies <- simulation_results$timing_policies[[sample[i]]]
#   for (j in 1:length(policies)) {
#     counterfactual_policies <- policies
# 
#     # make sure the policies start at the same dates that were used to generate the data in the first place
#     for (k in 1:length(counterfactual_policies)) {
#       counterfactual_policies[[k]]$start <- get_degenerate_random_variate_generator(timing_policies$start[k])
#       counterfactual_policies[[k]]$duration <- get_degenerate_random_variate_generator(timing_policies$end[k] - timing_policies$start[k])
#     }
# 
#     # remove policy i
#     counterfactual_policies[[j]]$prob <- 0
# 
#     # get the initial state of the simulation of the actual data up to the day before policy j came
#     # into effect so that the counterfactuals are identical to the actual data up to that point
#     initial_state <- simulation_results$snapshots[[sample[i]]][,,1:(timing_policies$start[j] - 1)]
# 
#     # compute the counterfactuals
#     real_counterfactual_simulation_results <- simulate_data(
#       nb_simulations = nb_counterfactual_simulations,
#       simulation_length = simulation_length,
#       population = population,
#       seed_rate = seed_rate,
#       seed_length = seed_length,
#       R0 = simulation_results$initial_conditions$R0[sample[i]],
#       gt = gt,
#       ip = ip,
#       ifr = simulation_results$initial_conditions$ifr[sample[i]],
#       i2d = i2d,
#       prob_detection_cases = prob_detection_cases,
#       prob_detection_deaths = prob_detection_deaths,
#       policies = counterfactual_policies,
#       initial_state = initial_state
#     )$end_results %>%
#       dplyr::rename(
#         id_simulation = id
#       ) %>%
#       dplyr::mutate(
#         id = sample[i],
#         removed_policy = policies[[j]]$long_name
#       ) %>%
#       dplyr::group_by(id_simulation) %>%
#       dplyr::mutate(
#         cumulative_cases = cumsum(cases)
#       ) %>%
#       dplyr::bind_rows(real_counterfactual_simulation_results)
# 
#     cat(sprintf("Simulation of counterfactuals for unit %d without policy %s done\r", i, policies[[j]]$long_name))
#     Sys.sleep(1)
#   }
# }

real_counterfactual_simulation_results <- tibble::tibble()
for (i in 1:length(sample)) {
  timing_policies <- simulation_results$timing_policies[[sample[i]]]
  for (j in 1:length(policies)) {
    counterfactual_policies <- policies
    
    # make sure the policies start at the same dates that were used to generate the data in the first place
    for (k in 1:length(counterfactual_policies)) {
      counterfactual_policies[[k]]$start <- get_degenerate_random_variate_generator(timing_policies$start[k])
      counterfactual_policies[[k]]$duration <- get_degenerate_random_variate_generator(timing_policies$end[k] - timing_policies$start[k])
    }
    
    # remove policy i
    counterfactual_policies[[j]]$prob <- 0
    
    # get the initial state of the simulation of the actual data up to the day before policy j came
    # into effect so that the counterfactuals are identical to the actual data up to that point
    initial_state <- simulation_results$snapshots[[sample[i]]][,,1:(timing_policies$start[j] - 1)]
    
    # compute the counterfactuals
    counterfactual_simulation_results <- simulate_data(
      nb_simulations = nb_counterfactual_simulations,
      simulation_length = simulation_length,
      population = population,
      seed_rate = seed_rate,
      seed_length = seed_length,
      R0 = get_degenerate_random_variate_generator(simulation_results$initial_conditions$R0[sample[i]]),
      ifr = get_degenerate_random_variate_generator(simulation_results$initial_conditions$ifr[sample[i]]),
      shape_gt = (mean_gt / sd_gt)^2,
      rate_gt = mean_gt / sd_gt^2,
      shape_ip = (mean_ip / sd_ip)^2,
      rate_ip = mean_ip / sd_ip^2,
      meanlog_i2d = 2 * log(mean_i2d) - log(sd_i2d^2 + mean_i2d^2) / 2,
      sdlog_i2d = - 2 * log(mean_i2d) + log(sd_i2d^2 + mean_i2d^2),
      prob_detection_cases = prob_detection_cases,
      prob_detection_deaths = prob_detection_deaths,
      policies = counterfactual_policies,
      initial_state = initial_state
    )$end_results %>%
      dplyr::rename(
        id_simulation = id
      ) %>%
      dplyr::mutate(
        id = sample[i],
        removed_policy = policies[[j]]$long_name
      ) %>%
      dplyr::group_by(id_simulation) %>%
      dplyr::mutate(
        cumulative_cases = cumsum(cases),
        cumulative_deaths = cumsum(deaths)
      )
    
    real_counterfactual_simulation_results <- real_counterfactual_simulation_results %>%
      dplyr::bind_rows(counterfactual_simulation_results)
    
    cat(sprintf("Simulation of counterfactuals for unit %d without policy %s done\r", i, policies[[j]]$long_name))
    Sys.sleep(1)
  }
}

real_counterfactual_simulation_summary <- real_counterfactual_simulation_results %>%
  dplyr::arrange(id, t) %>%
  dplyr::group_by(removed_policy, t, id) %>%
  dplyr::summarize(
    mean_cumulative_cases = mean(cumulative_cases),
    lower_cumulative_cases = quantile(cumulative_cases, probs = 0.025),
    upper_cumulative_cases = quantile(cumulative_cases, probs = 0.975),
    mean_cumulative_deaths = mean(cumulative_deaths),
    lower_cumulative_deaths = quantile(cumulative_deaths, probs = 0.025),
    upper_cumulative_deaths = quantile(cumulative_deaths, probs = 0.975)
  ) %>%
  dplyr::summarize(
    sum_mean_cumulative_cases = sum(mean_cumulative_cases),
    sum_lower_cumulative_cases = sum(lower_cumulative_cases),
    sum_upper_cumulative_cases = sum(upper_cumulative_cases),
    sum_mean_cumulative_deaths = sum(mean_cumulative_deaths),
    sum_lower_cumulative_deaths = sum(lower_cumulative_deaths),
    sum_upper_cumulative_deaths = sum(upper_cumulative_deaths)
  ) %>%
  dplyr::ungroup()

model_factual_simulation_results_cases <- perform_counterfactual_simulation(
  data_for_regression,
  data_for_regression,
  "dlogdc",
  unlist(purrr::map(policies, ~ .x$short_name)),
  unlist(purrr::map(policies, ~ .x$long_name)),
  round(mean_gt + mean_ip),
  difflag,
  nb_counterfactual_simulations
) %>%
  dplyr::mutate(
    dc = exp(logdc)
  ) %>%
  dplyr::group_by(id, id_simulation) %>%
  dplyr::mutate(
    cumulative_cases = as.vector(round(stats::filter(dc, c(rep(0, difflag - 1), 1), method = "recursive", init = rep(0, difflag)))),
    dc = round(dc),
    cases = cumulative_cases - dplyr::lag(cumulative_cases, 1)
  ) %>%
  dplyr::ungroup()

model_factual_simulation_results_deaths <- perform_counterfactual_simulation(
  data_for_regression,
  data_for_regression,
  "dlogdd",
  unlist(purrr::map(policies, ~ .x$short_name)),
  unlist(purrr::map(policies, ~ .x$long_name)),
  round(mean_i2d),
  difflag,
  nb_counterfactual_simulations
) %>%
  dplyr::mutate(
    dd = exp(logdd)
  ) %>%
  dplyr::group_by(id, id_simulation) %>%
  dplyr::mutate(
    cumulative_deaths = as.vector(round(stats::filter(dd, c(rep(0, difflag - 1), 1), method = "recursive", init = rep(0, difflag)))),
    dd = round(dd),
    deaths = cumulative_deaths - dplyr::lag(cumulative_deaths, 1)
  ) %>%
  dplyr::ungroup()

model_factual_simulation_summary_cases <- model_factual_simulation_results_cases %>%
  dplyr::group_by(t, id) %>%
  dplyr::summarize(
    mean_cumulative_cases = mean(cumulative_cases),
    lower_cumulative_cases = quantile(cumulative_cases, probs = 0.025),
    upper_cumulative_cases = quantile(cumulative_cases, probs = 0.975)
  ) %>%
  dplyr::summarize(
    sum_mean_cumulative_cases = sum(mean_cumulative_cases),
    sum_lower_cumulative_cases = sum(lower_cumulative_cases),
    sum_upper_cumulative_cases = sum(upper_cumulative_cases)
  ) %>%
  dplyr::ungroup()

model_factual_simulation_summary_deaths <- model_factual_simulation_results_deaths %>%
  dplyr::group_by(t, id) %>%
  dplyr::summarize(
    mean_cumulative_deaths = mean(cumulative_deaths),
    lower_cumulative_deaths = quantile(cumulative_deaths, probs = 0.025),
    upper_cumulative_deaths = quantile(cumulative_deaths, probs = 0.975)
  ) %>%
  dplyr::summarize(
    sum_mean_cumulative_deaths = sum(mean_cumulative_deaths),
    sum_lower_cumulative_deaths = sum(lower_cumulative_deaths),
    sum_upper_cumulative_deaths = sum(upper_cumulative_deaths)
  ) %>%
  dplyr::ungroup()

model_factual_simulation_summary <- model_factual_simulation_summary_cases %>%
  dplyr::inner_join(model_factual_simulation_summary_deaths, by = "t")

model_counterfactual_simulation_results_cases <- tibble::tibble()
model_counterfactual_simulation_results_deaths <- tibble::tibble()
for (i in 1:length(policies)) {
  counterfactual_data <- data_for_regression %>%
    dplyr::mutate(
      "{policies[[i]]$short_name}" := 0
    )
  
  model_counterfactual_simulation_results_cases <- perform_counterfactual_simulation(
    data_for_regression,
    counterfactual_data,
    "dlogdc",
    unlist(purrr::map(policies, ~ .x$short_name)),
    unlist(purrr::map(policies, ~ .x$long_name)),
    round(mean_gt + mean_ip),
    difflag,
    nb_counterfactual_simulations
  ) %>%
    dplyr::mutate(
      removed_policy = policies[[i]]$long_name,
      dc = exp(logdc)
    ) %>%
    dplyr::group_by(id, id_simulation) %>%
    dplyr::mutate(
      cumulative_cases = as.vector(round(stats::filter(dc, c(rep(0, difflag - 1), 1), method = "recursive", init = rep(0, difflag)))),
      dc = round(dc),
      cases = cumulative_cases - dplyr::lag(cumulative_cases, 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(
      model_counterfactual_simulation_results_cases
    )
  
  model_counterfactual_simulation_results_deaths <- perform_counterfactual_simulation(
    data_for_regression,
    counterfactual_data,
    "dlogdd",
    unlist(purrr::map(policies, ~ .x$short_name)),
    unlist(purrr::map(policies, ~ .x$long_name)),
    round(mean_i2d),
    difflag,
    nb_counterfactual_simulations
  ) %>%
    dplyr::mutate(
      removed_policy = policies[[i]]$long_name,
      dd = exp(logdd)
    ) %>%
    dplyr::group_by(id, id_simulation) %>%
    dplyr::mutate(
      cumulative_deaths = as.vector(round(stats::filter(dd, c(rep(0, difflag - 1), 1), method = "recursive", init = rep(0, difflag)))),
      dd = round(dd),
      deaths = cumulative_deaths - dplyr::lag(cumulative_deaths, 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(
      model_counterfactual_simulation_results_deaths
    )
}

model_counterfactual_simulation_summary_cases <- model_counterfactual_simulation_results_cases %>%
  dplyr::group_by(removed_policy, t, id) %>%
  dplyr::summarize(
    mean_cumulative_cases = mean(cumulative_cases),
    lower_cumulative_cases = quantile(cumulative_cases, probs = 0.025),
    upper_cumulative_cases = quantile(cumulative_cases, probs = 0.975)
  ) %>%
  dplyr::summarize(
    sum_mean_cumulative_cases = sum(mean_cumulative_cases),
    sum_lower_cumulative_cases = sum(lower_cumulative_cases),
    sum_upper_cumulative_cases = sum(upper_cumulative_cases)
  ) %>%
  dplyr::ungroup()

model_counterfactual_simulation_summary_deaths <- model_counterfactual_simulation_results_deaths %>%
  dplyr::group_by(removed_policy, t, id) %>%
  dplyr::summarize(
    mean_cumulative_deaths = mean(cumulative_deaths),
    lower_cumulative_deaths = quantile(cumulative_deaths, probs = 0.025),
    upper_cumulative_deaths = quantile(cumulative_deaths, probs = 0.975)
  ) %>%
  dplyr::summarize(
    sum_mean_cumulative_deaths = sum(mean_cumulative_deaths),
    sum_lower_cumulative_deaths = sum(lower_cumulative_deaths),
    sum_upper_cumulative_deaths = sum(upper_cumulative_deaths)
  ) %>%
  dplyr::ungroup()

model_counterfactual_simulation_summary <- model_counterfactual_simulation_summary_cases %>%
  dplyr::inner_join(model_counterfactual_simulation_summary_deaths, by = c("removed_policy", "t"))

actual_cumulative <- data_for_regression %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(
    cumulative_cases = cumsum(cases),
    cumulative_deaths = cumsum(deaths)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(t) %>%
  dplyr::summarize(
    sum_cumulative_cases = sum(cumulative_cases),
    sum_cumulative_deaths = sum(cumulative_deaths)
  )

real_counterfactual_reality_comparison <- real_counterfactual_simulation_summary %>%
  dplyr::group_by(removed_policy) %>%
  dplyr::mutate(
    difference_sum_mean_cumulative_cases = sum_mean_cumulative_cases - actual_cumulative$sum_cumulative_cases,
    difference_sum_lower_cumulative_cases = sum_lower_cumulative_cases - actual_cumulative$sum_cumulative_cases,
    difference_sum_upper_cumulative_cases = sum_upper_cumulative_cases - actual_cumulative$sum_cumulative_cases,
    difference_sum_mean_cumulative_deaths = sum_mean_cumulative_deaths - actual_cumulative$sum_cumulative_deaths,
    difference_sum_lower_cumulative_deaths = sum_lower_cumulative_deaths - actual_cumulative$sum_cumulative_deaths,
    difference_sum_upper_cumulative_deaths = sum_upper_cumulative_deaths - actual_cumulative$sum_cumulative_deaths
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    removed_policy,
    t,
    difference_sum_mean_cumulative_cases,
    difference_sum_lower_cumulative_cases,
    difference_sum_upper_cumulative_cases,
    difference_sum_mean_cumulative_deaths,
    difference_sum_lower_cumulative_deaths,
    difference_sum_upper_cumulative_deaths
  )

model_counterfactual_reality_comparison <- model_counterfactual_simulation_summary %>%
  dplyr::mutate(
    difference_sum_mean_cumulative_cases = sum_mean_cumulative_cases - model_factual_simulation_summary$sum_mean_cumulative_cases,
    difference_sum_lower_cumulative_cases = sum_lower_cumulative_cases - model_factual_simulation_summary$sum_lower_cumulative_cases,
    difference_sum_upper_cumulative_cases = sum_upper_cumulative_cases - model_factual_simulation_summary$sum_upper_cumulative_cases,
    difference_sum_mean_cumulative_deaths = sum_mean_cumulative_deaths - model_factual_simulation_summary$sum_mean_cumulative_deaths,
    difference_sum_lower_cumulative_deaths = sum_lower_cumulative_deaths - model_factual_simulation_summary$sum_lower_cumulative_deaths,
    difference_sum_upper_cumulative_deaths = sum_upper_cumulative_deaths - model_factual_simulation_summary$sum_upper_cumulative_deaths
  ) %>%
  dplyr::select(
    removed_policy,
    t,
    difference_sum_mean_cumulative_cases,
    difference_sum_lower_cumulative_cases,
    difference_sum_upper_cumulative_cases,
    difference_sum_mean_cumulative_deaths,
    difference_sum_lower_cumulative_deaths,
    difference_sum_upper_cumulative_deaths
  )

real_model_counterfactual_effect_comparison <- tibble::tibble(
  removed_policy = real_counterfactual_reality_comparison$removed_policy,
  t = real_counterfactual_reality_comparison$t,
  real_difference_sum_mean_cumulative_cases = real_counterfactual_reality_comparison$difference_sum_mean_cumulative_cases,
  real_difference_sum_lower_cumulative_cases = real_counterfactual_reality_comparison$difference_sum_lower_cumulative_cases,
  real_difference_sum_upper_cumulative_cases = real_counterfactual_reality_comparison$difference_sum_upper_cumulative_cases,
  model_difference_sum_mean_cumulative_cases = model_counterfactual_reality_comparison$difference_sum_mean_cumulative_cases,
  model_difference_sum_lower_cumulative_cases = model_counterfactual_reality_comparison$difference_sum_lower_cumulative_cases,
  model_difference_sum_upper_cumulative_cases = model_counterfactual_reality_comparison$difference_sum_upper_cumulative_cases,
  model_real_difference_sum_mean_cumulative_cases = model_difference_sum_mean_cumulative_cases / real_difference_sum_mean_cumulative_cases - 1,
  model_real_difference_sum_lower_cumulative_cases = model_difference_sum_lower_cumulative_cases / real_difference_sum_lower_cumulative_cases - 1,
  model_real_difference_sum_upper_cumulative_cases = model_difference_sum_upper_cumulative_cases / real_difference_sum_upper_cumulative_cases - 1,
  real_difference_sum_mean_cumulative_deaths = real_counterfactual_reality_comparison$difference_sum_mean_cumulative_deaths,
  real_difference_sum_lower_cumulative_deaths = real_counterfactual_reality_comparison$difference_sum_lower_cumulative_deaths,
  real_difference_sum_upper_cumulative_deaths = real_counterfactual_reality_comparison$difference_sum_upper_cumulative_deaths,
  model_difference_sum_mean_cumulative_deaths = model_counterfactual_reality_comparison$difference_sum_mean_cumulative_deaths,
  model_difference_sum_lower_cumulative_deaths = model_counterfactual_reality_comparison$difference_sum_lower_cumulative_deaths,
  model_difference_sum_upper_cumulative_deaths = model_counterfactual_reality_comparison$difference_sum_upper_cumulative_deaths,
  model_real_difference_sum_mean_cumulative_deaths = model_difference_sum_mean_cumulative_deaths / real_difference_sum_mean_cumulative_deaths - 1,
  model_real_difference_sum_lower_cumulative_deaths = model_difference_sum_lower_cumulative_deaths / real_difference_sum_lower_cumulative_deaths - 1,
  model_real_difference_sum_upper_cumulative_deaths = model_difference_sum_upper_cumulative_deaths / real_difference_sum_upper_cumulative_deaths - 1
) %>%
  dplyr::select(
    removed_policy,
    t,
    real_difference_sum_mean_cumulative_cases,
    real_difference_sum_lower_cumulative_cases,
    real_difference_sum_upper_cumulative_cases,
    model_difference_sum_mean_cumulative_cases,
    model_difference_sum_lower_cumulative_cases,
    model_difference_sum_upper_cumulative_cases,
    real_difference_sum_mean_cumulative_deaths,
    real_difference_sum_lower_cumulative_deaths,
    real_difference_sum_upper_cumulative_deaths,
    model_difference_sum_mean_cumulative_deaths,
    model_difference_sum_lower_cumulative_deaths,
    model_difference_sum_upper_cumulative_deaths
  )

real_model_counterfactual_effect_comparison_mean <- real_model_counterfactual_effect_comparison %>%
  tidyr::pivot_longer(
    cols = contains("mean"),
    names_to = "tmp",
    values_to = "mean"
  ) %>%
  dplyr::mutate(
    model = ifelse(stringr::str_detect(tmp, "real"), "Actual data generating\nprocess", "Model"),
    outcome = ifelse(stringr::str_detect(tmp, "cases"), "cases", "deaths")
  ) %>%
  dplyr::select(
    removed_policy,
    t,
    outcome,
    model,
    mean
  )

real_model_counterfactual_effect_comparison_lower <- real_model_counterfactual_effect_comparison %>%
  tidyr::pivot_longer(
    cols = contains("lower"),
    names_to = "tmp",
    values_to = "lower"
  ) %>%
  dplyr::mutate(
    model = ifelse(stringr::str_detect(tmp, "real"), "Actual data generating\nprocess", "Model"),
    outcome = ifelse(stringr::str_detect(tmp, "cases"), "cases", "deaths")
  ) %>%
  dplyr::select(
    removed_policy,
    t,
    outcome,
    model,
    lower
  )

real_model_counterfactual_effect_comparison_upper <- real_model_counterfactual_effect_comparison %>%
  tidyr::pivot_longer(
    cols = contains("upper"),
    names_to = "tmp",
    values_to = "upper"
  ) %>%
  dplyr::mutate(
    model = ifelse(stringr::str_detect(tmp, "real"), "Actual data generating\nprocess", "Model"),
    outcome = ifelse(stringr::str_detect(tmp, "cases"), "cases", "deaths")
  ) %>%
  dplyr::select(
    removed_policy,
    t,
    outcome,
    model,
    upper
  )

real_model_counterfactual_effect_comparison_for_plot <- real_model_counterfactual_effect_comparison_mean %>%
  dplyr::inner_join(real_model_counterfactual_effect_comparison_lower, by = c("removed_policy", "t", "outcome", "model")) %>%
  dplyr::inner_join(real_model_counterfactual_effect_comparison_upper, by = c("removed_policy", "t", "outcome", "model"))

last_cases_real_mean <- real_model_counterfactual_effect_comparison$real_difference_sum_mean_cumulative_cases[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_cases_real_lower <- real_model_counterfactual_effect_comparison$real_difference_sum_lower_cumulative_cases[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_cases_real_upper <- real_model_counterfactual_effect_comparison$real_difference_sum_upper_cumulative_cases[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_cases_model_mean <- real_model_counterfactual_effect_comparison$model_difference_sum_mean_cumulative_cases[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_cases_model_lower <- real_model_counterfactual_effect_comparison$model_difference_sum_lower_cumulative_cases[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_cases_model_upper <- real_model_counterfactual_effect_comparison$model_difference_sum_upper_cumulative_cases[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]

last_deaths_real_mean <- real_model_counterfactual_effect_comparison$real_difference_sum_mean_cumulative_deaths[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_deaths_real_lower <- real_model_counterfactual_effect_comparison$real_difference_sum_lower_cumulative_deaths[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_deaths_real_upper <- real_model_counterfactual_effect_comparison$real_difference_sum_upper_cumulative_deaths[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_deaths_model_mean <- real_model_counterfactual_effect_comparison$model_difference_sum_mean_cumulative_deaths[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_deaths_model_lower <- real_model_counterfactual_effect_comparison$model_difference_sum_lower_cumulative_deaths[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]
last_deaths_model_upper <- real_model_counterfactual_effect_comparison$model_difference_sum_upper_cumulative_deaths[seq(from = simulation_length, to = length(policies) * simulation_length, by = simulation_length)]

model_error_cases_annotations <- tibble::tibble(
  model = rep("", 4),
  removed_policy = c("Mandating masks", "Non-essential businesses closure", "Schools closure", "Stay-at-home order"),
  x = rep(simulation_length / 2, 4),
  y = rep(max(last_cases_real_upper, last_cases_model_upper) * 0.9, 4),
  label = paste0(
    "model's error = ",
    ifelse(
      last_cases_model_mean - last_cases_real_mean > 0,
      "+",
      ""
    ),
    round((last_cases_model_mean - last_cases_real_mean) / last_cases_real_mean * 100),
    "% [",
    ifelse(
      last_cases_model_lower - last_cases_real_upper > 0,
      "+",
      ""
    ),
    round((last_cases_model_lower - last_cases_real_upper) / last_cases_real_upper * 100),
    "%, ",
    ifelse(
      last_cases_model_upper - last_cases_real_lower > 0,
      "+",
      ""
    ),
    round((last_cases_model_upper - last_cases_real_lower) / last_cases_real_lower * 100),
    "%]"
  )
)

model_error_deaths_annotations <- tibble::tibble(
  model = rep("", 4),
  removed_policy = c("Mandating masks", "Non-essential businesses closure", "Schools closure", "Stay-at-home order"),
  x = rep(simulation_length / 2, 4),
  y = rep(max(last_deaths_real_upper, last_deaths_model_upper) * 0.9, 4),
  label = paste0(
    "model's error = ",
    ifelse(
      last_deaths_model_mean - last_deaths_real_mean > 0,
      "+",
      ""
    ),
    round((last_deaths_model_mean - last_deaths_real_mean) / last_deaths_real_mean * 100),
    "% [",
    ifelse(
      last_deaths_model_lower - last_deaths_real_upper > 0,
      "+",
      ""
    ),
    round((last_deaths_model_lower - last_deaths_real_upper) / last_deaths_real_upper * 100),
    "%, ",
    ifelse(
      last_deaths_model_upper - last_deaths_real_lower > 0,
      "+",
      ""
    ),
    round((last_deaths_model_upper - last_deaths_real_lower) / last_deaths_real_lower * 100),
    "%]"
  )
)

real_model_counterfactual_effect_comparison_cases_plot <- ggplot(real_model_counterfactual_effect_comparison_for_plot %>% dplyr::filter(outcome == "cases")) +
  geom_line(size = 1, mapping = aes(x = t, y = mean, group = model, color = model)) +
  geom_ribbon(aes(x = t, ymin = lower, ymax = upper, group = model, color = model), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_text(data = model_error_cases_annotations, aes(x = x, y = y, label = label)) +
  scale_color_discrete(name = "") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ removed_policy, ncol = 2) +
  theme_bw() +
  ggtitle("Effect of removing policies on the cumulative number of cases according to the actual data generating process and\na model similar to that used in Chernozhukov et al. (2021)") +
  xlab("Time") +
  ylab("Absolute increase in cumulative number of cases") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  ggsave(
    "Model on simulated data/Figures/Counterfactual experiment for cases with model on simulated data.png",
    width = 18,
    height = 12
  )

real_model_counterfactual_effect_comparison_deaths_plot <- ggplot(real_model_counterfactual_effect_comparison_for_plot %>% dplyr::filter(outcome == "deaths")) +
  geom_line(size = 1, mapping = aes(x = t, y = mean, group = model, color = model)) +
  geom_ribbon(aes(x = t, ymin = lower, ymax = upper, group = model, color = model), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_text(data = model_error_deaths_annotations, aes(x = x, y = y, label = label)) +
  scale_color_discrete(name = "") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ removed_policy, ncol = 2) +
  theme_bw() +
  ggtitle("Effect of removing policies on the cumulative number of deaths according to the actual data generating process and\na model similar to that used in Chernozhukov et al. (2021)") +
  xlab("Time") +
  ylab("Absolute increase in cumulative number of deaths") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

plot_grid(
  real_model_counterfactual_effect_comparison_cases_plot,
  real_model_counterfactual_effect_comparison_deaths_plot,
  labels = c("", ""),
  ncol = 1,
  rel_heights = c(1, 1)
) +
  ggsave(
    "Model on simulated data/Figures/Counterfactual experiment for cases and deaths with model on simulated data.png",
    width = 18,
    height = 24,
    limitsize = FALSE
  )

nb_placebo_tests <- 500

# placebo_analysis_results <- tibble::tibble()
# for (i in 1:nb_placebo_tests) {
#   policy_matrix <- abind(create_policy_matrix(policies, simulation_length, 50)$policy_matrix, along = 1)
#   
#   placebo_data <- data_for_regression
#   
#   for (j in 1:length(policies)) {
#     placebo_data <- placebo_data %>%
#       dplyr::mutate(
#         "{policies[[j]]$short_name}" := policy_matrix[,j]
#       )
#   }
#   
#   fit_summary <- estimate_model(
#     placebo_data,
#     "dlogdc",
#     unlist(purrr::map(policies, ~ .x$short_name)),
#     unlist(purrr::map(policies, ~ .x$long_name)),
#     round(mean_gt + mean_ip)
#   )$fit_summary
#   
#   placebo_analysis_results <- placebo_analysis_results %>%
#     dplyr::bind_rows(fit_summary)
# }

mean_parks_start <- 45
sd_parks_start <- 30
mean_parks_duration <- 60
sd_parks_duration <- 30
parks_start <- purrr::partial(rlnorm, meanlog = 2 * log(mean_parks_start) - log(sd_parks_start^2 + mean_parks_start^2) / 2, sdlog = - 2 * log(mean_parks_start) + log(sd_parks_start^2 + mean_parks_start^2))
parks_duration <- purrr::partial(rlnorm, meanlog = 2 * log(mean_parks_duration) - log(sd_parks_duration^2 + mean_parks_duration^2) / 2, sdlog = - 2 * log(mean_parks_duration) + log(sd_parks_duration^2 + mean_parks_duration^2))

parks <- list(
  long_name = "Closing parks",
  short_name = "parks",
  prob = 1,
  start = parks_start,
  duration = parks_duration,
  effect = 0
)

policies_with_placebo <- append(
  policies,
  list(parks)
)

progressr::with_progress({
  p <- progressr::progressor(steps = length(unique(data_for_regression$id)))
  
  placebo_analysis_results <- furrr::future_map_dfr(
    1:nb_placebo_tests,
    function( placebo_test_id ) {
      p()
      
      policy_matrix <- abind(create_policy_matrix(list(parks), simulation_length, 50)$policy_matrix, along = 1)
      
      placebo_data <- data_for_regression %>%
        dplyr::mutate(
          parks = policy_matrix[,1]
        )
      
      results_cases <- estimate_model(
        placebo_data,
        "dlogdc",
        unlist(purrr::map(policies_with_placebo, ~ .x$short_name)),
        unlist(purrr::map(policies_with_placebo, ~ .x$long_name)),
        round(mean_gt + mean_ip)
      )$fit_summary %>%
        dplyr::mutate(
          outcome = "Weekly cases growth",
          placebo_test_id = placebo_test_id
        )
      
      results_deaths <- estimate_model(
        placebo_data,
        "dlogdd",
        unlist(purrr::map(policies_with_placebo, ~ .x$short_name)),
        unlist(purrr::map(policies_with_placebo, ~ .x$long_name)),
        round(mean_i2d)
      )$fit_summary %>%
        dplyr::mutate(
          outcome = "Weekly deaths growth",
          placebo_test_id = placebo_test_id
        )
      
      dplyr::bind_rows(results_cases, results_deaths)
    },
    .options = furrr::furrr_options(
      seed = rng_seed,
      globals = c(
        "mean_parks_start",
        "sd_parks_start",
        "mean_parks_duration",
        "sd_parks_duration",
        "parks",
        "policies_with_placebo",
        "create_policy_matrix",
        "simulation_length",
        "data_for_regression",
        "estimate_model",
        "p",
        "mean_gt",
        "mean_ip",
        "mean_i2d",
        "abind",
        "felm"
        )
      )
  )
})

(placebo_analysis_cases_summary <- placebo_analysis_results %>%
    dplyr::filter(policy == "Closing parks" & outcome == "Weekly cases growth") %>%
    dplyr::summarize(
      mean_effect = mean(estimate),
      smallest = min(estimate),
      largest = max(estimate),
      percent_significantly_negative = sum(upper < 0) / dplyr::n(),
      percent_significantly_positive = sum(lower > 0) / dplyr::n()
    ))

(placebo_analysis_deaths_summary <- placebo_analysis_results %>%
    dplyr::filter(policy == "Closing parks" & outcome == "Weekly deaths growth") %>%
    dplyr::summarize(
      mean_effect = mean(estimate),
      smallest = min(estimate),
      largest = max(estimate),
      percent_significantly_negative = sum(upper < 0) / dplyr::n(),
      percent_significantly_positive = sum(lower > 0) / dplyr::n()
    ))
