#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List simulate_epidemic(
    const int id,
    const int simulation_length,
    const int population,
    const double seed_rate,
    const int seed_length,
    const double R0,
    const double ifr,
    const double shape_gt,
    const double rate_gt,
    const double shape_ip,
    const double rate_ip,
    const double meanlog_i2d,
    const double sdlog_i2d,
    const double prob_detection_cases,
    const double prob_detection_deaths,
    std::vector<std::string> col_names,
    const arma::mat features_matrix,
    const arma::vec effects,
    const Rcpp::Nullable<Rcpp::IntegerVector> initial_state_ = R_NilValue
) {
    unsigned int start;
    Rcpp::IntegerVector seed, infections, cases, deaths, recorded_deaths;
    arma::cube snapshots;

    arma::vec policy_effect_by_time_ = features_matrix * effects;
    Rcpp::NumericVector policy_effect_by_time = Rcpp::NumericVector(policy_effect_by_time_.begin(), policy_effect_by_time_.end());
    Rcpp::NumericVector R = R0 + policy_effect_by_time;
    Rcpp::NumericVector effective_R = Rcpp::clone(R);

    // I found a way to use Rcpp::Nullable with arma::cube by reading https://gallery.rcpp.org/articles/optional-null-function-arguments/
    // on how to have nullable optional arguments in Rcpp functions and https://www.mjdenny.com/Rcpp_Intro.html on how to create an
    // arma::cube from a NumericVector (see also https://stackoverflow.com/a/44231385 on this issue)
    if (initial_state_.isNotNull()) {
        Rcpp::NumericVector initial_state__ = Rcpp::as<Rcpp::NumericVector>(initial_state_);
        Rcpp::IntegerVector dim = initial_state__.attr("dim");
        arma::cube initial_state(initial_state__.begin(), dim[0], dim[1], dim[2], false);
        start = initial_state.n_slices;
        infections = Rcpp::IntegerVector(initial_state.slice(start - 1).col(0).begin(), initial_state.slice(start - 1).col(0).end());
        cases = Rcpp::IntegerVector(initial_state.slice(start - 1).col(1).begin(), initial_state.slice(start - 1).col(1).end());
        deaths = Rcpp::IntegerVector(initial_state.slice(start - 1).col(2).begin(), initial_state.slice(start - 1).col(2).end());
        recorded_deaths = Rcpp::NumericVector(initial_state.slice(start - 1).col(3).begin(), initial_state.slice(start - 1).col(3).end());
        if (start <= seed_length) {
            seed = Rcpp::rpois(seed_length - start, seed_rate);
            for (int i = start; i < seed_length; i++) {
                infections[i] = infections[i] + seed[i - start];
            }
        }
        arma::cube empty_snapshots(simulation_length, 4, simulation_length - start, arma::fill::zeros);
        snapshots = join_slices(initial_state, empty_snapshots);
    } else {
        start = 0;
        seed = Rcpp::rpois(seed_length, seed_rate);
        infections = Rcpp::rep(0, simulation_length);
        for (int i = 0; i < seed_length; i++) {
            infections[i] = seed[i];
        }
        cases = Rcpp::IntegerVector(simulation_length);
        deaths = Rcpp::IntegerVector(simulation_length);
        recorded_deaths = Rcpp::IntegerVector(simulation_length);
        snapshots = arma::cube(simulation_length, 4, simulation_length, arma::fill::zeros);
    }

    for (int i = start; i < simulation_length; i++) {
        // if I don't store the result of Rcpp::sum(infections[Rcpp::Range(0,i)]) into a variable declared as double,
        // but directly use the expression Rcpp::sum(infections[Rcpp::Range(0,i)]) / population to compute effective_R,
        // it's converted to an integer because C++ performs integer division when you use "/" with 2 integers, so
        // 1 - current_total_infections / population is always zero and effective_R doesn't have the right value
        double current_total_infections = Rcpp::sum(infections[Rcpp::Range(0,i)]);
        effective_R[i] *= 1 - current_total_infections / population;
        int nb_infections = infections[i];
        if (nb_infections > 0 && effective_R[i] > 0) {
            Rcpp::NumericVector nb_secondary_infections = Rcpp::rnbinom_mu(nb_infections, 0.1,  effective_R[i]);
            for (int j = 0; j < nb_infections; j++) {
                if (nb_secondary_infections[j] > 0) {
                    Rcpp::NumericVector generation_times = Rcpp::round(Rcpp::rgamma(nb_secondary_infections[j], shape_gt, 1 / rate_gt), 0);
                    Rcpp::NumericVector secondary_infections_times = i + generation_times;

                    for (int k = 0; k < nb_secondary_infections[j]; k++) {
                        if (secondary_infections_times[k] < simulation_length) {
                            infections[secondary_infections_times[k]]++;
                        }
                    }
                }
            }
        }

        snapshots.slice(i).col(0) = Rcpp::as<arma::vec>(infections);
    }

    for (int i = start; i < simulation_length; i++) {
        if (infections[i] > 0) {
            int nb_detected_infections = R::rbinom(infections[i], prob_detection_cases);
            if (nb_detected_infections > 0) {
                Rcpp::NumericVector incubation_periods = Rcpp::round(Rcpp::rgamma(nb_detected_infections, shape_ip, 1 / rate_ip), 0);
                Rcpp::NumericVector symptoms_onset_times = i + incubation_periods;
                for (int j = 0; j < nb_detected_infections; j++) {
                    if (symptoms_onset_times[j] < simulation_length) {
                        cases[symptoms_onset_times[j]]++;
                    }
                }
            }
        }

        snapshots.slice(i).col(1) = Rcpp::as<arma::vec>(cases);
    }

    for (int i = start; i < simulation_length; i++) {
        if (infections[i] >= 1) {
            int nb_fatal = R::rbinom(infections[i], ifr);
            if (nb_fatal > 0) {
                Rcpp::NumericVector death_lags = Rcpp::round(Rcpp::rlnorm(nb_fatal, meanlog_i2d, sdlog_i2d), 0);
                for (int j = 0; j < nb_fatal; j++) {
                    double death_time = i + death_lags[j];
                    if (death_time < simulation_length) {
                        deaths[death_time]++;
                    }
                }
            }
        }

        snapshots.slice(i).col(2) = Rcpp::as<arma::vec>(deaths);
    }

    for (int i = start; i < simulation_length; i++) {
        if (deaths[i] > 0) {
            recorded_deaths[i] = R::rbinom(deaths[i], prob_detection_deaths);
        }

        snapshots.slice(i).col(3) = Rcpp::as<arma::vec>(recorded_deaths);
    }

    Rcpp::NumericVector cumulative_infections = Rcpp::cumsum(Rcpp::NumericVector(infections.begin(), infections.end()));

    // there is no pounsigned int in creating a DataFrame since it would be turned unsigned into a List when I add
    // columns for policies below, see https://stackoverflow.com/a/59369233 on this pounsigned int
    Rcpp::List simulation = Rcpp::List::create(
        Rcpp::Named("id") = Rcpp::IntegerVector(simulation_length, id),
        Rcpp::Named("t") = Rcpp::seq(1, simulation_length),
        Rcpp::Named("R0") = R,
        Rcpp::Named("R") = effective_R,
        Rcpp::Named("infections") = infections,
        Rcpp::Named("cases") = cases,
        Rcpp::Named("deaths") = deaths,
        Rcpp::Named("recorded_deaths") = recorded_deaths
        );

    for (int i = 0; i < col_names.size(); i++) {
        simulation.push_back(features_matrix.col(i), col_names[i]);
    }

    return Rcpp::List::create(Rcpp::Named("end_result") = simulation, Rcpp::Named("snapshots") = snapshots);
}