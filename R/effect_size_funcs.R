#' @title Calculate posterior effect sizes from ALSFRS-r model fits
#' @description Function to calculate posterior effect sizes from ALSFRS-r model fits given start and stop times and summation parameters
#' @param fit a brmsfit object - currently must be a cumulative logit or mv-probit fit
#' @param time1 start time for effect size calculation. Can be single time for all subjects or individualised start time
#' @param time2 stop time for effect size calculation. Can be single time for all subjects or individualised start time
#' @param timevar name for the time variable in fit
#' @param method Determines method used to perform posterio predictions. Default value is "expected". "expected" or "predicted" are allowed.
#' @param summation This tells `effect_from_sim_study()` how to sum across ALSFRS-r subscores. Default is "sum_allQ" which performs effect size calculation using total ALSFRS-r. Other allowed values are "sum_indQ", "sum_subscores". "sum_indQ" calculates effect size for individual questions. "sum_subscores" calculates effect size for ALSFRS-r subscores in 4 domains - bulbar, fine motor, gross motor and respiratory.
#' @param model This tells `effect_from_sim_study()` whether the fit id a logit or mv-probit model
#'
#' @return Returns calculated effect size as a vector when `summation = "sum_allQ"` or as a matrix for other options
#'

effect_from_sim_study <- function(fit, time1, time2, timevar,
                                  method = c("expected", "predicted"),
                                  summation = c("sum_allQ", "sum_indQ", "sum_subscores"),
                                  model = c("mv-probit", "logit")) {
    # resolve summing method
    method <- match.arg(method)
    summation <- match.arg(summation)

    funcs <-list(
        sum_allQ = if(method == "expected") {expected_sum_allQ} else if(method == "predicted") {rowSums},
        sum_indQ = if(method == "expected") {expected_sum_indQ} else if(method == "predicted") {identity},
        sum_subscores = if(method == "expected") {expected_sum_ALSFRSdomains} else if(method == "predicted") {sum_alsfrs_domains}
    )

    # make data for prediction and rename time variable
    predict_data <- data.frame(group = c("Control", "Control", "Treatment", "Treatment"),
                               time = c(time1, time2, time1, time2)) |>
        rename(!!timevar := "time")

    if(method == "expected"){
        # predict
        pred <- posterior_epred(fit, newdata = predict_data, re_formula = NA)

    } else if (method == "predicted"){
        for (q in paste0("Q", sprintf("%02d", 1:12))) {
            predict_data[[q]] <- 1   # or any valid category
        }

        pred <- if(model == "mv-probit") {
            posterior_predict_mv_probit(fit, newdata = predict_data, re_formula = NA)
        } else if (model == "logit") {
            posterior_predict(fit, newdata = predict_data, re_formula = NA)
        }
        pred <- pred - 1
    }

    # sums based on method and summation
    start_control <- funcs[[summation]](pred[,1,])
    end_control <- funcs[[summation]](pred[,2,])
    start_treatment <- funcs[[summation]](pred[,3,])
    end_treatment <- funcs[[summation]](pred[,4,])

    # calculate effect
    (end_treatment - start_treatment) - (end_control - start_control)
}


#### effect_from_sim_study() helper functions ####
# core expected function - sums across thresholds of individual question ####
expected_sum <- function(x) {
    col_values <- as.integer(colnames(x))
    sweep(x, MARGIN = 2, STATS = col_values, FUN = "*")
}

# sum across all questions
expected_sum_allQ <- function(x){
    rowSums(expected_sum(x))
}

# sum individual questions
expected_sum_indQ <- function(x){

    levs <- length(unique(colnames(x)))
    qs <- dim(x)[2] / levs
    expected_sum(x) |>
        array(dim = c(nrow(x), levs, qs)) |> # note only valid if each question has same number of levels
        aperm(c(1,3,2)) |>
        rowSums(dim = 2) - 1
}

# sum typical alsfrs subscores domains
sum_alsfrs_domains <- function(x){
    bulbar <- rowSums(x[, 1:3])
    finemotor <- rowSums(x[, 4:6])
    grossmotor <- rowSums(x[, 7:9])
    resp <- rowSums(x[, 10:12])
    data.frame(bulbar, grossmotor, finemotor, resp)
}

# sum typical alsfrs subscores domains (for expected calculations only)
expected_sum_ALSFRSdomains <- function(x){
    indsums <- expected_sum_indQ(x)
    sum_alsfrs_domains(indsums)
}

