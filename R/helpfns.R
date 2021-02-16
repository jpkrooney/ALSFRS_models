# define helper functions to extract correlations
# helper functions
predcor_dim1 <- function( df_long, fitobj, timevar, times = 6:36, nsamp = 2000){

    all_questions <- sort(unique(df_long$question))
    n_reps <- 100

    new_data <- crossing(ID = paste0("new_patient", 1:n_reps),
                         question = all_questions, time = times)
    new_data <- new_data %>% rename(!!timevar := "time")
    # get linear predictions
    linpreds <- posterior_linpred(fitobj,
                                  newdata = new_data,
                                  allow_new_levels = TRUE,
                                  nsamples = nsamp)

    # make empty array to hold correlations
    per_sample_cor <- array(NA_real_, c(nsamp, length(all_questions), length(all_questions)  ))
    #
    for(s in 1:nsamp) {
        predictors_wide <- array(NA_real_, c(n_reps * length(times), length(all_questions)))
        for(q in 1:length(all_questions)) {
            predictors_wide[ , q] <- linpreds[s,
                                              new_data$question == all_questions[q]]
        }
        per_sample_cor[s, , ] <- cor(predictors_wide) #cov2cor(cov(predictors_wide))
    }

    apply(per_sample_cor, MARGIN = c(2,3), FUN = mean)
}

predcor_multv <- function( fitobj, timevar, times = 6:36, nsamp = 2000){

    n_reps <- 100
    new_data <- crossing(ID = paste0("new_patient", 1:n_reps), time = times)
    new_data <- new_data %>% rename(!!timevar := "time")

    # get linear predictions
    linpreds <- posterior_linpred(fitobj,
                                  newdata = new_data,
                                  allow_new_levels = TRUE,
                                  nsamples = nsamp)

    # Create a list of datasets for each sample
    post_datasets <- lapply(1:dim(linpreds)[1], function(i) {
        data.frame(linpreds[ i, ,])
    })

    # calculate correlation of outcomes overall
    cor_all <- lapply(1:nsamp, function(i)
        cor(post_datasets[[i]] ))
    apply(simplify2array(cor_all), 1:2, mean)
}


predcor_cratio<- function( df_long, fitobj, timevar, times = 6:36, nsamp = 2000){

    all_questions <- sort(unique(df_long$question))
    n_reps <- 100

    new_data <- crossing(ID = paste0("new_patient", 1:n_reps),
                         question = all_questions, time = times)
    new_data <- new_data %>% rename(!!timevar := "time")

    # get linear predictions
    linpreds <- posterior_linpred(fitobj,
                                  newdata = new_data,
                                  allow_new_levels = TRUE,
                                  nsamples = nsamp)

    # make empty array to hold correlations
    per_sample_cor <- array(NA_real_, c(nsamp, length(all_questions), length(all_questions)  ))
    #
    for(s in 1:nsamp) {
        predictors_wide <- array(NA_real_, c(n_reps * length(times), length(all_questions)))
        for(q in 1:length(all_questions)) {
            predictors_wide[ , q] <- linpreds[s,
                                              new_data$question == all_questions[q], 1]
            # NOTE - the only difference with dim1 function is ,1]. cratio model outputs a 3D array with a 3rd dimensions of size 4. On inspection each eleemnt 1 to 4 appears to be the same and so I simply subset it. This may need revision.

        }
        per_sample_cor[s, , ] <- cor(predictors_wide) #cov2cor(cov(predictors_wide))
    }

    apply(per_sample_cor, MARGIN = c(2,3), FUN = mean)
}



