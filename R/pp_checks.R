pp_check_cor_array <- function(predicted, actual, group = NULL, width = c(0.5, 0.8, 0.9,0.95),
                               actual_point_size = 3) {

    ##### Argument checks #####
    if(any(is.na(actual)) || any(is.na(predicted))) {
        stop("NAs in predicted/actual not supported (yet)")
    }
    if(length(dim(predicted)) != 3) {
        stop("Predicted has to be a samples x observations x questions array")
    }
    n_samples <- dim(predicted)[1]
    n_obs <- dim(predicted)[2]
    n_questions <- dim(predicted)[3]
    if(n_obs != dim(actual)[1]) {
        stop("Predicted has to be a samples x observations x questions array and actual
             has to be observations x questions array - number of OBSERVATIONS in model data and provided data differ.")
    }
    if(n_questions != dim(actual)[2]) {
        stop("Predicted has to be a samples x observations x questions array and actual
             has to be observations x questions array - number of QUESTIONS in model data and provided data differ.")
    }

    if(!is.null(group) && length(group) != n_obs) {
        stop("Group has to be NULL or the same size as number of observations")
    }

    if(is.null(group)) {
        group <- rep(1, n_obs)
        show_group <- FALSE
    } else {
        show_group <- TRUE
    }

    ##### Check group names ####
    group_names <- unique(group)
    n_groups <- length(group_names)

    if(is.null(colnames(actual))) {
        question_names <- paste0("q", 1:n_questions)
    } else {
        question_names <- colnames(actual)
    }

    #### Build matrices ####
    # Building matrices of the same size as `cor_actual_matrix` and then converting them to vectors
    # To ensure that things match
    cor_names <- as.character(outer(question_names, question_names, FUN = paste, sep = "_"))
    cors_to_use <- as.logical(outer(question_names, question_names, FUN = "<"))
    cor_names_to_use <- cor_names[cors_to_use]

    ##### Will store per-group summaries ####
    all_pred_dfs <- list()
    all_actual_dfs <- list()

    for(g in 1:n_groups) {

        cor_actual_matrix <- cor(actual[group == group_names[g], ])


        cor_actual <- as.numeric(cor_actual_matrix)[cors_to_use]
        names(cor_actual) <- cor_names_to_use

        if(any(is.na(cor_actual))) {
            warning(paste0("Group = ", group_names[g], ": ", sum(is.na(cor_actual)),
                           " actual correlations are NA, possibly because no within-group variability is observed\n"))
        }

        cor_pred_all <- suppressWarnings(
            apply(predicted[ , group == group_names[g], ], MARGIN = 1, FUN = cor)
        )
        if(any(is.na(cor_pred_all))) {
            warning(paste0("Group = ", group_names[g], ": ", sum(is.na(cor_pred_all)),
                           " predicted correlations are NA, possibly because no within-group variability is predicted\n"))
        }

        cor_pred <- cor_pred_all[cors_to_use, ]
        rownames(cor_pred) <- cor_names_to_use

        # A bit of a hack to easily get multiple quantile summaries with summarise.
        fn_list <- lapply(c(paste0("~ quantile(.x, probs = ",  (1 - width) / 2, ", na.rm = TRUE)"),
                            paste0("~ quantile(.x, probs = ", 1 - (1 - width) / 2, ", na.rm = TRUE)"))
                          , as.formula)
        names(fn_list) <- c(paste0("lower_|_", width), paste0("upper_|_", width ))
        fn_list[[]]

        pred_df <- as.data.frame(t(cor_pred))
        pred_df <- dplyr::summarise(pred_df, dplyr::across(dplyr::everything(), .fns = fn_list, .names = "{.col}_|_{.fn}"))
        pred_df <- tidyr::pivot_longer(pred_df, dplyr::everything(), names_sep = "_\\|_", names_to = c("pair", "bound", ".width"), values_to = "value")
        pred_df <- tidyr::pivot_wider(pred_df, names_from = "bound", values_from = "value")
        pred_df$fake_median = 0.5 #This is just to make ggdist::geom_interval work
        pred_df$group <- group_names[g]

        all_pred_dfs[[g]] <- pred_df

        all_actual_dfs[[g]] <- data.frame(pair = names(cor_actual), actual = cor_actual, group = group_names[g])
    }

    pred_df <- do.call(rbind, all_pred_dfs)
    actual_df <- do.call(rbind, all_actual_dfs)

    pp_plot <- ggplot(pred_df, aes(x = pair, ymin = lower, ymax = upper)) +
        ggdist::geom_interval(aes(y = fake_median)) +
        geom_point(aes(x = pair, y = actual), inherit.aes = FALSE, data = actual_df, shape = 18, size = actual_point_size, position = position_nudge(x = 0.1)) + ggdist::theme_ggdist() +
        scale_color_brewer("PP interval") + scale_y_continuous("Correlation") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    if(show_group) {
        pp_plot <- pp_plot + facet_wrap(~group)
    }

    pp_plot
}


pp_check_cor_wide <- function(predicted, data, answer_cols, group = NULL,
                              actual_point_size = 3) {
    actual <- data %>% dplyr::select(all_of(answer_cols)) %>% as.matrix()
    pp_check_cor_array(predicted, actual, group, actual_point_size = actual_point_size)
}

pp_check_cor_long <- function(predicted, data, answer_col, question_col, obs_id_col, group = NULL,
                              actual_point_size = 3) {
    all_questions <- unique(data[[question_col]])
    n_questions <- length(all_questions)
    n_samples <- dim(predicted)[1]
    all_obs <- unique(data[[obs_id_col]])
    n_obs <- length(all_obs)

    #obs_map <- 1:n_obs
    #names(obs_map) <- all_obs

    predicted_wide <- array(NA_real_, dim = c(n_samples, n_obs, n_questions),
                            dimnames = list(NULL, NULL, all_questions))
    actual_wide <- array(NA_real_, dim = c(n_obs, n_questions),
                         dimnames = list(NULL, all_questions))
    group_wide <- NULL

    for(q_id in 1:n_questions) {
        #q_data <- data[data[[question_col]] == all_questions[q_id],]
        q_indices_l <- data[[question_col]] == all_questions[q_id]
        if(!identical(data[[obs_id_col]][q_indices_l], all_obs)) {
            stop(paste0("Question '", all_questions[q_id] , "': Need to implement reordering"))
        }
        q_indices <- which(q_indices_l)

        predicted_wide[,,q_id] <- predicted[, q_indices]
        actual_wide[,q_id] <- data[[answer_col]][q_indices]

        if(is.null(group_wide)) {
            group_wide <- group[q_indices]
        } else {
            if(!identical(group_wide, group[q_indices])) {
                stop("Inconsistent grouping")
            }
        }
    }

    pp_check_cor_array(predicted_wide, actual_wide, group_wide, actual_point_size = actual_point_size)
}


pp_check_cor_long2 <- function(model, data = model$data, answer_col, question_col,
                               obs_id_col, timevar, groupvar=NULL, n_samples=4000,
                               actual_point_size = 3) {
    # Get needed data from arguments ####
    #data <- model$data
    all_questions <- unique(data[[question_col]])
    n_questions <- length(all_questions)
    all_obs <- unique(data[[obs_id_col]])
    n_obs <- length(all_obs)

    if( !is.null(groupvar) ){
        if( !groupvar %in% names(data) ){
            stop("specified group variable is not in the data")
        }
    }

    # Make predictions ####
    predictions <- posterior_predict(model, nsamples = n_samples)

    # Create a list of datasets for each sample
    post_datasets <- lapply(1:dim(predictions)[1], function(i) {
        data.frame(data, pred = predictions[i, ])
    })

    # Reshape each datasets predictions from long to wide
    preds_wide <- list()
    for(i in 1:n_samples){
        preds_wide[[i]] <- post_datasets[[i]] %>% pivot_wider(id_cols = c(obs_id_col, all_of(timevar)),
                                       names_from = question_col,
                                       values_from = c("pred"))
    }
    predicted_wide <- array(NA_real_, dim = c(n_samples, nrow(data)/n_questions, n_questions),
                            dimnames = list(NULL, NULL, all_questions))
    for(i in 1:n_samples){
        predicted_wide[i, , ] <- as.matrix(preds_wide[[i]][, all_questions])
    }

    # Reshape just once for actual_wide
    if( is.null(groupvar) ){
        actual_wide <- post_datasets[[i]] %>% pivot_wider(id_cols = c(obs_id_col, all_of(timevar)),
                                       names_from = question_col,
                                       values_from = c(answer_col)) %>%
            data.frame()
        group <- NULL
    } else {
        actual_wide <- post_datasets[[i]] %>% pivot_wider(id_cols = c(obs_id_col, all_of(timevar),
                                                                      all_of(groupvar)),
                                                          names_from = question_col,
                                                          values_from = c(answer_col)) %>%
            data.frame()
        group <- as.character(actual_wide[, groupvar])
    }

    # Call the pp_check_cor_array function
    pp_check_cor_array(predicted_wide, actual_wide[, all_questions], group=group,
                       actual_point_size = actual_point_size)
}










