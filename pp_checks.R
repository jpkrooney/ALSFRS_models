pp_check_cor_array <- function(predicted, actual, group = NULL, width = c(0.5, 0.8, 0.9,0.95)) {
    if(any(is.na(actual)) || any(is.na(predicted))) {
        stop("NAs in predicted/actual not supported (yet)")
    }
    if(length(dim(predicted)) != 3) {
        stop("Predicted has to be a samples x observations x questions array")
    }
    n_samples <- dim(predicted)[1]
    n_obs <- dim(predicted)[2]
    n_questions <- dim(predicted)[3]
    if(n_obs != dim(actual)[1] || n_questions != dim(actual)[2]) {
        stop("Predicted has to be a samples x observations x questions array and actual
             has to be observations x questions array")
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

    group_names <- unique(group)
    n_groups <- length(group_names)

    if(is.null(colnames(actual))) {
        question_names <- paste0("q", 1:n_questions)
    } else {
        question_names <- colnames(actual)
    }

    # Building matrices of the same size as `cor_actual_matrix` and then converting them to vectors
    # To ensure that things match
    cor_names <- as.character(outer(question_names, question_names, FUN = paste, sep = "_"))
    cors_to_use <- as.logical(outer(question_names, question_names, FUN = "<"))
    cor_names_to_use <- cor_names[cors_to_use]

    # Will store per-group summaries
    all_pred_dfs <- list()
    all_actual_dfs <- list()

    for(g in 1:n_groups) {

        cor_actual_matrix <- cor(actual[group == group_names[g], ])


        cor_actual <- as.numeric(cor_actual_matrix)[cors_to_use]
        names(cor_actual) <- cor_names_to_use

        cor_pred_all <- apply(predicted[ , group == group_names[g], ], MARGIN = 1, FUN = cor)

        cor_pred <- cor_pred_all[cors_to_use, ]
        rownames(cor_pred) <- cor_names_to_use

        fn_list <- lapply(c(paste0("~ quantile(.x, probs = ",  (1 - width) / 2, ")"),
                            paste0("~ quantile(.x, probs = ", 1 - (1 - width) / 2, ")"))
                          , as.formula)
        names(fn_list) <- c(paste0("lower_|_", width), paste0("upper_|_", width ))
        fn_list[[]]

        pred_df <- as.data.frame(t(cor_pred_to_use))
        pred_df <- dplyr::summarise(pred_df, dplyr::across(dplyr::everything(), .fns = fn_list, .names = "{.col}_|_{.fn}"))
        pred_df <- tidyr::pivot_longer(pred_df, dplyr::everything(), names_sep = "_\\|_", names_to = c("pair", "bound", ".width"), values_to = "value")
        pred_df <- tidyr::pivot_wider(pred_df, names_from = "bound", values_from = "value")
        pred_df$fake_median = 0.5 #This is just to make ggdist::geom_interval work
        pred_df$group <- group_names[g]

        all_pred_dfs[[g]] <- pred_df

        all_actual_dfs[[g]] <- data.frame(pair = names(cor_actual_to_use), actual = cor_actual_to_use, group = group_names[g])
    }

    pred_df <- do.call(rbind, all_pred_dfs)
    actual_df <- do.call(rbind, all_actual_dfs)

    pp_plot <- ggplot(pred_df, aes(x = pair, ymin = lower, ymax = upper)) + ggdist::geom_interval(aes(y = fake_median)) +
        geom_point(aes(x = pair, y = actual), inherit.aes = FALSE, data = actual_df, shape = 18, size = 10, position = position_nudge(x = 0.1)) + ggdist::theme_ggdist() +
        scale_color_brewer("PP interval") + scale_y_continuous("Correlation")

    if(show_group) {
        pp_plot <- pp_plot + facet_wrap(~group)
    }

    pp_plot
}


pp_check_cor_wide <- function(predicted, data, answer_cols, group = NULL) {
    actual <- data %>% select(all_of(answer_cols)) %>% as.matrix()
    pp_check_cor_array(predicted, actual, group)
}

pp_check_cor_long <- function(predicted, data, answer_col, question_col, obs_id_col, group = NULL) {
    all_questions <- unique(data[[question_col]])
    n_questions <- length(all_questions)
    n_samples <- dim(predicted)[1]
    all_obs <- unique(data[[obs_id_col]])
    n_obs <- length(all_obs)

    obs_map <- 1:n_obs
    names(obs_map) <- all_obs

    predicted_wide <- array(NA_real_, dim = c(n_samples, n_obs, n_questions), dimnames = list(NULL, NULL, all_questions))
    actual_wide <- array(NA_real_, dim = c(n_obs, n_questions), dimnames = list(NULL, all_questions))
    group_wide <- NULL

    for(q_id in 1:n_questions) {
        #q_data <- data[data[[question_col]] == all_questions[q_id],]
        q_indices_l <- data[[question_col]] == all_questions[q_id]
        if(!identical(data[[obs_id_col]][q_indices_l], all_obs)) {
            stop("Need to implement reordering")
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

    pp_check_cor_array(predicted_wide, actual_wide, group_wide)
}
