simulate_data_from_registry <- function(standardized_data,
                                         n_subjects_per_group,
                                         min_measurements_per_subject,
                                         max_measurements_per_subject,
                                         max_duration,
                                         effect_prob){
    n_subjects <- n_subjects_per_group * 2
    #TODO add time from diagnosis, don't always filter on dly_mnts > 0
    selected_data <- standardized_data %>%
        filter(alsfrs_dly_mnths >= 0, alsfrs_dly_mnths < max_duration) %>%
        group_by(subject_id) %>%
        filter(n() >= min_measurements_per_subject) %>%
        ungroup() %>%
        filter(subject_id %in% sample(unique(subject_id), n_subjects, replace = FALSE)) %>%
        mutate(subject_id = as.integer(factor(subject_id)))

    #Subsample to max_measurements_per_subject
    subsampled_data <- selected_data %>%
        group_by(subject_id) %>%
        mutate(            is_first = alsfrs_dly_mnths == min(alsfrs_dly_mnths),
                           is_last = alsfrs_dly_mnths == min(alsfrs_dly_mnths)
        ) %>%
        group_by(subject_id, is_first, is_last) %>%
        sample_n(min(max_measurements_per_subject, n())) %>%
        ungroup()


    if(length(unique(subsampled_data$subject_id)) != n_subjects) {
        stop("Bad selection")
    }

    # Assign subjects into pairs with equal or almost equal ALSFRS at start
    # and with prob = effect_prob pick the one with
    # better outcome (higher slope of begin -> end trajectory)
    # into the treatment group
    patient_stats <- subsampled_data %>%
        group_by(subject_id) %>%
        summarise(
            alsfrs_start = alsfrs_total[which.min(alsfrs_dly_mnths)],
            alsfrs_end = alsfrs_total[which.max(alsfrs_dly_mnths)],
            alsfrs_change = alsfrs_end - alsfrs_start,
            duration = max(alsfrs_dly_mnths) - min(alsfrs_dly_mnths),
            slope = alsfrs_change / duration,
        ) %>%
        arrange(alsfrs_start) %>%
        mutate(pair = rep(1:n_subjects_per_group, each = 2)) %>%
        group_by(pair) %>%
        mutate(pick_better = rbinom(1, size = 1, prob = effect_prob) == 1,
               treatment_subject_id = if_else(pick_better,
                                              unique(subject_id[which.max(slope)]),
                                              unique(subject_id[which.min(slope)])),
               group = if_else(subject_id == treatment_subject_id, "Treatment", "Control")) %>%
        ungroup()


    grouped_subjects <- subsampled_data %>%
        inner_join(patient_stats, by = "subject_id") %>%
        ungroup()

    if(length(unique(grouped_subjects$subject_id)) != n_subjects) {
        stop("Bad selection")
    }

    group_sizes <- grouped_subjects  %>%
        group_by(group) %>%
        summarise(count = length(unique(subject_id)), .groups = "drop") %>%
        pull(count)

    if(length(group_sizes) != 2 || !all(group_sizes == n_subjects_per_group)) {
        stop("Bad grouping")
    }


    grouped_subjects
}
