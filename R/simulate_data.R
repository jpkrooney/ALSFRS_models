simulate_data_from_registry <- function(standardized_data,
                                         n_patients_per_group,
                                         min_measurements_per_patient,
                                         max_duration,
                                         effect_prob){
    n_patients <- n_patients_per_group * 2
    selected_patients <- standardized_data %>%
        filter(alsfrs_dly_mnths < max_duration) %>%
        group_by(subject_id) %>%
        filter(n() >= min_measurements_per_patient) %>%
        ungroup() %>%
        filter(subject_id %in% sample(unique(subject_id), n_patients, replace = FALSE))

    if(length(unique(selected_patients$subject_id)) != n_patients) {
        stop("Bad selection")
    }

    # Assign subjects into random pairs and with prob = effect_prob pick the one with
    # better outcome (lower decrease in score at the end)
    # into the treatment group
    grouped_patients <- selected_patients %>%
        mutate(subject_id = as.integer(factor(subject_id))) %>%
        group_by(subject_id) %>%
        mutate(
            alsfrs_start = alsfrs_total[which.min(alsfrs_dly_mnths)],
            alsfrs_end = alsfrs_total[which.max(alsfrs_dly_mnths)],
            total_decrease = alsfrs_start - alsfrs_end
        ) %>%
        ungroup() %>%
        mutate(pair = sample(rep(1:n_patients_per_group, each = 2),
                             size = n_patients_per_group * 2)[subject_id]) %>%
        group_by(pair) %>%
        mutate(pick_better = rbinom(1, size = 1, prob = effect_prob) == 1,
               treatment_subject_id = if_else(pick_better,
                                              unique(subject_id[which.min(total_decrease)]),
                                              unique(subject_id[which.max(total_decrease)])),
               group = if_else(subject_id == treatment_subject_id, "Treatment", "Control")) %>%
        ungroup()

    if(length(unique(grouped_patients$subject_id)) != n_patients) {
        stop("Bad selection")
    }

    group_sizes <- grouped_patients  %>%
        group_by(group) %>%
        summarise(count = length(unique(subject_id)), .groups = "drop") %>%
        pull(count)

    if(length(group_sizes) != 2 || !all(group_sizes == n_patients_per_group)) {
        stop("Bad grouping")
    }


    grouped_patients
}
