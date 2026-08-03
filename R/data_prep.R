data_prep <- function(proact) {


    # DATA PREPARATION
    # unify Q5
    # omit vars that are not useful
    # unify question names
    # keep only complete obs

    days_per_month <- 365.25 / 12

    any(!is.na(proact$Q5a_Cutting_without_Gastrostomy) & !is.na(proact$Q5b_Cutting_with_Gastrostomy) &
            proact$Q5a_Cutting_without_Gastrostomy != proact$Q5b_Cutting_with_Gastrostomy)
    #TODO: need to resolve this

    proact_full <- proact %>%
        mutate(Q5_Cutting_Combined = if_else(is.na(Q5a_Cutting_without_Gastrostomy),
                                             Q5b_Cutting_with_Gastrostomy, Q5a_Cutting_without_Gastrostomy)) %>%
        select(-Q5b_Cutting_with_Gastrostomy, -Q5a_Cutting_without_Gastrostomy,
               -Q10_Respiratory, -ALSFRS_Total, -Mode_of_Administration, -ALSFRS_Responded_By) %>%
        filter(if_all(c(starts_with("Q"), starts_with("R")), function(x) { !is.na(x) } ))


    question_col_names_proact <- sort(names(proact_full)[grepl("^[QR]", names(proact_full))])
    question_col_names <- sprintf("Q%02d", 1:12)
    question_col_names_translation <- sprintf("Q%02d", 1:12)
    names(question_col_names_translation) <- question_col_names_proact

    proact_standardized <- proact_full %>%
        rename_with(function(x) { question_col_names_translation[x]}, .cols = all_of(question_col_names_proact)) %>%
        rename(alsfrs_total = ALSFRS_R_Total) %>%
        mutate(alsfrs_dly_mnths = ALSFRS_Delta / days_per_month) %>%
        select(subject_id, all_of(question_col_names), alsfrs_total, alsfrs_dly_mnths)


    proact_standardized
}

