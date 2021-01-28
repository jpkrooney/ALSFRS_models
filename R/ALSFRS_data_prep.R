# Load required packages
library(tidyverse)

set.seed(3)

# Load ALS register data

df1 <- readRDS("../../Irish ALS Register Data/Data/Formatted data/IRL_cases_1995_to_16Jan2020_updated_16Jan2020.rds")
df_frs <- readRDS("../../Irish ALS Register Data/Data/Formatted data/IRL_ALSFRS_updated_29Oct2020.rds")

# Make a Q5 joining 5A and 5B
df_frs$ALSFRS_Q5 <- ifelse( is.na(df_frs$ALSFRS_Q5a), df_frs$ALSFRS_Q5b, df_frs$ALSFRS_Q5a)
df_frs$ALSFRS_Q5a <- NULL
df_frs$ALSFRS_Q5b <- NULL

# Make var with 12 Q names
frs_qs <- c("ALSFRS_Q1", "ALSFRS_Q2", "ALSFRS_Q3", "ALSFRS_Q4", "ALSFRS_Q5", "ALSFRS_Q6",
            "ALSFRS_Q7", "ALSFRS_Q8", "ALSFRS_Q9", "ALSFRS_Q10", "ALSFRS_Q11", "ALSFRS_Q12")

# Limit cohort to only those with complete ALSFRS scores
df_frs <- df_frs[ rowSums(is.na(df_frs[, frs_qs])) == 0 ,  ]
df1 <- df1[df1$ID %in% df_frs$ID, ]

# Remove any aged under 40 to protect privacy
df1 <- df1[ !is.na(df1$age_dx) & df1$age_dx > 40, ]
# Remove any aged over 85 also to protect privacy
df1 <- df1[ !is.na(df1$age_dx) & df1$age_dx <= 85, ]


write.csv(df1, file = "Data/casedata.csv", row.names = FALSE)
write.csv(df_frs, file = "Data/longdata.csv", row.names = FALSE)
