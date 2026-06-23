#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    24-Mar-2025
# Author:  Robert Fletcher
# Purpose: Temporary script file to import pregnancy data from CCU018
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))
library(DBI)


# Establish Databricks connection -----------------------------------------

con <- dbConnect(
  odbc::odbc(), dsn = data_source, HTTPPath = path, PWD = password
)


# Import data and write to project folder ---------------------------------

# Connect to Databricks

# Load HES APC MAT deliveries dataset
hes_apc_mat_tbl <- tbl(
  con, 
  dbplyr::in_schema(
    "", "ccu018_02_tmp_hes_apc_mat_del_clean"
  )
)

hes_apc_mat <- hes_apc_mat_tbl %>% collect()

# Load cohort dataset
cohort_tbl <- tbl(
  con, 
  dbplyr::in_schema("", "ccu018_02_out_cohort")
)

cohort <- cohort_tbl %>% collect()

# Filter out records with implausible gestational age
hes_apc_mat_2 <- hes_apc_mat %>% select(EPIKEY, PROCODE3)

cohort_tmp <- cohort %>% left_join(hes_apc_mat_2, by = "EPIKEY")

cohort_filtered <- cohort_tmp %>%
  mutate(
    implausible_gest_preg = ifelse(
      source == "HES_MAT" & 
        PROCODE3 %in% c('R0A', 'RDU', 'RXU', 'RTD', 'RAE', 'R36'),
      TRUE,
      FALSE
    )
  ) %>%
  filter(implausible_gest_preg == FALSE)

# Identify the number of filtered out records
check <- anti_join(cohort, cohort_filtered, by = "PERSON_ID")

# Update cohort dataset
cohort_2 <- cohort_filtered |>
  distinct(person_id = PERSON_ID)


# Save --------------------------------------------------------------------

save(cohort_2, "pregnant_women")
