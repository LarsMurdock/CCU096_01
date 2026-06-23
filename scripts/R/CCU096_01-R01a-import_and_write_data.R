#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Import data from database and write to project folder
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))
library(DBI)


# Define helper functions -------------------------------------------------

exports <- function(x) {
  paste0("export_", x)
}


# Define variables --------------------------------------------------------

# Tables
# tbl <- exports(
#   c(
#     # cohort_tables,
#     # ipd,
#     # incidence_tables,
#     # prevalence_tables
#   )
# )
tbl <- c(
#   "tmp_incidence_numerator_comb1_18grp_30kgm2_std",
#   "tmp_incidence_denominator_comb1_18grp_30kgm2_std",
#   "export_prevalence_comb1_30kgm2_std"
#  "export_prevalence_comb1_18grp_30kgm2_std"
  
  # "export_prevalence_comb1_30kgm2_182_std",
  # "export_prevalence_comb1_30kgm2_365_std",
  # "export_prevalence_comb1_30kgm2_730_std",
  # "export_prevalence_comb1_30kgm2_1095_std",
  # "export_prevalence_comb1_30kgm2_1460_std",
  # "export_prevalence_comb1_30kgm2_1825_std"
 
  # "tmp_incidence_numerator_comb1_30kgm2_182_std",
  # "tmp_incidence_numerator_comb1_30kgm2_365_std",
  # "tmp_incidence_numerator_comb1_30kgm2_730_std", 
  # "tmp_incidence_numerator_comb1_30kgm2_1095_std", 
  # "tmp_incidence_numerator_comb1_30kgm2_1460_std",
  # "tmp_incidence_numerator_comb1_30kgm2_1825_std",
  # 
  # "tmp_incidence_denominator_comb1_30kgm2_182_std",
  # "tmp_incidence_denominator_comb1_30kgm2_365_std",
  # "tmp_incidence_denominator_comb1_30kgm2_730_std",
  # "tmp_incidence_denominator_comb1_30kgm2_1095_std",
  # "tmp_incidence_denominator_comb1_30kgm2_1460_std",
  # "tmp_incidence_denominator_comb1_30kgm2_1825_std",
  
  "export_missing_data_summary"
)

# New sub-directory in the data folder
dir.create(here::here("data"), showWarnings = FALSE)


# Establish Databricks connection -----------------------------------------

con <- dbConnect(
  odbc::odbc(), dsn = data_source, HTTPPath = path, PWD = password
)


# Import data and write to project folder ---------------------------------

for (t in tbl) {
  import_data(
    .connection = con, .db_name = db, .proj = prj, .table_name = t, 
    .start_date = stringr::str_remove_all(start_date, "-")
  )
}
