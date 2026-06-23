#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Prepare baseline data
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Load data ---------------------------------------------------------------

bl <- load_data_locally(.table = "cohort")
# ct <- load_data_locally(.table = "controls")

# Obesity sources
ob <- load_data_locally(.table = "obesity_sources")


# Recode variables --------------------------------------------------------

bl_processed <- 
  prepare_baseline_variables(bl, "cov_hx_first_obesity_30kgm2_date")
# ct_processed <- prepare_baseline_variables(ct, "entry_date")


# Join --------------------------------------------------------------------

# Get relevant variables for obesity sources
ob <- ob |> 
  drop_na(cov_hx_first_obesity_30kgm2_date) |> 
  dplyr::select(
    person_id, study_start_date = cov_hx_first_obesity_30kgm2_date, 
    obesity_source = cov_hx_first_obesity_30kgm2_source
  )

# Join to the baseline data
bl_processed <- bl_processed |>
  left_join(ob, by = join_by(person_id, study_start_date)) |>
  relocate(obesity_source, .after = "person_id")


# Write data --------------------------------------------------------------

save(bl_processed, "baseline")
# save(ct_processed, "controls")
