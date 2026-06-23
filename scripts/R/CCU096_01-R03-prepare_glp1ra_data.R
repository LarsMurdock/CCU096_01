#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Prepare GLP1-RA dispensing data
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Load data ---------------------------------------------------------------

# GLP-1RA data
glp_0 <- load_data_locally(.table = "glp1ra") |> 
  rename_with(tolower) |> 
  # Code obesity flag
  mutate(
    obesity_flag = case_when(
      cov_hx_first_obesity_date <= date ~ 1,
      .default = 0
    ),
    cov_chron_diabetes_all_flag = case_when(
      if_any(matches("^cov_chron_diabetes_.*_flag$"), \(x) x == 1) ~ 1,
      .default = 0
    )
  ) |>
  select(
    person_id, dispensing_date = date, first_date, date_of_birth,
    dispensing_age, obesity_flag, cov_hx_first_obesity_date, 
    not_dispensed = notdispensedindicator, raw_name, 
    standardised_name, product_name, imd_quintile, imd_decile, sex, 
    ethnicity_5_group, prescriber_type = prescribertype,
    private_prescription = privateprescriptionindicator,
    cov_chron_diabetes_type2_flag, cov_chron_diabetes_type2_date,
    cov_chron_diabetes_all_flag, region, in_gdppr
  )

# Quality control data
qc <- load_data_locally(.table = "quality_assurance") |>
  select(PERSON_ID, `_rule_total`)


# Clean -------------------------------------------------------------------

# Remove individuals with missing date of birth or sex
glp_1 <- glp_0 |>
  drop_na(date_of_birth, sex) |>
  filter(sex %in% c("M", "F"))

# Remove individuals ages <18 years or >=80 years
glp_2 <- glp_1 |>
  filter(dispensing_age >= 18 & dispensing_age < 81)

# Remove individuals not in GDPPR
glp_3 <- glp_2 |> filter(in_gdppr == 1)

# Remove individuals with only one dispensed GLP1-RA prescription
glp_4 <- glp_3 |> 
  group_by(person_id) |>
  filter(n() >= 2) |>
  ungroup() |>
  # Re-order
  arrange(person_id, dispensing_date)

# Remove individuals outside of England
glp_5 <- glp_4 |>
  filter(!str_detect(region, "Wales|Scotland"))

# Remove dispensing before the "2018-01-01"
glp_6 <- glp_5 |>
  filter(dispensing_date >= ymd("2018-01-01"))


# Save data ---------------------------------------------------------------

save(glp_6, "glp1ra_dispensing")
