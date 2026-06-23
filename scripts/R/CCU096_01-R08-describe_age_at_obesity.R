#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Describe characteristics of individuals with obesity by age group
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Load data ---------------------------------------------------------------

# Baseline data
bl <- load_data_locally(.from = "output", .table = "baseline")
# ct <- load_data_locally(.from = "output", .table = "controls")

# General IPD data for assessing population-level differences in age
ipd <- load_data_locally(.table = "ipd_prevalence_data_30kgm2_apr_2025")


# Reduce data -------------------------------------------------------------

# Function to code missing indicator variables
code_missing_indicator <- function(var) {
  var_name <- deparse(substitute(var))
  if_else(is.na(var), paste("missing", var_name), paste("present", var_name))
}

bl <- bl |>
  dplyr::select(
    person_id, study_start_age, age_grp, bmi_5yr, sbp, hba1c, tchol_hdl_ratio, 
    egfr, sex, imd_quintile, imd_decile, ethnicity, ethnicity_18, region, 
    smoking, cov_chron_diabetes_all_flag, cov_chron_diabetes_type1_flag, 
    cov_chron_diabetes_type2_flag, cov_chron_hypertension_flag, 
    cov_chron_ascvd_flag, cov_chron_heart_failure_flag, 
    cov_chron_ckd_flag, cov_chron_liver_disease_flag,
    cov_chron_cancer_flag, cov_chron_asthma_copd_flag, 
    cov_chron_sleep_apnoea_flag, cov_chron_depression_flag
  ) |>
  mutate(
    across(
      c(egfr, sbp, bmi_5yr, hba1c, tchol_hdl_ratio),
      \(x) code_missing_indicator(x), .names = "{.col}_na"
    )
  )



# Describe mean age at obesity in different groups ------------------------

mean_sd_age_by_ethnicity_deprivation <- function(data, .sex = NULL) {
  
  if (!missing(.sex)) {
    
    data <- data |> filter(sex == .sex)
    
  }
  
  # Mean age for everyone
  mean_age_all <- data |>
    summarise(
      mean_age = mean(study_start_age), sd = sd(study_start_age),
      median_age = median(study_start_age),
      lcl = quantile(study_start_age, 0.25), 
      ucl = quantile(study_start_age, 0.75),
      n = rnd(n_distinct(person_id) / 5, 0) * 5
    ) |>
    mutate(imd_quintile = "All", ethnicity = "All", .before = 1)
  
  # Mean age by IMD quintiles for all ethnicities
  mean_age_imd5 <- data |>
    group_by(imd_quintile) |> 
    summarise(
      mean_age = mean(study_start_age), sd = sd(study_start_age),
      median_age = median(study_start_age),
      lcl = quantile(study_start_age, 0.25), 
      ucl = quantile(study_start_age, 0.75),
      n = rnd(n_distinct(person_id) / 5, 0) * 5, .groups = "drop"
    ) |>
    mutate(ethnicity = "All", .after = "imd_quintile")
  
  # Mean age for all 5 ethnicities
  mean_age_eth5 <- data |>
    group_by(ethnicity) |> 
    summarise(
      mean_age = mean(study_start_age), sd = sd(study_start_age),
      median_age = median(study_start_age),
      lcl = quantile(study_start_age, 0.25), 
      ucl = quantile(study_start_age, 0.75),
      n = rnd(n_distinct(person_id) / 5, 0) * 5
    ) |>
    ungroup() |>
    mutate(imd_quintile = "All", .before = "ethnicity")
  
  # Mean age for all 18 ethnicities
  mean_age_eth18 <- data |>
    group_by(ethnicity_18) |> 
    summarise(
      mean_age = mean(study_start_age), sd = sd(study_start_age),
      median_age = median(study_start_age),
      lcl = quantile(study_start_age, 0.25), 
      ucl = quantile(study_start_age, 0.75),
      n = rnd(n_distinct(person_id) / 5, 0) * 5, .groups = "drop"
    ) |>
    rename(ethnicity = ethnicity_18) |>
    mutate(imd_quintile = "All", .before = "ethnicity")
  
  # Mean age by IMD quintiles and 5 ethnicity groups
  mean_age_imd5_eth5 <- data |>
    group_by(ethnicity, imd_quintile) |> 
    summarise(
      mean_age = mean(study_start_age), sd = sd(study_start_age),
      median_age = median(study_start_age),
      lcl = quantile(study_start_age, 0.25), 
      ucl = quantile(study_start_age, 0.75),
      n = rnd(n_distinct(person_id) / 5, 0) * 5, .groups = "drop"
    )
  
  # Mean age by IMD quintiles and 18 ethnicity groups
  mean_age_imd5_eth18 <- data |>
    group_by(ethnicity_18, imd_quintile) |> 
    summarise(
      mean_age = mean(study_start_age), sd = sd(study_start_age),
      median_age = median(study_start_age),
      lcl = quantile(study_start_age, 0.25), 
      ucl = quantile(study_start_age, 0.75),
      n = rnd(n_distinct(person_id) / 5, 0) * 5, .groups = "drop"
    ) |>
    rename(ethnicity = ethnicity_18)
  
  bind_rows(
    mean_age_all, mean_age_imd5, mean_age_eth5, mean_age_eth18, 
    mean_age_imd5_eth5, mean_age_imd5_eth18
  ) |>
    mutate(sex = .sex)
}

mean_age_all     <- mean_sd_age_by_ethnicity_deprivation(bl) |> mutate(sex = "All")
mean_age_females <- mean_sd_age_by_ethnicity_deprivation(bl, "Female")
mean_age_males   <- mean_sd_age_by_ethnicity_deprivation(bl, "Male")

bind_rows(mean_age_all, mean_age_females, mean_age_males) |>
  write_csv(
    here::here(
      "figs", "age_by_ethnicity_deprivation", 
      "mean_age_obesity_by_ethnicity_deprivation_11aug2025.csv"
    )
  )

# IPD



# Summarise characteristics -----------------------------------------------

summarise_by_age_bands <- function(data, 
                                   by = NULL,
                                   stratum_1 = NULL,
                                   stratum_2 = NULL) {
  
  if (!missing(by) && by %in% c("sex", "imd_decile")) {
    
    data <- data |> filter(str_detect(!!sym(by), stratum_1))
    
  } else if (!missing(by) && by == "sex and imd_decile") {
    
    # This is sadly hard coded (sex stratum must go first...)
    data <- data |> 
      filter(
        str_detect(sex, stratum_1), 
        str_detect(imd_decile, stratum_2)
      )
    
  }
  
  split_data <- data %>%
    split(.$age_grp)
  
  # Continuous variables
  con <- split_data |>
    map(
      \(.x) summarise_continuous(
        .x, c(bmi_5yr, sbp, hba1c, tchol_hdl_ratio, egfr)
      )
    ) %>%
    map2(
      .y = names(.), \(.x, .y)
      rename_with(
        .x, \(x) case_when(str_detect(x, "N =") ~ paste0(.y, x), .default = x)
      )
    ) |>
    reduce(left_join, by = join_by(Col, Characteristic))
  
  # Categorical variables
  cat <- split_data |>
    map(
      \(.x) summarise_categorical(
        .x, c(
          sex, imd_quintile, imd_decile, ethnicity, ethnicity_18, region,
          smoking, cov_chron_diabetes_all_flag, cov_chron_diabetes_type1_flag, 
          cov_chron_diabetes_type2_flag, cov_chron_hypertension_flag, 
          cov_chron_ascvd_flag, cov_chron_heart_failure_flag, 
          cov_chron_ckd_flag, cov_chron_liver_disease_flag,
          cov_chron_cancer_flag, cov_chron_asthma_copd_flag, 
          cov_chron_sleep_apnoea_flag, cov_chron_depression_flag,
          egfr_na, sbp_na, bmi_5yr_na, hba1c_na, tchol_hdl_ratio_na
        )
      )
    ) %>%
    map2(
      .y = names(.), \(.x, .y)
      rename_with(
        .x, \(x) case_when(str_detect(x, "N =") ~ paste0(.y, x), .default = x)
      )
    ) |>
    reduce(left_join, by = join_by(Col, Characteristic)) |>
    mutate(
      across(everything(), \(x) if_else(str_detect(x, "^5 "), "REDACTED", x))
    )

  return(bind_rows(con, cat))
}

# All sexes, all deprivation
case_overall <- summarise_by_age_bands(bl)
# ctrl_overall <- summarise_by_age_bands(ct)

# Male sex, all deprivation
case_male <- summarise_by_age_bands(bl, by = "sex", stratum_1 = "Male")
# ctrl_male <- summarise_by_age_bands(ct, by = "sex", stratum_1 = "Male")

# Female sex, all deprivation
case_female <- summarise_by_age_bands(bl, by = "sex", stratum_1 = "Female")
# ctrl_female <- summarise_by_age_bands(ct, by = "sex", stratum_1 = "Female")

# All sexes, highest deprivation
dep1_case <- summarise_by_age_bands(bl, by = "imd_decile", stratum_1 = "Most")
# dep1_ctrl <- summarise_by_age_bands(ct, by = "imd_decile", stratum_1 = "Most")

# Male sex, highest deprivation
dep1_case_male <- summarise_by_age_bands(bl, by = "sex and imd_decile", stratum_1 = "Male", stratum_2 = "Most")
# dep1_ctrl_male <- summarise_by_age_bands(ct, by = "sex and imd_decile", stratum_1 = "Male", stratum_2 = "Most")

# Female sex, highest deprivation
dep1_case_female <- summarise_by_age_bands(bl, by = "sex and imd_decile", stratum_1 = "Female", stratum_2 = "Most")
# dep1_ctrl_female <- summarise_by_age_bands(ct, by = "sex and imd_decile", stratum_1 = "Female", stratum_2 = "Most")

# All sexes, lowest deprivation
dep10_case <- summarise_by_age_bands(bl, by = "imd_decile", stratum_1 = "Least")
# dep10_ctrl <- summarise_by_age_bands(ct, by = "imd_decile", stratum_1 = "Least")

# Male sex, lowest deprivation
dep10_case_male <- summarise_by_age_bands(bl, by = "sex and imd_decile", stratum_1 = "Male", stratum_2 = "Least")
# dep10_ctrl_male <- summarise_by_age_bands(ct, by = "sex and imd_decile", stratum_1 = "Male", stratum_2 = "Least")

# Female sex, lowest deprivation
dep10_case_female <- summarise_by_age_bands(bl, by = "sex and imd_decile", stratum_1 = "Female", stratum_2 = "Least")
# dep10_ctrl_female <- summarise_by_age_bands(ct, by = "sex and imd_decile", stratum_1 = "Female", stratum_2 = "Least")


# Output ------------------------------------------------------------------

# All sexes, all deprivation
write_csv(case_overall, here::here("figs", "characteristics_by_age", "01_overall_obesity_characteristics_by_age.csv"))
# write_csv(ctrl_overall, here::here("figs", "characteristics_by_age", "02_overall_controls_characteristics_by_age.csv"))

# Male sex, all deprivation
write_csv(case_male, here::here("figs", "characteristics_by_age", "03_males_obesity_characteristics_by_age.csv"))
# write_csv(ctrl_male, here::here("figs", "characteristics_by_age", "04_males_controls_characteristics_by_age.csv"))

# Female sex, all deprivation
write_csv(case_female, here::here("figs", "characteristics_by_age", "05_females_obesity_characteristics_by_age.csv"))
# write_csv(ctrl_female, here::here("figs", "characteristics_by_age", "06_females_controls_characteristics_by_age.csv"))

# ---

# All sexes, highest deprivation
write_csv(dep1_case, here::here("figs", "characteristics_by_age", "07_dep1_overall_obesity_characteristics_by_age.csv"))
# write_csv(dep1_ctrl, here::here("figs", "characteristics_by_age", "08_dep1_overall_controls_characteristics_by_age.csv"))

# Male sex, highest deprivation
write_csv(dep1_case_male, here::here("figs", "characteristics_by_age", "09_dep1_males_obesity_characteristics_by_age.csv"))
# write_csv(dep1_ctrl_male, here::here("figs", "characteristics_by_age", "10_dep1_males_controls_characteristics_by_age.csv"))

# Female sex, highest deprivation
write_csv(dep1_case_female, here::here("figs", "characteristics_by_age", "11_dep1_females_obesity_characteristics_by_age.csv"))
# write_csv(dep1_ctrl_female, here::here("figs", "characteristics_by_age", "12_dep1_females_controls_characteristics_by_age.csv"))

# ---

# All sexes, lowest deprivation
write_csv(dep10_case, here::here("figs", "characteristics_by_age", "13_dep10_overall_obesity_characteristics_by_age.csv"))
# write_csv(dep10_ctrl, here::here("figs", "characteristics_by_age", "14_dep10_overall_controls_characteristics_by_age.csv"))

# Male sex, lowest deprivation
write_csv(dep10_case_male, here::here("figs", "characteristics_by_age", "15_dep10_males_obesity_characteristics_by_age.csv"))
# write_csv(dep10_ctrl_male, here::here("figs", "characteristics_by_age", "16_dep10_males_controls_characteristics_by_age.csv"))

# Female sex, lowest deprivation
write_csv(dep10_case_female, here::here("figs", "characteristics_by_age", "17_dep10_females_obesity_characteristics_by_age.csv"))
# write_csv(dep10_ctrl_female, here::here("figs", "characteristics_by_age", "18_dep10_females_controls_characteristics_by_age.csv"))
