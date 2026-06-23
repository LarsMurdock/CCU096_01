#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Describe baseline characteristics of individuals in the cohort
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Load data ---------------------------------------------------------------

bl <- load_data_locally(.from = "output", .table = "baseline")
# ct <- load_data_locally(.from = "output", .table = "controls")


# Recode variables --------------------------------------------------------

# Parameters for re-coding chronic diseases
chron <- tribble(
  ~new, ~old, ~label,
  "ckd", "cov_chron_ckd_flag", "CKD",
  "htn", "cov_chron_hypertension_flag", "HTN",
  "hf", "cov_chron_heart_failure_flag", "HF",
  "livd", "cov_chron_liver_disease_flag", "LIVD",
  "ihd", "cov_chron_ihd_flag", "IHD",
  "pad", "cov_chron_pad_flag", "PAD",
  "stroke", "cov_chron_stroke_flag", "STROKE",
  "ascvd", "cov_chron_ascvd_flag", "ASCVD",
  "diab", "cov_chron_diabetes_all_flag", "DIAB",
  "diab1", "cov_chron_diabetes_type1_flag", "DIAB1",
  "diab2", "cov_chron_diabetes_type2_flag", "DIAB2",
  "dep", "cov_chron_depression_flag", "DEP",
  "asthm", "cov_chron_asthma_flag", "ASTHM",
  "copd", "cov_chron_copd_flag", "COPD",
  "asthm_copd", "cov_chron_asthma_copd_flag", "ASTHM_COPD",
  "osa", "cov_chron_sleep_apnoea_flag", "OSA",
  "canc", "cov_chron_cancer_flag", "CANC",
) |> 
  mutate(new = paste0("hist_", new))

# Function to re-code disease flags
code_flags <- function(flag, disease) {
  factor(
    flag, levels = c(1, 0), 
    labels = c(paste(disease, "yes"), paste(disease, "no"))
  )
}

# Function to code missing indicator variables
code_missing_indicator <- function(var) {
  var_name <- deparse(substitute(var))
  if_else(is.na(var), paste("missing", var_name), paste("present", var_name))
}

bl2 <- bl |>
  mutate(
    # Code year of obesity diagnosis
    study_start_year = as.numeric(format(study_start_date, "%Y")),
    #  Code time period
    time_period = case_when(
      study_start_year <  2022 ~ "1",
      study_start_year >= 2022 ~ "2",
      .default = NA
    ),
    time_period = factor(
      time_period, levels = c(1:2), labels = c("2019 to 2021", "2022 to 2025") 
    ),
    #  Code binary age
    age_binary = case_when(
      study_start_age <  40 ~ "1",
      study_start_age >= 40 ~ "2",
      .default = NA
    ),
    age_binary = factor(
      age_binary, levels = c(1:2), labels = c("<40 years", ">=40 years") 
    ),
    # Re-code chronic diseases for baseline table
    across(
      all_of(chron$old), 
      \(x) code_flags(x, chron$label[match(cur_column(), chron$old)]),
      .names = "tmp_{.col}"
    ),
    # Code missing indicator variables
    across(
      c(egfr, sbp, bmi_2yr, bmi_5yr, hba1c, tchol_hdl_ratio),
      \(x) code_missing_indicator(x), .names = "{.col}_na"
    )
  ) |>
  rename_with(\(x) chron$new, starts_with("tmp_")) |>
  select(-starts_with("tmp_"))

# ct <- ct |>
#   mutate(
#     across(
#       all_of(chron$old), 
#       \(x) code_flags(x, chron$label[match(cur_column(), chron$old)]),
#       .names = "tmp_{.col}"
#     ),
#     # Code missing indicator variables
#     across(
#       c(egfr, sbp, bmi_2yr, bmi_5yr, hba1c, tchol_hdl_ratio),
#       \(x) code_missing_indicator(x), .names = "{.col}_na"
#     )
#   ) |>
#   rename_with(\(x) chron$new, starts_with("tmp_")) |>
#   select(-starts_with("tmp_"))


# Reduce data -------------------------------------------------------------

# For describing the most and least deprived individuals
bl2_imd <- bl2 |>
  filter(str_detect(imd_quintile, "1|5")) |>
  mutate(
    imd_quintile = fct_recode(
      imd_quintile, "Least deprived" = "5 (Least deprived)", 
      "Most deprived" = "1 (Most deprived)"
    )
  )


# Summarise characteristics -----------------------------------------------

generate_baseline_table <- function() {
  
  # Function to remove things from a vector
  remove_characteristic <- function(.vector, .characteristic) {
    .vector[.vector != .characteristic]
  }
  
  # Define variables to summarise
  age_median_var <- "study_start_age"
  con_var <- c(
    "study_start_age", "bmi_2yr", "bmi_5yr", "sbp", "egfr", "hba1c", 
    "tchol_hdl_ratio"
  )
  cat_var <- c(
    "age_grp", "sex", "imd_quintile", "ethnicity", "smoking", "hist_htn", 
    "hist_hf", "hist_ckd", "hist_ascvd", "hist_diab", "hist_diab1", 
    "hist_diab2", "hist_livd", "hist_dep", "hist_osa", "hist_canc", 
    "hist_asthm_copd", "egfr_na", "sbp_na", "bmi_2yr_na", "bmi_5yr_na", 
    "hba1c_na", "tchol_hdl_ratio_na"
  )
  no_sex <- remove_characteristic(cat_var, "sex")
  no_imd <- remove_characteristic(cat_var, "imd_quintile")
  no_age <- remove_characteristic(cat_var, "age_grp")
  
  # Function to remove p value column
  remove_p <- function(.data) .data |> select(-`P-value`)
  
  # Generate baseline table columns
  
  # Median age
  message("Summarising median age...")
  med_01_all <- summarise_continuous(bl2, con_var, .normally_distributed = FALSE)
  med_02_sex <- remove_p(summarise_continuous(bl2, con_var, sex, .normally_distributed = FALSE))
  med_03_imd <- remove_p(summarise_continuous(bl2_imd, con_var, imd_quintile, .normally_distributed = FALSE))
  med_04_tmp <- remove_p(summarise_continuous(bl2, con_var, time_period, .normally_distributed = FALSE))
  med_05_age <- remove_p(summarise_continuous(bl2, con_var, age_binary, .normally_distributed = FALSE))
  
  med_lst    <- mget(ls(pattern = "^med_[0-9]{2}"))
  med     <- reduce(med_lst, left_join, by = join_by(Characteristic, Col))
  message("Completed summarising median age.")
  
  # Continuous variables
  message("Summarising continuous variables...")
  con_01_all <- summarise_continuous(bl2, con_var)
  con_02_sex <- remove_p(summarise_continuous(bl2, con_var, sex))
  con_03_imd <- remove_p(summarise_continuous(bl2_imd, con_var, imd_quintile))
  con_04_tmp <- remove_p(summarise_continuous(bl2, con_var, time_period))
  con_05_age <- remove_p(summarise_continuous(bl2, con_var, age_binary))
  # con_06_ctl <- summarise_continuous(ct, con_var)
  con_lst    <- mget(ls(pattern = "^con_[0-9]{2}"))
  con     <- reduce(con_lst, left_join, by = join_by(Characteristic, Col))
  message("Completed summarising continuous variables.")
  
  # Categorical variables
  message("Summarising categorical variables...")
  cat_01_all <- summarise_categorical(bl2, cat_var)
  cat_02_sex <- remove_p(summarise_categorical(bl2, no_sex, sex))
  cat_03_imd <- remove_p(summarise_categorical(bl2_imd, no_imd, imd_quintile))
  cat_04_tmp <- remove_p(summarise_categorical(bl2, cat_var, time_period))
  cat_05_imd <- remove_p(summarise_categorical(bl2, no_age, age_binary))
  # cat_05_ctl <- summarise_categorical(ct, cat_var)
  cat_lst <- mget(ls(pattern = "^cat_[0-9]{2}"))
  cat     <- reduce(cat_lst, left_join, by = join_by(Characteristic, Col))
  message("Completed summarising categorical variables.")

  # Separate
  select_ <- function(data, var, rename_to, title = NULL) {
    out <- filter(data, str_detect(Col, var)) |>
    mutate(
      Characteristic = if_else(
        str_detect(Characteristic, var), rename_to, Characteristic
      )
    )
    if (!missing(title)) {
      out <- out |>
        add_row(Characteristic = title, .before = 1)
    }
    return(out)
  }
  message("Re-ordering characteristics...")
  tbl_01_age_med <- select_(med, "age", "Age, years (median [IQR])")
  tbl_02_age     <- select_(con, "age", "Age, years (mean [SD])")
  tbl_03_age_grp <- select_(cat, "^age", "Age, years (groups)")
  tbl_04_sex     <- select_(cat, "^sex", "Sex")
  tbl_05_eth     <- select_(cat, "^ethnic", "Ethnicity")
  tbl_06_imd     <- select_(cat, "^imd_quint", "Socioeconomic status quintile")
  tbl_07_smok    <- select_(cat, "^smoking", "Smoking")
  tbl_08_bmi_2yr     <- select_(con, "^bmi_2yr", "Mean", "Body-mass index (2 years), kg/m2")
  tbl_09_bmi_2yr_grp <- select_(cat, "^bmi_2yr", "Missing")
  tbl_10_bmi_5yr     <- select_(con, "^bmi_5yr", "Mean", "Body-mass index (5 years), kg/m2")
  tbl_11_bmi_5yr_grp <- select_(cat, "^bmi_5yr", "Missing")
  tbl_12_sbp     <- select_(con, "^sbp", "Mean", "Systolic blood pressure, mm Hg")
  tbl_13_sbp_grp <- select_(cat, "^sbp", "Missing")
  tbl_14_hba1c     <- select_(con, "^hba1c", "Mean", "Glycated haemoglobin, %")
  tbl_15_hba1c_grp <- select_(cat, "^hba1c", "Missing")
  tbl_16_chol     <- select_(con, "^tchol", "Mean", "Cholesterol (total cholesterol to HDL ratio)")
  tbl_17_chol_grp <- select_(cat, "^tchol", "Missing")
  tbl_18_gfr     <- select_(con, "^egfr", "Mean", "Estimated GFR, mL/min per 1.73 m2")
  tbl_19_gfr_grp <- select_(cat, "^egfr", "Missing")
  tbl_20_diab    <- select_(cat, "diab$", "Diabetes", "Co-existing chronic medical conditions")
  tbl_21_diab1   <- select_(cat, "diab1$", "  Type 1")
  tbl_22_diab2   <- select_(cat, "diab2$", "  Type 2")
  tbl_23_htn     <- select_(cat, "htn$", "Hypertension")
  tbl_24_ascvd   <- select_(cat, "ascvd", "Atherosclerotic cardiovascular disease")
  tbl_25_hf      <- select_(cat, "_hf$", "Heart failure")
  tbl_26_ckd     <- select_(cat, "ckd", "Chronic kidney disease")
  tbl_27_livd    <- select_(cat, "livd", "Liver disease")
  tbl_28_canc    <- select_(cat, "_canc$", "Cancer")
  tbl_29_asth_co <- select_(cat, "_asthm_copd$", "Asthma or COPD")
  tbl_30_osa     <- select_(cat, "_osa$", "Obstructive sleep apnoea")
  tbl_31_dep     <- select_(cat, "_dep$", "Depression")

  # Combine
  message("Combining...")
  tbl_lst <- mget(ls(pattern = "^tbl_"))
  tbl <- bind_rows(tbl_lst) |> select(-Col)
    
  return(tbl)
  message("Baseline table completed.")
}

table2 <- generate_baseline_table()
write_csv(
  table2, 
  here::here("figs", "baseline_table", 
             "01_obesity_baseline_table_11aug2025.csv")
)
