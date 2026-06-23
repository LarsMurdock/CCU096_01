#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Describe GLP1-RA dispensing
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Variables ---------------------------------------------------------------

# New figure sub-directory
dir.create(here::here("figs", "therapy_trends"), showWarnings = FALSE)


# Load data ---------------------------------------------------------------

# GLP-1RA data
glp <- load_data_locally(.from = "output", .table = "glp1ra_dispensing")

# Bariatric surgery data
bs  <- load_data_locally(.from = "output", .table = "bariatric_surgery")


# Filter for obese individuals --------------------------------------------

glp_obesity <- glp |> filter(obesity_flag == 1)


# Overall summary ---------------------------------------------------------

glp_overall <- 
  summarise_medications(
    glp,
    by_vars = c("obesity_flag", "dispensing_date"),
    pivot_vars = list(id = "dispensing_date", names = "obesity_flag"),
    rename_map = c(`1` = "n_obese", `0` = "n_non_obese")
  ) |>
  mutate(n_total = rowSums(across(c(n_obese, n_non_obese)), na.rm = TRUE))


# Drug type and bariatric surgery -----------------------------------------

# Overall
glp_overall_drug_type_month <- glp |> 
  mutate(dispensing_date = floor_date(dispensing_date, "month")) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, dispensing_date)
  ) |>
  pivot_wider(
    id_cols = dispensing_date, names_from = product_name, values_from = n_total
  ) |>
  rename(Other = Missing) |> 
  arrange(dispensing_date) |>
  left_join(
    bs |>
      mutate(dispensing_date = floor_date(date, "month")) |>
      summarise(
        bariatric_surgery = n_distinct(person_id), .by = dispensing_date
      ) |>
      arrange(dispensing_date), 
    by = join_by(dispensing_date)
  ) |>
  mutate(
    Total = rowSums(
      across(c(Rybelsus, Other, Ozempic, Mounjaro, Wegovy)), na.rm = TRUE
    ),
    across(!matches("dispensing_date"), \(x) rnd(x / 5, 0) * 5),
    across(
      !matches("dispensing_date"), 
      \(x) if_else(x == 5, "REDACTED", as.character(x))
    )
  ) |>
  write_csv(here::here("figs", "therapy_trends", "glp1ra_bariatric_surgery_overall.csv"))

# By sex
glp_overall_drug_type_sex_month <- glp |> 
  mutate(dispensing_date = floor_date(dispensing_date, "month")) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, sex, dispensing_date)
  ) |>
  pivot_wider(
    id_cols = c(sex, dispensing_date), names_from = product_name,
    values_from = n_total
  ) |>
  rename(Other = Missing) |> 
  arrange(sex, dispensing_date) |>
  left_join(
    bs |>
      mutate(dispensing_date = floor_date(date, "month")) |>
      summarise(
        bariatric_surgery = n_distinct(person_id), .by = c(sex, dispensing_date)
      ) |>
      arrange(dispensing_date), 
    by = join_by(sex, dispensing_date)
  ) |>
  mutate(
    Total = rowSums(
      across(c(Rybelsus, Other, Ozempic, Mounjaro, Wegovy)), na.rm = TRUE
    ),
    across(!matches("dispensing_date|sex"), \(x) rnd(x / 5, 0) * 5),
    across(
      !matches("dispensing_date"), 
      \(x) if_else(x == 5, "REDACTED", as.character(x))
    )
  ) |>
  write_csv(here::here("figs", "therapy_trends", "glp1ra_bariatric_surgery_sex.csv"))



# Summary by sex ----------------------------------------------------------

glp_sex <- 
  summarise_medications(
    glp,
    by_vars = c("obesity_flag", "sex", "dispensing_date"),
    pivot_vars = list(id = "dispensing_date", names = c("obesity_flag", "sex")),
    rename_map = c(
      `1_M` = "n_obese_males", `0_M` = "n_non_obese_males",
      `1_F` = "n_obese_females", `0_F` = "n_non_obese_females"
    )
  ) |>
  mutate(
    n_total_males = rowSums(across(ends_with("_males")), na.rm = TRUE),
    n_total_females = rowSums(across(ends_with("females")), na.rm = TRUE)
  )


# Summary by IMD quintile -------------------------------------------------

glp_imd <- 
  summarise_medications(
    glp,
    by_vars = c("obesity_flag", "imd_quintile", "dispensing_date"),
    pivot_vars = list(
      id = "dispensing_date", names = c("obesity_flag", "imd_quintile")
    ),
    rename_map = set_names(
      paste0("n_", c("non_obese_", "obese_"), rep(1:5, each = 2)),
      paste0(c("0_", "1_"), rep(1:5, each = 2))
    )
  ) |> 
  mutate(
    n_total_1 = rowSums(across(ends_with("_1")), na.rm = TRUE),
    n_total_2 = rowSums(across(ends_with("_2")), na.rm = TRUE),
    n_total_3 = rowSums(across(ends_with("_3")), na.rm = TRUE),
    n_total_4 = rowSums(across(ends_with("_4")), na.rm = TRUE),
    n_total_5 = rowSums(across(ends_with("_5")), na.rm = TRUE)
  )


# Summarise drug types in obese people ------------------------------------

obese_by_drug <- glp |> 
  filter(obesity_flag == 1) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, dispensing_date)
  ) |>
  pivot_wider(
    id_cols = dispensing_date, names_from = product_name, values_from = n_total
  ) |>
  rename(Other = Missing) |> 
  arrange(dispensing_date) |>
  mutate(
    Total = rowSums(
      across(c(Rybelsus, Other, Ozempic, Mounjaro, Wegovy)), na.rm = TRUE
    ),
    across(!matches("dispensing_date"), \(x) rnd(x / 5, 0) * 5),
    across(!matches("dispensing_date"), \(x) if_else(x == 5, "REDACTED", as.character(x)))
  )

write_csv(obese_by_drug, here::here("figs", "therapy_trends", "glp1ra_obesity_drug_type.csv"))

obese_by_drug_diabetes <- glp |> 
  filter(obesity_flag == 1) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, cov_chron_diabetes_type2_flag, dispensing_date)
  ) |>
  arrange(dispensing_date, cov_chron_diabetes_type2_flag, product_name) |>
  pivot_wider(
    id_cols = dispensing_date, 
    names_from = c(cov_chron_diabetes_type2_flag, product_name), 
    values_from = n_total
  ) |>
  rename_with(
    \(x) case_when(
      str_detect(x, "^0_") ~ str_replace(x, "^0_", "n_non_diab_"),
      str_detect(x, "^1_") ~ str_replace(x, "^1_", "n_diab_"),
      str_detect(x, "_Missing$") ~ str_replace(x, "_Missing$", "_Other"),
      .default = x
    )
  ) |>
  rename_with(
    \(x) case_when(
      str_detect(x, "_Missing$") ~ str_replace(x, "_Missing$", "_Other"),
      .default = x
    )
  ) |>
  mutate(
    across(!matches("dispensing_date"), \(x) rnd(x / 5, 0) * 5),
    across(!matches("dispensing_date"), \(x) if_else(x == 5, "REDACTED", as.character(x)))
  )

write_csv(obese_by_drug_diabetes, here::here("figs", "therapy_trends", "glp1ra_obesity_diabetes_drug_type.csv"))


# Summarise off-label -----------------------------------------------------

off_label <- glp |> 
  filter(obesity_flag == 1 & cov_chron_diabetes_type2_flag == 0) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, dispensing_date)
  ) |>
  pivot_wider(
    id_cols = dispensing_date, names_from = product_name, 
    values_from = n_total
  ) |>
  rename(Other = Missing) |> 
  arrange(dispensing_date) |>
  mutate(
    Total = rowSums(
      across(c(Rybelsus, Other, Ozempic, Mounjaro, Wegovy)), na.rm = TRUE
    )
  ) |>
  mutate(
    across(!matches("dispensing_date"), \(x) rnd(x / 5, 0) * 5),
    across(!matches("dispensing_date"), \(x) if_else(x == 5, "REDACTED", as.character(x)))
  )

write_csv(off_label, here::here("figs", "therapy_trends", "glp1ra_off_label_drug_type.csv"))

off_label_sex <- glp |> 
  filter(obesity_flag == 1 & cov_chron_diabetes_type2_flag == 0) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, sex, dispensing_date)
  ) |>
  pivot_wider(
    id_cols = dispensing_date, names_from = c(sex, product_name), 
    values_from = n_total
  ) |>
  arrange(dispensing_date) |>
  rename_with(
    \(x) case_when(
      str_detect(x, "_Missing$") ~ str_replace(x, "_Missing$", "_Other"),
      .default = x
    )
  ) |>
  mutate(
    across(!matches("dispensing_date"), \(x) rnd(x / 5, 0) * 5),
    across(!matches("dispensing_date"), \(x) if_else(x == 5, "REDACTED", as.character(x)))
  )

write_csv(off_label_sex, here::here("figs", "therapy_trends", "glp1ra_off_label_sex_drug_type.csv"))

off_label_imd <- glp |> 
  filter(obesity_flag == 1 & cov_chron_diabetes_all_flag == 0) |>
  summarise(
    n_total = n_distinct(person_id), .by = c(product_name, imd_quintile, dispensing_date)
  ) |>
  pivot_wider(
    id_cols = dispensing_date, names_from = c(imd_quintile, product_name), 
    values_from = n_total
  ) |>
  arrange(dispensing_date) |>
  rename_with(
    \(x) case_when(
      str_detect(x, "_Missing$") ~ str_replace(x, "_Missing$", "_Other"),
      .default = x
    )
  ) |>
  mutate(
    across(!matches("dispensing_date"), \(x) rnd(x / 5, 0) * 5),
    across(!matches("dispensing_date"), \(x) if_else(x == 5, "REDACTED", as.character(x)))
  )

write_csv(off_label_imd, here::here("figs", "therapy_trends", "glp1ra_off_label_imd_drug_type.csv"))


# Write data --------------------------------------------------------------

write_csv(glp_overall, here::here("figs", "therapy_trends", "glp1ra_dispensing_total.csv"))
write_csv(glp_sex, here::here("figs", "therapy_trends", "glp1ra_dispensing_sex.csv"))
write_csv(glp_imd, here::here("figs", "therapy_trends", "glp1ra_dispensing_imd.csv"))

