#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Prepare prevalence data
#
#*******************************************************************************


# Libraries and functions -------------------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Variables ---------------------------------------------------------------

# New figure sub-directory
dir.create(here::here("figs", "prevalence_tables"), showWarnings = FALSE)


# Load data ---------------------------------------------------------------

# Define all files in the data folder
dt <- list.files(here::here("data"))
# Define prevalence tables
prevalence_tables <- gsub("(^ccu096_01_\\d+_|\\.rds)", "", dt[grep("prev", dt)])
# Remove IPD data
prevalence_tables <- prevalence_tables[!grepl("ipd", prevalence_tables)]
# Remove MSOA tables as we want to load them separately as they are in a 
# different format
prevalence_no_msoa <- prevalence_tables[!grepl("msoa", prevalence_tables)]

for (i in 1:length(prevalence_no_msoa)) {
  assign(
    gsub("prevalence_", "", prevalence_no_msoa[i]),
    load_data_locally(.table = prevalence_no_msoa[i]) |>
    mutate(
      across(matches("numerat|denominat"), \(x) as.integer(as.character(x)))
    ) |>
    arrange(month_date) |>
    filter(month_date >= ymd(start_date)) |>
    mutate(row = row_number()),
    envir = .GlobalEnv
  )
}

# Load MSOA tables separately (crude)
msoa <- load_data_locally(.table = "prevalence_msoa") |>
  mutate(across(matches("numerat|denominat"), \(x) as.integer(as.character(x))))
msoa_30kgm2 <- load_data_locally(.table = "prevalence_msoa_30kgm2") |>
  mutate(across(matches("numerat|denominat"), \(x) as.integer(as.character(x))))

msoa_std <- load_data_locally(.table = "prevalence_msoa_std") |>
  mutate(across(matches("numerat|denominat"), \(x) as.integer(as.character(x))))
msoa_30kgm2_std <- load_data_locally(.table = "prevalence_msoa_30kgm2_std") |>
  mutate(across(matches("numerat|denominat"), \(x) as.integer(as.character(x))))


# Prepare -----------------------------------------------------------------

# all |>
#   mutate(across(matches("numera|denomi"), \(.x) rnd(.x / 5, 0) * 5)) |>
#   write_csv(here::here("figs", "prevalence_tables", "obesity_prevalence_ethnicity_specific_all.csv"))
# 
# all_30kgm2 |>
#   mutate(across(matches("numera|denomi"), \(.x) rnd(.x / 5, 0) * 5)) |>
#   write_csv(here::here("figs", "prevalence_tables", "obesity_prevalence_30kgm2_all.csv"))
# 
# sex |> 
#   pivot_wider(
#     id_cols = month_date, names_from = sex, 
#     values_from = c(numerator, denominator, prevalence)
#   ) |>
#   mutate(across(matches("numera|denomi"), \(.x) rnd(.x / 5, 0) * 5)) |>
#   write_csv(here::here("figs", "prevalence_tables", "obesity_prevalence_ethnicity_specific_sex.csv"))
# 
# eth18 |> 
#   pivot_wider(
#     id_cols = month_date, names_from = ethnicity_18_recoded , 
#     values_from = c(numerator, denominator, prevalence)
#   ) |>
#   mutate(across(matches("numera|denomi"), \(.x) rnd(.x / 5, 0) * 5)) |>
#   write_csv(here::here("figs", "prevalence_tables", "obesity_prevalence_ethnicity_specific_eth18.csv"))
# 
# imd |> 
#   pivot_wider(
#     id_cols = month_date, names_from = imd_quintile, 
#     values_from = c(numerator, denominator, prevalence)
#   ) |>
#   mutate(across(matches("numera|denomi"), \(.x) rnd(.x / 5, 0) * 5)) |>
#   write_csv(here::here("figs", "prevalence_tables", "obesity_prevalence_ethnicity_specific_imd.csv"))
# 
# imd10 |>
#   pivot_wider(
#     id_cols = month_date, names_from = imd_decile, 
#     values_from = c(numerator, denominator, prevalence)
#   ) |>
#   mutate(across(matches("numera|denomi"), \(.x) rnd(.x / 5, 0) * 5)) |>
#   write_csv(here::here("figs", "prevalence_tables", "obesity_prevalence_ethnicity_specific_imd10.csv"))


# Prepare standardised prevalence data ------------------------------------

prepare_prev_std <- function(data, 
                             std = european_standard_popn_2013) {
  
  data |>
    rename(age = age_grp_5yr) |>
    left_join(
      std |>
        mutate(
          age = str_remove(AgeGroup, " years"),
          age = str_replace(age, "plus", "+"),
          sex = str_extract(Sex, "^[M|F]")
        ),
      by = join_by(age, sex)
    ) |>
    rename(date_month = month_date) |>
    arrange(date_month, age, sex) |> 
    filter(date_month <= ymd("2025-05-01")) |>  
    mutate(date_quarter = floor_date(date_month, "quarter"))
  
}

all_std_2   <- prepare_prev_std(all_std)
eth5_std_2  <- prepare_prev_std(eth5_std)
eth18_std_2 <- prepare_prev_std(eth18_std)
imd_std_2   <- prepare_prev_std(imd_std)
imd10_std_2 <- prepare_prev_std(imd10_std)
reg_std_2   <- prepare_prev_std(reg_std)
lacl_std_2  <- prepare_prev_std(lacl_std)
gdp_std_2   <- prepare_prev_std(gdp_std)

all_30kgm2_std_2 <- prepare_prev_std(all_30kgm2_std)
eth5_30kgm2_std_2  <- prepare_prev_std(eth5_30kgm2_std)
eth18_30kgm2_std_2 <- prepare_prev_std(eth18_30kgm2_std)
imd_30kgm2_std_2   <- prepare_prev_std(imd_30kgm2_std)
imd10_30kgm2_std_2 <- prepare_prev_std(imd10_30kgm2_std)
reg_30kgm2_std_2   <- prepare_prev_std(reg_30kgm2_std)
lacl_30kgm2_std_2 <- prepare_prev_std(lacl_30kgm2_std)
gdp_30kgm2_std_2  <- prepare_prev_std(gdp_30kgm2_std)

# MSOA
msoa_std_2019_2 <- msoa_std |>
  filter(month_date == ymd("2019-11-01")) |>
  prepare_prev_std()
msoa_30kgm2_std_2019_2 <- msoa_30kgm2_std |>
  filter(month_date == ymd("2019-11-01")) |>
  prepare_prev_std()

msoa_std_2 <- msoa_std |>
  filter(month_date == ymd("2025-04-01")) |>
  prepare_prev_std()
msoa_30kgm2_std_2 <- msoa_30kgm2_std |>
  filter(month_date == ymd("2025-04-01")) |>
  prepare_prev_std()


# Calculate standardised prevalence (for each sex) ------------------------

by = c("date_month", "sex")
calculate_prevalence_std <- function(data, by_cols) {
  
  if ("age" %in% by_cols) {
    
    data <- data |>
      mutate(
        age = case_when(
          str_detect(age, "^1") ~ "18-19",
          str_detect(age, "^2") ~ "20-29",
          str_detect(age, "^3") ~ "30-39",
          str_detect(age, "^4") ~ "40-49",
          str_detect(age, "^5") ~ "50-59",
          str_detect(age, "^6") ~ "60-69",
          str_detect(age, "^7") ~ "70-79",
          str_detect(age, "^8") ~ "80-89",
          str_detect(age, "^9") ~ "90+",
        )
      )
    
  }
  
  # By month
  std_month <- data |>
    mutate(
      numerator = if_else(is.na(numerator), 0, numerator),
      prevalence = if_else(is.na(prevalence), 0, prevalence),
      weighted_prevalence = prevalence * EuropeanStandardPopulation
    ) |>
    summarise(
      prevalence_std = sum(weighted_prevalence) / sum(EuropeanStandardPopulation),
      prevalence_crude = sum(numerator) / sum(denominator) * 100,
      numerator_rnd = rnd(sum(numerator) / 5, 0) * 5, 
      denominator_rnd = rnd(sum(denominator) / 5, 0) * 5,
      .by = all_of(by_cols)
    )
  
  return(std_month)
}

std_all   <- calculate_prevalence_std(all_std_2, by)
std_age   <- calculate_prevalence_std(all_std_2, c(by, "age"))
std_eth5  <- calculate_prevalence_std(eth5_std_2, c(by, "ethnicity_5_recoded"))
std_eth18 <- calculate_prevalence_std(eth18_std_2, c(by, "ethnicity_18_recoded"))
std_imd   <- calculate_prevalence_std(imd_std_2, c(by, "imd_quintile"))
std_imd10 <- calculate_prevalence_std(imd10_std_2, c(by, "imd_decile"))
std_reg   <- calculate_prevalence_std(reg_std_2, c(by, "region_recoded"))
std_lacl  <- calculate_prevalence_std(lacl_std_2, c(by, "LACLASS"))
std_gdp   <- calculate_prevalence_std(gdp_std_2, c(by, "gdp_grp"))

std_all_30kgm2   <- calculate_prevalence_std(all_30kgm2_std_2, by)
std_age_30kgm2   <- calculate_prevalence_std(all_30kgm2_std_2, c(by, "age"))
std_eth5_30kgm2  <- calculate_prevalence_std(eth5_30kgm2_std_2, c(by, "ethnicity_5_recoded"))
std_eth18_30kgm2 <- calculate_prevalence_std(eth18_30kgm2_std_2, c(by, "ethnicity_18_recoded"))
std_imd_30kgm2   <- calculate_prevalence_std(imd_30kgm2_std_2, c(by, "imd_quintile"))
std_imd10_30kgm2 <- calculate_prevalence_std(imd10_30kgm2_std_2, c(by, "imd_decile"))
std_reg_30kgm2   <- calculate_prevalence_std(reg_30kgm2_std_2, c(by, "region_recoded"))
std_lacl_30kgm2  <- calculate_prevalence_std(lacl_30kgm2_std_2, c(by, "LACLASS"))
std_gdp_30kgm2  <- calculate_prevalence_std(gdp_30kgm2_std_2, c(by, "gdp_grp"))

# MSOA (ethnicity-specific)
std_msoa   <- calculate_prevalence_std(msoa_std_2, c(by, "msoa"))
std_msoa_f <- std_msoa |> filter(sex == "F")
std_msoa_m <- std_msoa |> filter(sex == "M")

# MSOA (30kg/m2 obesity threshold)
std_msoa_30kgm2   <- calculate_prevalence_std(msoa_30kgm2_std_2, c(by, "msoa"))
std_msoa_30kgm2_f <- std_msoa_30kgm2 |> filter(sex == "F")
std_msoa_30kgm2_m <- std_msoa_30kgm2 |> filter(sex == "M")


# Calculate standardised prevalence (overall) -----------------------------

by_overall <- "date_month"

std_all_overall   <- calculate_prevalence_std(all_std_2, by_overall)

std_all_30kgm2_overall  <- calculate_prevalence_std(all_30kgm2_std_2, by_overall) 
std_age_30kgm2_overall  <- calculate_prevalence_std(all_30kgm2_std_2, c(by_overall, "age"))
std_eth5_30kgm2_overall <- calculate_prevalence_std(eth5_30kgm2_std_2, c(by_overall, "ethnicity_5_recoded"))
std_imd_30kgm2_overall  <- calculate_prevalence_std(imd_30kgm2_std_2, c(by_overall, "imd_quintile"))
std_lacl_30kgm2_overall <- calculate_prevalence_std(lacl_30kgm2_std_2, c(by_overall, "LACLASS"))
std_gdp_30kgm2_overall  <- calculate_prevalence_std(gdp_30kgm2_std_2, c(by_overall, "gdp_grp"))

# MSOA
std_msoa_2019_overall <- 
  calculate_prevalence_std(msoa_std_2019_2, c(by_overall, "msoa")) 
std_msoa_30kgm2_2019_overall <- 
  calculate_prevalence_std(msoa_30kgm2_std_2019_2, c(by_overall, "msoa")) 

std_msoa_overall <- 
  calculate_prevalence_std(msoa_std_2, c(by_overall, "msoa")) 
std_msoa_30kgm2_overall <- 
  calculate_prevalence_std(msoa_30kgm2_std_2, c(by_overall, "msoa")) 


# Save --------------------------------------------------------------------

write_chunk <- function(data, file_name) {
  
  if (nrow(data) > 1000 && nrow(data) <= 2000) {
    
    data_1 <- data |> slice(1:1000)
    data_2 <- data |> slice(1001:2000)
    
    write_csv(
      data_1, 
      here::here("figs", "prevalence_tables", paste0(file_name, "_chunk_01.csv"))
    )
    
    write_csv(
      data_2, 
      here::here("figs", "prevalence_tables", paste0(file_name, "_chunk_02.csv"))
    )
    
  } else if (nrow(data) > 2000 && nrow(data) <= 3000) {
    
    data_1 <- data |> slice(1:1000)
    data_2 <- data |> slice(1001:2000)
    data_3 <- data |> slice(2001:3000)
    
    write_csv(
      data_1, 
      here::here("figs", "prevalence_tables", paste0(file_name, "_chunk_01.csv"))
    )
    
    write_csv(
      data_2, 
      here::here("figs", "prevalence_tables", paste0(file_name, "_chunk_02.csv"))
    )
    
    write_csv(
      data_3, 
      here::here("figs", "prevalence_tables", paste0(file_name, "_chunk_03.csv"))
    )
    
  } else if (nrow(data) <= 1000) {
    
    write_csv(
      data, here::here("figs", "prevalence_tables", paste0(file_name, ".csv"))
    )
    
  } else {
    
    message("You need more segments.")
    
  }
  
}

write_chunk(std_all, "obesity_prevalence_ethnicity_specific_all_std")
write_chunk(std_age, "obesity_prevalence_ethnicity_specific_age_std")
write_chunk(std_eth5, "obesity_prevalence_ethnicity_specific_eth5_std")
write_chunk(std_eth18, "obesity_prevalence_ethnicity_specific_eth18_std")
write_chunk(std_imd, "obesity_prevalence_ethnicity_specific_imd5_std")
write_chunk(std_imd10, "obesity_prevalence_ethnicity_specific_imd10_std")
write_chunk(std_reg, "obesity_prevalence_ethnicity_specific_reg_std")
write_chunk(std_lacl, "obesity_prevalence_ethnicity_specific_lacl_std")
write_chunk(std_gdp, "obesity_prevalence_ethnicity_specific_gdp_std")

write_chunk(std_lacl_30kgm2, "obesity_prevalence_30kgm2_lacl_std")

write_chunk(std_lacl_30kgm2_overall, "obesity_prevalence_30kgm2_lacl_overall_std")


write_csv(std_gdp_30kgm2, here::here("figs", "prevalence_tables", "obesity_prevalence_30kgm2_gdp_sex_std.csv"))
write_chunk(std_gdp_30kgm2_overall, "obesity_prevalence_30kgm2_gdp_overall_std")


# Prevalence grid ---------------------------------------------------------

# Write to CSV to get around the databricks/R integer encoding issue in the future
comb1_30kgm2_std |>
  write_csv(here::here("data", "ccu096_01_20191101_prevalence_comb1_30kgm2_std.csv"))

comb1_18grp_30kgm2_std |>
  write_csv(here::here("data", "ccu096_01_20191101_prevalence_comb1_18grp_30kgm2_std.csv"))

comb1_std |>
  write_csv(here::here("data", "ccu096_01_20191101_prevalence_comb1_std.csv"))

prepped_data2019 <- comb1_30kgm2_std |> 
  filter(month_date == ymd("2019-11-01")) |>
  rename(age = age_grp_5yr) |> 
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
  left_join(
    european_standard_popn_2013 |>
      mutate(
        age = str_remove(AgeGroup, " years"),
        age = str_replace(age, "plus", "+"),
        sex = str_extract(Sex, "^[M|F]")
      ),
    by = join_by(age, sex)
  ) |>
  mutate(
    numerator = if_else(is.na(numerator), 0, numerator),
    prevalence = if_else(is.na(prevalence), 0, prevalence),
    age_5yr = age,
    age = case_when(
      str_detect(age, "^1") ~ "18-19",
      str_detect(age, "^2") ~ "20-29",
      str_detect(age, "^3") ~ "30-39",
      str_detect(age, "^4") ~ "40-49",
      str_detect(age, "^5") ~ "50-59",
      str_detect(age, "^6") ~ "60-69",
      str_detect(age, "^7") ~ "70-79",
      str_detect(age, "^8") ~ "80-89",
      str_detect(age, "^9") ~ "90-99",
    )
  ) 

prepped_data <- comb1_30kgm2_std |> 
  filter(month_date == ymd("2025-04-01")) |>
  rename(age = age_grp_5yr) |> 
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
  left_join(
    european_standard_popn_2013 |>
      mutate(
        age = str_remove(AgeGroup, " years"),
        age = str_replace(age, "plus", "+"),
        sex = str_extract(Sex, "^[M|F]")
      ),
    by = join_by(age, sex)
  ) |>
  mutate(
    numerator = if_else(is.na(numerator), 0, numerator),
    prevalence = if_else(is.na(prevalence), 0, prevalence),
    age_5yr = age,
    age = case_when(
      str_detect(age, "^1") ~ "18-19",
      str_detect(age, "^2") ~ "20-29",
      str_detect(age, "^3") ~ "30-39",
      str_detect(age, "^4") ~ "40-49",
      str_detect(age, "^5") ~ "50-59",
      str_detect(age, "^6") ~ "60-69",
      str_detect(age, "^7") ~ "70-79",
      str_detect(age, "^8") ~ "80-89",
      str_detect(age, "^9") ~ "90-99",
    )
  ) 

prepped_data_18grp <- comb1_18grp_30kgm2_std |> 
  filter(month_date == ymd("2025-04-01")) |>
  rename(age = age_grp_5yr) |> 
  mutate(
    ethnicity = str_remove(combined_1lvl_18grp, "_[1-5]$"),
    imd = str_extract(combined_1lvl_18grp, "[1-5]$")
  ) |>
  left_join(
    european_standard_popn_2013 |>
      mutate(
        age = str_remove(AgeGroup, " years"),
        age = str_replace(age, "plus", "+"),
        sex = str_extract(Sex, "^[M|F]")
      ),
    by = c("age", "sex")
  ) |>
  mutate(
    numerator = if_else(is.na(numerator), 0, as.double(numerator)),
    prevalence = if_else(is.na(prevalence), 0, prevalence),
    age_5yr = age,
    age = case_when(
      str_detect(age, "^1") ~ "18-19",
      str_detect(age, "^2") ~ "20-29",
      str_detect(age, "^3") ~ "30-39",
      str_detect(age, "^4") ~ "40-49",
      str_detect(age, "^5") ~ "50-59",
      str_detect(age, "^6") ~ "60-69",
      str_detect(age, "^7") ~ "70-79",
      str_detect(age, "^8") ~ "80-89",
      str_detect(age, "^9") ~ "90-99",
    )
  ) 

prepped_data_ETHNICITY_SPECIFIC <- comb1_std |> 
  filter(month_date == ymd("2025-04-01")) |>
  rename(age = age_grp_5yr) |> 
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
  left_join(
    european_standard_popn_2013 |>
      mutate(
        age = str_remove(AgeGroup, " years"),
        age = str_replace(age, "plus", "+"),
        sex = str_extract(Sex, "^[M|F]")
      ),
    by = join_by(age, sex)
  ) |>
  mutate(
    numerator = if_else(is.na(numerator), 0, numerator),
    prevalence = if_else(is.na(prevalence), 0, prevalence),
    age_5yr = age,
    age = case_when(
      str_detect(age, "^1") ~ "18-19",
      str_detect(age, "^2") ~ "20-29",
      str_detect(age, "^3") ~ "30-39",
      str_detect(age, "^4") ~ "40-49",
      str_detect(age, "^5") ~ "50-59",
      str_detect(age, "^6") ~ "60-69",
      str_detect(age, "^7") ~ "70-79",
      str_detect(age, "^8") ~ "80-89",
      str_detect(age, "^9") ~ "90-99",
    )
  ) 

grid_crude <- prepare_prevalence_grid_crude(prepped_data) |>
  write_csv(
    here::here("figs", "prevalence_grid", "prevalence_grid_crude_14aug2025.csv")
  )

grid_std2019 <- prepare_prevalence_grid_standardised(prepped_data2019) |>
  mutate(
    across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))
  ) |>
  write_csv(
    here::here("figs", "prevalence_grid", "prevalence_grid_standardised_2019_14aug2025.csv")
  )


grid_std <- prepare_prevalence_grid_standardised(prepped_data) |>
  mutate(
    across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))
  ) |>
  write_csv(
    here::here("figs", "prevalence_grid", "prevalence_grid_standardised_14aug2025.csv")
  )

grid_std_ETHNICITY_SPECIFIC <- prepare_prevalence_grid_standardised(prepped_data_ETHNICITY_SPECIFIC) |>
  mutate(
    across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))
  ) |>
  write_csv(
    here::here("figs", "prevalence_grid", "prevalence_grid_standardised_ethnicity_specific_14aug2025.csv")
  )

grid_grp18 <- prepare_prevalence_grid_standardised(prepped_data_18grp) |>
  mutate(
    across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))
  ) |>
  write_csv(
    here::here("figs", "prevalence_grid", "prevalence_grid_standardised_18_ethnicity_groups.csv")
  )


# Sensitivity analysis by visit frequency ---------------------------------

prepare_grid_data <- function(data) {
  data |> 
    filter(month_date == ymd("2025-04-01")) |>
    rename(age = age_grp_5yr) |> 
    separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
    left_join(
      european_standard_popn_2013 |>
        mutate(
          age = str_remove(AgeGroup, " years"),
          age = str_replace(age, "plus", "+"),
          sex = str_extract(Sex, "^[M|F]")
        ),
      by = join_by(age, sex)
    ) |>
    mutate(
      numerator = if_else(is.na(numerator), 0, numerator),
      prevalence = if_else(is.na(prevalence), 0, prevalence),
      age_5yr = age,
      age = case_when(
        str_detect(age, "^1") ~ "18-19",
        str_detect(age, "^2") ~ "20-29",
        str_detect(age, "^3") ~ "30-39",
        str_detect(age, "^4") ~ "40-49",
        str_detect(age, "^5") ~ "50-59",
        str_detect(age, "^6") ~ "60-69",
        str_detect(age, "^7") ~ "70-79",
        str_detect(age, "^8") ~ "80-89",
        str_detect(age, "^9") ~ "90-99",
      )
    ) 
  
}

comb1_30kgm2_182_std %>% 
  prepare_grid_data() %>% 
  prepare_prevalence_grid_standardised() %>% 
  mutate(attendance_frequency = "0.5 years") %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  write_csv(here::here("figs", "visit_frequency_prevalence_sensitivity_analysis", "01_prevalence_grid_std_05years.csv"))
comb1_30kgm2_365_std %>% 
  prepare_grid_data() %>% 
  prepare_prevalence_grid_standardised() %>%
  mutate(attendance_frequency = "1 year") %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  write_csv(here::here("figs", "visit_frequency_prevalence_sensitivity_analysis", "02_prevalence_grid_std_1year.csv"))
comb1_30kgm2_730_std %>% 
  prepare_grid_data() %>% 
  prepare_prevalence_grid_standardised() %>%
  mutate(attendance_frequency = "2 years") %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  write_csv(here::here("figs", "visit_frequency_prevalence_sensitivity_analysis", "03_prevalence_grid_std_2years.csv"))
comb1_30kgm2_1095_std %>% 
  prepare_grid_data() %>% 
  prepare_prevalence_grid_standardised() %>%
  mutate(attendance_frequency = "3 years") %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  write_csv(here::here("figs", "visit_frequency_prevalence_sensitivity_analysis", "04_prevalence_grid_std_3years.csv"))
comb1_30kgm2_1460_std %>% 
  prepare_grid_data() %>%
  prepare_prevalence_grid_standardised() %>%
  mutate(attendance_frequency = "4 years") %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  write_csv(here::here("figs", "visit_frequency_prevalence_sensitivity_analysis", "05_prevalence_grid_std_4years.csv"))
comb1_30kgm2_1825_std %>% 
  prepare_grid_data() %>% 
  prepare_prevalence_grid_standardised() %>%
  mutate(attendance_frequency = "5 years") %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  write_csv(here::here("figs", "visit_frequency_prevalence_sensitivity_analysis", "06_prevalence_grid_std_5years.csv"))
comb1_30kgm2_std %>% 
  prepare_grid_data() %>%
  prepare_prevalance_grid_standardised() %>%
  mutate(across(c(numerator, denominator), \(x) if_else(x < 10, "REDACTED", as.character(x)))) |>
  mutate(attendence_frequency = "All years")




# Regional data -----------------------------------------------------------

# socio_demographics <- prev23 |>
#   mutate(imd_decile = as.numeric(str_extract(imd_decile, "^[0-9]{1,2}"))) |>
#   summarise(
#     median_imd = median(imd_decile),
#     median_age = median(study_start_age),
#     non_white = n_distinct(person_id[ethnicity_5_recoded != "White"]),
#     total_individuals = n_distinct(person_id),
#     percentage_non_white = non_white / total_individuals * 100,
#     .by = msoa
#   ) |>
#   mutate(across(matches("^non_white$|total"), \(.x) rnd(.x / 5, 0) * 5))

# MSOA (ethnicity-specific)
prev_25 <- prepare_msoa_prevalence_chunks(
  std_msoa_overall, 
  # socio_demographics, 
  "25",
  "msoa_obesity_prevalence"
)

# MSOA (30kg/m2 obesity threshold) overall
write_csv(std_msoa_30kgm2_overall, here::here("figs", "prevalence_tables_msoa", "msoa_obesity_prevalence_30kgm2_apr2025_all_rows.csv"))
# prev_30kgm2_25 <- prepare_msoa_prevalence_chunks(
#   std_msoa_30kgm2_overall, 
#   # socio_demographics,
#   "25", 
#   "msoa_obesity_prevalence_30kgm2"
# )

# November 2019
write_csv(std_msoa_30kgm2_2019_overall, here::here("figs", "prevalence_tables_msoa", "msoa_obesity_prevalence_30kgm2_nov2025_all_rows.csv"))

write_csv(std_msoa_30kgm2_f, here::here("figs", "prevalence_tables_msoa", "msoa_obesity_prevalence_30kgm2_females_apr2025_all_rows.csv"))
# prev_30kgm2_f_25 <- prepare_msoa_prevalence_chunks(
#   std_msoa_30kgm2_f, 
#   # socio_demographics,
#   "25", 
#   "msoa_obesity_prevalence_30kgm2_females"
# )

write_csv(std_msoa_30kgm2_m, here::here("figs", "prevalence_tables_msoa", "msoa_obesity_prevalence_30kgm2_males_apr2025_all_rows.csv"))
# prev_30kgm2_m_25 <- prepare_msoa_prevalence_chunks(
#   std_msoa_30kgm2_m, 
#   # socio_demographics,
#   "25", 
#   "msoa_obesity_prevalence_30kgm2_males"
# )

