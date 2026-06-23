#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    12-Mar-2025
# Author:  Robert Fletcher
# Purpose: Prepare bariatric surgery data
#
#*******************************************************************************


# Load libraries and source functions -------------------------------------

source(here::here("src", "initialise_workspace.R"))


# Load data ---------------------------------------------------------------

# GLP-1RA data
bariatric_surgery <- load_data_locally(.table = "out_bariatric_surgery") 


# Clean -------------------------------------------------------------------

# Remove individuals ages <18 years or >=80 years
bariatric_surgery_1 <- bariatric_surgery |>
  rename_with(tolower) |> 
  # Code age at surgery
  mutate(surgery_age = as.numeric(date - date_of_birth) / 365.25) |> 
  filter(surgery_age >= 18 & surgery_age < 81)

# Remove individuals not in GDPPR if their surgery was after 2019
bariatric_surgery_2 <- bariatric_surgery_1 |> 
  filter(
    date >= ymd("2019-01-01"), 
    # in_gdppr == 1 | date < ymd(start_date),
    # Only keep those in the 2019 data if they died before Nov 2019
    # !(date < ymd(start_date) & is.na(date_of_death))
  )

# Remove individuals outside of England
bariatric_surgery_3 <- bariatric_surgery_2 |>
  filter(!str_detect(region, "Wales|Scotland"))

# Remove revision procedures and gastric balloons
bariatric_surgery_4 <- bariatric_surgery_3 |> 
  filter(procedure_type == "primary")

# Take the very first operation for each individual
bariatric_surgery_5 <- bariatric_surgery_4 |> 
  group_by(person_id) |> 
  arrange(date) |>
  slice(1L) |>
  ungroup()


# Check/plot --------------------------------------------------------------

plot_trends <- function(data, strata, upper_limit) {
  data |>
    mutate(date = floor_date(date, "quarter")) |>
    summarise(n = n_distinct(person_id), .by = c(date, {{ strata }})) |>
    ggplot(aes(x = date, y = n, group = {{ strata }})) +
    geom_line(aes(colour = {{ strata }})) +
    scale_y_continuous(limits = c(0, upper_limit))
}

# Sex
plot_trends(bariatric_surgery_5, sex, 1750)

# Ethnicity
plot_trends(bariatric_surgery_5, ethnicity_18_group, 2000)

# Deprivation
plot_trends(bariatric_surgery_5, imd_decile, 300)
  

# Save data ---------------------------------------------------------------

save(bariatric_surgery_5, "bariatric_surgery")
