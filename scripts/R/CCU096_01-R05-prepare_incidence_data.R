#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Prepare incidence data
#
#*******************************************************************************


# Libraries and functions -------------------------------------------------

source(here::here("src", "initialise_workspace.R"))
library(MASS)


# Variables ---------------------------------------------------------------

# New figure sub-directory
dir.create(here::here("figs", "incidence_tables"), showWarnings = FALSE)


# Load data ---------------------------------------------------------------

# Define all files in the data folder
dt <- list.files(here::here("data"))
# Define incidence tables
incidence_tables <- gsub("(^ccu096_01_\\d+_|\\.rds)", "", dt[grep("inc", dt)])
# Remove IPD data
incidence_tables <- incidence_tables[!grepl("ipd", incidence_tables)]
# Remove MSOA tables as we want to load them separately as they are in a 
# different format
incidence_no_msoa <- incidence_tables[!grepl("msoa", incidence_tables)]

for (i in 1:length(incidence_no_msoa)) {
  assign(
    gsub("incidence_", "", incidence_no_msoa[i]),
    load_data_locally(.table = incidence_no_msoa[i]) |>
    mutate(
      across(
        !matches("month|_py|inciden|age|sex|eth|imd|region|LA|comb"),
        \(x) as.integer(as.character(x))
      ),
      date = as.Date(paste0(years, "-", months, "-01"), format = "%Y-%b-%d")
    ) |>
    arrange(date) |>
    filter(date >= ymd(start_date)) |>
    mutate(row = row_number()) |>
    write_csv(here::here("output", paste0(incidence_no_msoa[i], ".csv"))),
    envir = .GlobalEnv
  )
}

# Load MSOA tables separately
# msoa <- load_data_locally(.table = "incidence_msoa") |>
#   mutate(across(matches("year|^n$"), \(x) as.integer(as.character(x))))
# msoa_30kgm2 <- load_data_locally(.table = "incidence_msoa_30kgm2") |>
#   mutate(across(matches("year|^n$"), \(x) as.integer(as.character(x))))


E# Prepare standardised incidence ------------------------------------------

prepare_inc_std <- function(data, 
                            std = european_standard_popn_2013) {
  
  data |>
    left_join(
      std |>
        mutate(
          age = str_remove(AgeGroup, " years"),
          age = str_replace(age, "plus", "+"),
          sex = str_extract(Sex, "^[M|F]")
        ),
      by = join_by(age, sex)
    ) |>
    arrange(date, age, sex) |> 
    filter(date <= ymd("2025-04-01")) |>  
    mutate(
      date_month = floor_date(date, "month"),
      # Edit as needed
      # date_month = if_else(
      #   date_month == ymd("2020-03-01"), ymd("2020-04-01"), date_month
      # ),
      date_quarter = floor_date(date_month, "quarter"),
      incidence = if_else(is.na(incidence), 0, incidence)
    )
  
}

all_std_2   <- prepare_inc_std(all_std)
eth5_std_2  <- prepare_inc_std(eth5_std)
eth18_std_2 <- prepare_inc_std(eth18_std)
imd_std_2   <- prepare_inc_std(imd_std)
imd10_std_2 <- prepare_inc_std(imd10_std)
lacl_std_2  <- prepare_inc_std(lacl_std)
reg_std_2   <- prepare_inc_std(reg_std)
comb1_std_2 <- prepare_inc_std(comb1_std)

all_30kgm2_std_2   <- prepare_inc_std(all_30kgm2_std)
eth5_30kgm2_std_2  <- prepare_inc_std(eth5_30kgm2_std)
eth18_30kgm2_std_2 <- prepare_inc_std(eth18_30kgm2_std)
imd_30kgm2_std_2   <- prepare_inc_std(imd_30kgm2_std)
imd10_30kgm2_std_2 <- prepare_inc_std(imd10_30kgm2_std)
lacl_30kgm2_std_2  <- prepare_inc_std(lacl_30kgm2_std)
reg_30kgm2_std_2   <- prepare_inc_std(reg_30kgm2_std)
comb1_30kgm2_std_2 <- prepare_inc_std(comb1_30kgm2_std)

comb1_preg_std_2 <- prepare_inc_std(comb1_preg_std)


# Calculate standardised trends overall -----------------------------------

calculate_incidence_std <- function(data, by_cols) {
  
  # By month
  std_month <- data |> 
    mutate(weighted_incidence = incidence * EuropeanStandardPopulation) |>
    summarise(
      incidence_std = sum(weighted_incidence) / sum(EuropeanStandardPopulation),
      incidence_crude = sum(n) / sum(py) * 100000,
      events = rnd(sum(n) / 5, 0) * 5, 
      py = sum(py),
      .by = all_of(by_cols)
    )
  
  # By quarter
  by_quarter <- c(by_cols, c("age", "EuropeanStandardPopulation"))
  std_quarter <- data |>
    dplyr::select(-date_month) |>
    rename(date_month = date_quarter) |>
    summarise(
      incidence_quarter = sum(n) / sum(py) * 100000,
      events = sum(n), 
      py = sum(py),
      .by = all_of(by_quarter)
    ) |>
    mutate(
      incidence_quarter = if_else(
        is.na(incidence_quarter), 0, incidence_quarter
      ),
      weighted_incidence_quarter = incidence_quarter * EuropeanStandardPopulation
    ) |>
    summarise(
      incidence_std_quarter = 
        sum(weighted_incidence_quarter) / sum(EuropeanStandardPopulation),
      incidence_crude_quarter = sum(events) / sum(py) * 100000,
      events_quarter = rnd(sum(events) / 5, 0) * 5,  
      py_quarter = sum(py),
      .by = all_of(by_cols)
    ) |>
    mutate(
      incidence_crude_quarter = if_else(
        is.na(incidence_crude_quarter), 0, incidence_crude_quarter
      )
    )
  
  out <-
    full_join(std_month, std_quarter, by = by_cols) |>
    arrange(!!!syms(by_cols))
  
  return(out)
}

# Overall
by_overall = "date_month"

std_all_overall   <- calculate_incidence_std(all_std_2, by_overall)
std_age_overall   <- calculate_incidence_std(all_std_2, c(by_overall, "age"))
std_eth5_overall  <- calculate_incidence_std(eth5_std_2, c(by_overall, "ethnicity_5_recoded"))
std_eth18_overall <- calculate_incidence_std(eth18_std_2, c(by_overall, "ethnicity_18_recoded"))
std_imd5_overall  <- calculate_incidence_std(imd_std_2, c(by_overall, "imd_quintile"))
std_imd10_overall <- calculate_incidence_std(imd10_std_2, c(by_overall, "imd_decile"))
std_reg_overall   <- calculate_incidence_std(reg_std_2, c(by_overall, "region_recoded"))
std_lacl_overall  <- calculate_incidence_std(lacl_std_2, c(by_overall, "LACLASS"))

std_all_30kgm2_overall   <- calculate_incidence_std(all_30kgm2_std_2, by_overall)
std_age_30kgm2_overall   <- calculate_incidence_std(all_30kgm2_std_2, c(by_overall, "age"))
std_eth5_30kgm2_overall  <- calculate_incidence_std(eth5_30kgm2_std_2, c(by_overall, "ethnicity_5_recoded"))
std_eth18_30kgm2_overall <- calculate_incidence_std(eth18_30kgm2_std_2, c(by_overall, "ethnicity_18_recoded"))
std_imd5_30kgm2_overall  <- calculate_incidence_std(imd_30kgm2_std_2, c(by_overall, "imd_quintile"))
std_imd10_30kgm2_overall <- calculate_incidence_std(imd10_30kgm2_std_2, c(by_overall, "imd_decile"))
std_reg_30kgm2_overall   <- calculate_incidence_std(reg_30kgm2_std_2, c(by_overall, "region_recoded"))
std_lacl_30kgm2_overall  <- calculate_incidence_std(lacl_30kgm2_std_2, c(by_overall, "LACLASS"))

# Main trends (if dplyr not updated)
std_all_not_preg_overall = all_30kgm2_std_2 |>
  mutate(weighted_incidence = incidence * EuropeanStandardPopulation) |>
  group_by(sex, date_month) |>
  summarise(
    incidence_std = sum(weighted_incidence) / sum(EuropeanStandardPopulation),
    incidence_crude = sum(n) / sum(py) * 100000,
    events = rnd(sum(n) / 5, 0) * 5, 
    py = sum(py)
  ) |> 
  filter(sex == "F")

# Pregnancy trends
std_all_preg_overall <- comb1_preg_std_2 |>
  group_by(months, years, date_month, age, sex, AgeGroup, EuropeanStandardPopulation) |>
  summarise(n = sum(n), py = sum(py), incidence = n / py * 100000) |>
  mutate(weighted_incidence = incidence * EuropeanStandardPopulation) |>
  group_by(date_month) |>
  summarise(
    incidence_std = sum(weighted_incidence) / sum(EuropeanStandardPopulation),
    incidence_crude = sum(n) / sum(py) * 100000,
    events = rnd(sum(n) / 5, 0) * 5, 
    py = sum(py)
  )

# Combine for the sensitivity analysis plot
pregnancy_sensitivity_analysis_trends <- std_all_preg_overall |> 
  dplyr::select(
    date_month, incidence_preg = incidence_std,
    events_preg = events, py_preg = py
  ) |>
  left_join(std_all_not_preg_overall)


# By sex
by = c("date_month", "sex")

std_all   <- calculate_incidence_std(all_std_2, by)
std_age   <- calculate_incidence_std(all_std_2, c(by, "age"))
std_eth5  <- calculate_incidence_std(eth5_std_2, c(by, "ethnicity_5_recoded"))
std_eth18 <- calculate_incidence_std(eth18_std_2, c(by, "ethnicity_18_recoded"))
std_imd5  <- calculate_incidence_std(imd_std_2, c(by, "imd_quintile"))
std_imd10 <- calculate_incidence_std(imd10_std_2, c(by, "imd_decile"))
std_reg   <- calculate_incidence_std(reg_std_2, c(by, "region_recoded"))
std_lacl  <- calculate_incidence_std(lacl_std_2, c(by, "LACLASS"))

std_all_30kgm2   <- calculate_incidence_std(all_30kgm2_std_2, by)
std_age_30kgm2   <- calculate_incidence_std(all_30kgm2_std_2, c(by, "age"))
std_eth5_30kgm2  <- calculate_incidence_std(eth5_30kgm2_std_2, c(by, "ethnicity_5_recoded"))
std_eth18_30kgm2 <- calculate_incidence_std(eth18_30kgm2_std_2, c(by, "ethnicity_18_recoded"))
std_imd5_30kgm2  <- calculate_incidence_std(imd_30kgm2_std_2, c(by, "imd_quintile"))
std_imd10_30kgm2 <- calculate_incidence_std(imd10_30kgm2_std_2, c(by, "imd_decile"))
std_reg_30kgm2   <- calculate_incidence_std(reg_30kgm2_std_2, c(by, "region_recoded"))
std_lacl_30kgm2  <- calculate_incidence_std(lacl_30kgm2_std_2, c(by, "LACLASS"))


# Calculate incidence grid ------------------------------------------------

inc_grid <- comb1_30kgm2_std_2 |>
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
  summarise(
    n = sum(n), 
    py = sum(py), 
    incidence = n/py * 1000, 
    .by = c(age, sex, EuropeanStandardPopulation, ethnicity, imd)
  ) |>
  mutate(
    weighted_incidence = incidence * EuropeanStandardPopulation,
      # age = case_when(
      #   str_detect(age, "^1") ~ "18-19",
      #   str_detect(age, "^2") ~ "20-29",
      #   str_detect(age, "^3") ~ "30-39",
      #   str_detect(age, "^4") ~ "40-49",
      #   str_detect(age, "^5") ~ "50-59",
      #   str_detect(age, "^6") ~ "60-69",
      #   str_detect(age, "^7") ~ "70-79",
      #   str_detect(age, "^8") ~ "80-89",
      #   str_detect(age, "^9") ~ "90+",
      # )
  )

# Output
inc_grid |>
  dplyr::select(-n) |>
  write_csv(here::here("figs", "incidence_grid", "incidence_grid_derivation_data.csv"))
  

get_incidence_grid <- function(data, by = NULL) {
  
  if(!is.null(by)) {
    
  data |>
    summarise(
      incidence_std = sum(weighted_incidence) / sum(EuropeanStandardPopulation),
      incidence_crude = sum(n) / sum(py) * 1000,
      events = rnd(sum(n) / 5, 0) * 5, 
      py = sum(py),
      .by = all_of(by)
    )
    
  } else {
    
    data |>
      summarise(
        incidence_std = sum(weighted_incidence) / sum(EuropeanStandardPopulation),
        incidence_crude = sum(n) / sum(py) * 1000,
        events = rnd(sum(n) / 5, 0) * 5, 
        py = sum(py)
      )
    
  }
}

incidence_grid <- bind_rows(
  get_incidence_grid(inc_grid) |> mutate(sex = "all", imd = "all", ethnicity = "all"),
  
  get_incidence_grid(inc_grid, "sex")       |> mutate(imd = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, "imd")       |> mutate(sex = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, "ethnicity") |> mutate(sex = "all", imd = "all"),
  
  get_incidence_grid(inc_grid, c("sex", "imd"))       |> mutate(ethnicity = "all"),
  get_incidence_grid(inc_grid, c("sex", "ethnicity")) |> mutate(imd = "all"),
  get_incidence_grid(inc_grid, c("imd", "ethnicity")) |> mutate(sex = "all"),
  
  get_incidence_grid(inc_grid, c("sex", "imd", "ethnicity"))
)

incidence_grid |>
  write_csv(here::here("figs", "incidence_grid", "incidence_grid_without_age_14aug2025.csv"))


incidence_grid_with_age <- bind_rows(
  # None
  get_incidence_grid(inc_grid) |> mutate(sex = "all", imd = "all", ethnicity = "all", age = "all"),
  
  get_incidence_grid(inc_grid, "age")       |> mutate(sex = "all", imd = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, "sex")       |> mutate(age = "all", imd = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, "imd")       |> mutate(age = "all", sex = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, "ethnicity") |> mutate(age = "all", sex = "all", imd = "all"),
  
  get_incidence_grid(inc_grid, c("age", "sex"))       |> mutate(imd = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, c("age", "imd"))       |> mutate(sex = "all", ethnicity = "all"),
  get_incidence_grid(inc_grid, c("age", "ethnicity")) |> mutate(sex = "all", imd = "all"),
  get_incidence_grid(inc_grid, c("sex", "imd"))       |> mutate(ethnicity = "all", age = "all"),
  get_incidence_grid(inc_grid, c("sex", "ethnicity")) |> mutate(imd = "all", age = "all"),
  get_incidence_grid(inc_grid, c("imd", "ethnicity")) |> mutate(age = "all", sex = "all"),
  
  get_incidence_grid(inc_grid, c("age", "sex", "imd"))       |> mutate(ethnicity = "all"),
  get_incidence_grid(inc_grid, c("age", "sex", "ethnicity")) |> mutate(imd = "all"),
  get_incidence_grid(inc_grid, c("age", "imd", "ethnicity")) |> mutate(sex = "all"),
  get_incidence_grid(inc_grid, c("sex", "imd", "ethnicity")) |> mutate(age = "all"),
  
  get_incidence_grid(inc_grid, c("age", "sex", "imd", "ethnicity"))
)

incidence_grid_with_age |>
  mutate(
    events = case_when(
      events == 5 ~ "REDACTED",
      .default = as.character(events)
    )
  ) |> 
  write_csv(here::here("figs", "incidence_grid", "incidence_grid_with_age_14aug2025.csv"))


# Calculate incidence rate ratios (overall) -------------------------------

calculate_irr_overall <- function(data, period_col, by = NULL) {
  
  by_1 <- c("period", "sex", "age", "EuropeanStandardPopulation")
  by_2 <- c("period")
  
  if (!missing(by)) {
    
    by_1 <- c(by, by_1)
    by_2 <- c(by, by_2)
    
  }
  
  data <- data |>
    filter(
      date_month %in% c(
        ymd("2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", 
            "2020-11-01", "2020-12-01", "2021-01-01", "2021-02-01",
            "2021-11-01", "2021-12-01", "2022-01-01", "2022-02-01",
            "2022-11-01", "2022-12-01", "2023-01-01", "2023-02-01",
            "2023-11-01", "2023-12-01", "2024-01-01", "2024-02-01",
            "2024-11-01", "2024-12-01", "2025-01-01", "2025-02-01")
      )
    ) |>
    mutate(
      period_winter = case_when(
        date_month < ymd("2020-03-01") ~ 0,
        date_month >= ymd("2020-11-01") & date_month <= ymd("2021-02-01") ~ 1,
        date_month >= ymd("2021-11-01") & date_month <= ymd("2022-02-01") ~ 2,
        date_month >= ymd("2022-11-01") & date_month <= ymd("2023-02-01") ~ 3,
        date_month >= ymd("2023-11-01") & date_month <= ymd("2024-02-01") ~ 4,
        date_month >= ymd("2024-11-01") & date_month <= ymd("2025-02-01") ~ 5,
      ),
      period_winter = factor(period_winter, levels = c(0:5), labels = c(0:5)),
      period_pre_during_post = case_when(
        date_month < ymd("2020-03-01") ~ 0, 
        date_month >= ymd("2020-03-01") & date_month < ymd("2022-11-01") ~ 1,
        .default = 2
      ),
      period_pre_during_post = factor(
        period_pre_during_post, levels = c(0:2), labels = c(0:2)
      ),
      period = {{ period_col }}
    ) 
  
  data_1 <- data |>
    summarise(
      incidence_period = sum(n) / sum(py) * 100000,
      events = sum(n), 
      py = sum(py),
      .by = all_of(by_1)
    )
  
  std_period <- data_1 |>
    mutate(
      incidence_period = if_else(is.na(incidence_period), 0, incidence_period),
      weighted_incidence_period = incidence_period * EuropeanStandardPopulation
    ) |>
    summarise(
      incidence_std_period = 
        sum(weighted_incidence_period) / sum(EuropeanStandardPopulation),
      incidence_crude_period = sum(events) / sum(py) * 100000,
      events_period = rnd(sum(events) / 5, 0) * 5,  
      py_period = sum(py),
      .by = all_of(by_2)
    ) |>
    mutate(
      incidence_crude_period = if_else(
        is.na(incidence_crude_period), 0, incidence_crude_period
      )
    ) |>
    arrange(period)
  
  if (missing(by)) {
    
    irr <- data_1 %>%
      glm(
        events ~ period + age + sex, offset = log(py), family = poisson(), 
        data = .
      ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE) |>
      filter(str_detect(term, "period")) |>
      dplyr::select(
        period = term, estimate, lci = conf.low, uci = conf.high, std.error
      ) |>
      mutate(period = str_remove(period, "period"))
    
    std_period <- std_period |> 
      left_join(irr, by = "period") # |>
      # left_join(
      #   data |> 
      #     summarise(n_total = sum(n))  |> 
      #     mutate(n_total = rnd(n_total / 5, 0) * 5)
      # )
    
  } else if (by == "age") {
    
    std_period <- std_period |>
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
          str_detect(age, "^9") ~ "90+"
        )
      ) |>
      summarise(
        events_period = sum(events_period),
        py_period = sum(py_period),
        incidence_std_period = events_period / py_period * 100000,
        .by = c(age, period, sex)
      ) |>
      arrange(age, sex, period)
    
    irr_f <- std_period |>
      filter(sex == "F") %>%
      split(.$age) |>
      map(
        \(x) glm(
          events_period ~ period, offset = log(py_period), family = poisson(), 
          data = x
        ) |>
          broom::tidy(exponentiate = TRUE, conf.int = TRUE)
      )
    names_f <- names(irr_f)
    irr_f <- irr_f |>
      map2(names_f, \(x, y) mutate(x, age = y)) |>
      bind_rows() |>
      dplyr::select(
        age, period = term, estimate, lci = conf.low, uci = conf.high, 
        std.error
      ) |>
      mutate(sex = "F", period = str_remove(period, "period"))
    
    irr_m <- std_period |>
      filter(sex == "M") %>%
      split(.$age) |>
      map(
        \(x) glm(
          events_period ~ period, offset = log(py_period), family = poisson(), 
          data = x
        ) |>
          broom::tidy(exponentiate = TRUE, conf.int = TRUE)
      )
    names_m <- names(irr_m)
    irr_m <- irr_m |>
      map2(names_m, \(x, y) mutate(x, age = y)) |>
      bind_rows() |>
      dplyr::select(
        age, period = term, estimate, lci = conf.low, uci = conf.high, 
        std.error
      ) |>
      mutate(sex = "M", period = str_remove(period, "period"))
    irr <- bind_rows(irr_f, irr_m)
    
    std_period <- std_period |> 
      left_join(irr, by = c("age", "sex", "period")) |>
      left_join(
        std_period |> 
          summarise(n_total = sum(events_period), .by = c("sex", "age")) |>
          mutate(n_total = rnd(n_total / 5, 0) * 5)
      )
    
  } else if (!missing(by) && by != "age") {
    
    irr_f <- data_1 |> 
      filter(sex == "F") %>%
      split(.[[by]]) |>
      map(
        \(x) glm(
          events ~ period + age, offset = log(py), 
          family = poisson(), data = x
        ) |>
          broom::tidy(exponentiate = TRUE, conf.int = TRUE)
      )
    names_f <- names(irr_f)
    irr_f <- irr_f |>
      map2(names_f, \(x, y) mutate(x, !!sym(by) := y)) |>
      bind_rows() |>
      dplyr::select(
        all_of(by), period = term, estimate, lci = conf.low, uci = conf.high, 
        std.error
      ) |>
      mutate(sex = "F", period = str_remove(period, "period"))
    
    irr_m <- data_1 |> 
      filter(sex == "M") %>%
      split(.[[by]]) |>
      map(
        \(x) glm(
          events ~ period + age, offset = log(py), 
          family = poisson(), data = x
        ) |>
          broom::tidy(exponentiate = TRUE, conf.int = TRUE)
      )
    names_m <- names(irr_m)
    irr_m <- irr_m |>
      map2(names_m, \(x, y) mutate(x, !!sym(by) := y)) |>
      bind_rows() |>
      dplyr::select(
        all_of(by), period = term, estimate, lci = conf.low, uci = conf.high, 
        std.error
      ) |>
      mutate(sex = "M", period = str_remove(period, "period"))
    
    irr <- bind_rows(irr_f, irr_m)
    std_period <- std_period |> 
      mutate(!!sym(by) := factor(std_period[[by]])) |>
      left_join(irr, by = c(by, "sex", "period")) |>
      arrange(!!sym(by), sex, period) |>
      left_join(
        data |> 
          mutate(!!sym(by) := factor(data[[by]])) |>
          summarise(n_total = sum(n), .by = c("sex", by)) |>
          mutate(n_total = rnd(n_total / 5, 0) * 5)
      )
    
    
  }
  
  return(data_1)
}


# Calculate incidence rate ratios (by sex) --------------------------------

calculate_irr <- function(data, period_col, by = NULL) {
  
  by_1 <- c("period", "sex", "age", "EuropeanStandardPopulation")
  by_2 <- c("period", "sex")
  
  if (!missing(by)) {
    
    by_1 <- c(by, by_1)
    by_2 <- c(by, by_2)
    
  }
  
  data_1 <- data |>
    filter(
      date_month %in% c(
        ymd("2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", 
            "2020-11-01", "2020-12-01", "2021-01-01", "2021-02-01",
            "2021-11-01", "2021-12-01", "2022-01-01", "2022-02-01",
            "2022-11-01", "2022-12-01", "2023-01-01", "2023-02-01",
            "2023-11-01", "2023-12-01", "2024-01-01", "2024-02-01",
            "2024-11-01", "2024-12-01", "2025-01-01", "2025-02-01")
      )
    ) |>
    mutate(
      period_winter = case_when(
        date_month < ymd("2020-03-01") ~ 0,
        date_month >= ymd("2020-11-01") & date_month <= ymd("2021-02-01") ~ 1,
        date_month >= ymd("2021-11-01") & date_month <= ymd("2022-02-01") ~ 2,
        date_month >= ymd("2022-11-01") & date_month <= ymd("2023-02-01") ~ 3,
        date_month >= ymd("2023-11-01") & date_month <= ymd("2024-02-01") ~ 4,
        date_month >= ymd("2024-11-01") & date_month <= ymd("2025-02-01") ~ 5,
      ),
      period_winter = factor(period_winter, levels = c(0:5), labels = c(0:5)),
      period_pre_during_post = case_when(
        date_month < ymd("2020-03-01") ~ 0, 
        date_month >= ymd("2020-03-01") & date_month < ymd("2022-11-01") ~ 1,
        .default = 2
      ),
      period_pre_during_post = factor(
        period_pre_during_post, levels = c(0:2), labels = c(0:2)
      ),
      period = {{ period_col }}
    ) |>
    summarise(
      incidence_period = sum(n) / sum(py) * 100000,
      events = sum(n), 
      py = sum(py),
      .by = all_of(by_1)
    )
  
  std_period <- data_1 |>
    mutate(
      incidence_period = if_else(is.na(incidence_period), 0, incidence_period),
      weighted_incidence_period = incidence_period * EuropeanStandardPopulation
    ) |>
    summarise(
      incidence_std_period = 
        sum(weighted_incidence_period) / sum(EuropeanStandardPopulation),
      incidence_crude_period = sum(events) / sum(py) * 100000,
      events_period = rnd(sum(events) / 5, 0) * 5,  
      py_period = sum(py),
      .by = all_of(by_2)
    ) |>
    mutate(
      incidence_crude_period = if_else(
        is.na(incidence_crude_period), 0, incidence_crude_period
      )
    ) |>
    arrange(sex, period)
  
  if (missing(by)) {
    
    irr_f <- data_1 |>
      filter(sex == "F") %>%
      glm(
        events ~ period + age, offset = log(py), family = poisson(), 
        data = .
      ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE) |>
      filter(str_detect(term, "period")) |>
      dplyr::select(
        period = term, estimate, lci = conf.low, uci = conf.high, std.error
      ) |>
      mutate(sex = "F", period = str_remove(period, "period"))
    
    irr_m <- data_1 |>
      filter(sex == "M") %>%
      glm(
        events ~ period + age, offset = log(py), family = poisson(), 
        data = .
      ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE) |>
      filter(str_detect(term, "period")) |>
      dplyr::select(
        period = term, estimate, lci = conf.low, uci = conf.high, std.error
      ) |>
      mutate(sex = "M", period = str_remove(period, "period"))
    
    irr <- bind_rows(irr_f, irr_m)
    std_period <- std_period |> 
      left_join(irr, by = c("period", "sex")) |>
      left_join(
        data |> 
          summarise(n_total = sum(n), .by = "sex")  |> 
          mutate(n_total = rnd(n_total / 5, 0) * 5)
      )
  
  } else if (by == "age") {
    
     std_period <- std_period |>
       mutate(
         age = case_when(
           str_detect(age, "^1") ~ "18-19",
           str_detect(age, "^2") ~ "20-29",
           str_detect(age, "^3") ~ "30-39",
           str_detect(age, "^4") ~ "40-49",
           str_detect(age, "^5") ~ "50-59",
           str_detect(age, "^6") ~ "60-69",
           str_detect(age, "^7") ~ "70-79",
           str_detect(age, "^8") ~ "80-89"
         )
       ) |>
       summarise(
         events_period = sum(events_period),
         py_period = sum(py_period),
         incidence_std_period = events_period / py_period * 100000,
         .by = c(age, period, sex)
       ) |>
       arrange(age, sex, period)
     
     irr_f <- std_period |>
       filter(sex == "F") %>%
       split(.$age) |>
       map(
         \(x) glm(
           events_period ~ period, offset = log(py_period), family = poisson(), 
           data = x
         ) |>
           broom::tidy(exponentiate = TRUE, conf.int = TRUE)
       )
     names_f <- names(irr_f)
     irr_f <- irr_f |>
       map2(names_f, \(x, y) mutate(x, age = y)) |>
       bind_rows() |>
       dplyr::select(
         age, period = term, estimate, lci = conf.low, uci = conf.high, 
         std.error
       ) |>
       mutate(sex = "F", period = str_remove(period, "period"))
     
     irr_m <- std_period |>
       filter(sex == "M") %>%
       split(.$age) |>
       map(
         \(x) glm(
           events_period ~ period, offset = log(py_period), family = poisson(), 
           data = x
         ) |>
           broom::tidy(exponentiate = TRUE, conf.int = TRUE)
       )
     names_m <- names(irr_m)
     irr_m <- irr_m |>
       map2(names_m, \(x, y) mutate(x, age = y)) |>
       bind_rows() |>
       dplyr::select(
         age, period = term, estimate, lci = conf.low, uci = conf.high, 
         std.error
       ) |>
       mutate(sex = "M", period = str_remove(period, "period"))
     irr <- bind_rows(irr_f, irr_m)
     
     std_period <- std_period |> 
       left_join(irr, by = c("age", "sex", "period")) |>
       left_join(
         std_period |> 
           summarise(n_total = sum(events_period), .by = c("sex", "age")) |>
           mutate(n_total = rnd(n_total / 5, 0) * 5)
       )
     
  } else if (!missing(by) && by != "age") {
    
    irr_f <- data_1 |> 
      filter(sex == "F") %>%
      split(.[[by]]) |>
      map(
        \(x) glm(
          events ~ period + age, offset = log(py), 
          family = poisson(), data = x
        ) |>
          broom::tidy(exponentiate = TRUE, conf.int = TRUE)
      )
    names_f <- names(irr_f)
    irr_f <- irr_f |>
      map2(names_f, \(x, y) mutate(x, !!sym(by) := y)) |>
      bind_rows() |>
      dplyr::select(
        all_of(by), period = term, estimate, lci = conf.low, uci = conf.high, 
        std.error
      ) |>
      mutate(sex = "F", period = str_remove(period, "period"))
    
    irr_m <- data_1 |> 
      filter(sex == "M") %>%
      split(.[[by]]) |>
      map(
        \(x) glm(
          events ~ period + age, offset = log(py), 
          family = poisson(), data = x
        ) |>
          broom::tidy(exponentiate = TRUE, conf.int = TRUE)
      )
    names_m <- names(irr_m)
    irr_m <- irr_m |>
      map2(names_m, \(x, y) mutate(x, !!sym(by) := y)) |>
      bind_rows() |>
      dplyr::select(
        all_of(by), period = term, estimate, lci = conf.low, uci = conf.high, 
        std.error
      ) |>
      mutate(sex = "M", period = str_remove(period, "period"))
    
    irr <- bind_rows(irr_f, irr_m)
    std_period <- std_period |> 
      mutate(!!sym(by) := factor(std_period[[by]])) |>
      left_join(irr, by = c(by, "sex", "period")) |>
      arrange(!!sym(by), sex, period) |>
      left_join(
        data |> 
          mutate(!!sym(by) := factor(data[[by]])) |>
          summarise(n_total = sum(n), .by = c("sex", by)) |>
          mutate(n_total = rnd(n_total / 5, 0) * 5)
      )
    
    
  }
  
  return(data_1)
}

all_irr_winter     <- calculate_irr(all_std_2,   period_winter)
age_irr_winter     <- calculate_irr(all_std_2,   period_winter, "age")
eth5_irr_winter    <- calculate_irr(eth5_std_2,  period_winter, "ethnicity_5_recoded")
eth18_irr_winter   <- calculate_irr(eth18_std_2, period_winter, "ethnicity_18_recoded")
imd5_irr_winter    <- calculate_irr(imd_std_2,   period_winter, "imd_quintile")
imd10_irr_winter   <- calculate_irr(imd10_std_2, period_winter, "imd_decile")
reg_irr_winter     <- calculate_irr(reg_std_2,   period_winter, "region_recoded")
lacl_irr_winter    <- calculate_irr(lacl_std_2,  period_winter, "LACLASS")

all_irr_pre_post   <- calculate_irr(all_std_2,   period_pre_during_post)
age_irr_pre_post   <- calculate_irr(all_std_2,   period_pre_during_post, "age")
eth5_irr_pre_post  <- calculate_irr(eth5_std_2,  period_pre_during_post, "ethnicity_5_recoded")
eth18_irr_pre_post <- calculate_irr(eth18_std_2, period_pre_during_post, "ethnicity_18_recoded")
imd5_irr_pre_post  <- calculate_irr(imd_std_2,   period_pre_during_post, "imd_quintile")
imd10_irr_pre_post <- calculate_irr(imd10_std_2, period_pre_during_post, "imd_decile")
reg_irr_pre_post   <- calculate_irr(reg_std_2,   period_pre_during_post, "region_recoded")
lacl_irr_pre_post  <- calculate_irr(lacl_std_2,  period_pre_during_post, "LACLASS")


# Calculate incidence rate ratios by NBM ----------------------------------

mutual <- comb1_30kgm2_std_2 |>
  filter(
    date_month %in% c(
      ymd("2019-11-01", "2019-12-01", "2020-01-01", "2020-02-01", 
          "2020-11-01", "2020-12-01", "2021-01-01", "2021-02-01",
          "2021-11-01", "2021-12-01", "2022-01-01", "2022-02-01",
          "2022-11-01", "2022-12-01", "2023-01-01", "2023-02-01",
          "2023-11-01", "2023-12-01", "2024-01-01", "2024-02-01",
          "2024-11-01", "2024-12-01", "2025-01-01", "2025-02-01")
    )
  ) |>
  mutate(
    period_winter = case_when(
      date_month < ymd("2020-03-01") ~ 0,
      date_month >= ymd("2020-11-01") & date_month <= ymd("2021-02-01") ~ 1,
      date_month >= ymd("2021-11-01") & date_month <= ymd("2022-02-01") ~ 2,
      date_month >= ymd("2022-11-01") & date_month <= ymd("2023-02-01") ~ 3,
      date_month >= ymd("2023-11-01") & date_month <= ymd("2024-02-01") ~ 4,
      date_month >= ymd("2024-11-01") & date_month <= ymd("2025-02-01") ~ 5,
    ),
    period_winter = factor(period_winter, levels = c(0:5), labels = c(0:5)),
    period_pre_during_post = case_when(
      date_month < ymd("2020-03-01") ~ 0, 
      date_month >= ymd("2020-03-01") & date_month < ymd("2022-11-01") ~ 1,
      .default = 2
    ),
    period_pre_during_post = factor(
      period_pre_during_post, levels = c(0:2), labels = c(0:2)
    ),
    period = period_pre_during_post
  ) |>
  summarise(
    incidence_period = sum(n) / sum(py) * 100000,
    events = sum(n), 
    py = sum(py),
    .by = all_of(c("combined_1lvl", "period", "sex", "age", "EuropeanStandardPopulation"))
  ) |>
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_")


## Overall ------------

# Overall trends (not stratified by anything)
irr_overall_all <- mutual %>%
  glm.nb(
    events ~ period + age + sex + ethnicity + imd + offset(log(py)),
      data = .
  ) |>
  broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
  filter(str_detect(term, "period")) |>
  dplyr::select(
    period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "Overall", sex = "All", .before = 1)

irr_overall_sex <- mutual %>%
  split(.[["sex"]]) |>
  map(
    \(x) glm.nb(
      events ~ period + age + ethnicity + imd + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
names_irr_overall_sex <- names(irr_overall_sex)
irr_overall_sex <- irr_overall_sex |>
  map2(names_irr_overall_sex, \(x, y) mutate(x, sex = y, .before = 1)) |>
  bind_rows() |>
  dplyr::select(
    sex, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "Overall", .before = 1)

overall_trends <- bind_rows(irr_overall_all, irr_overall_sex) 


## Age ----------------

irr_age_all <- mutual %>%
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
      str_detect(age, "^9") ~ "90+"
    )
  ) %>% 
  split(.[["age"]]) |>
  map(
    \(x) glm.nb(
      events ~ period + sex + ethnicity + imd + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
names_irr_age_all <- names(irr_age_all)
irr_age_all <- irr_age_all |>
  map2(names_irr_age_all, \(x, y) mutate(x, group = y, .before = 1)) |>
  bind_rows() |>
  dplyr::select(
    group, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "Age", sex = "All", .before = 1)

irr_age_sex <- mutual %>%
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
      str_detect(age, "^9") ~ "90+"
    )
  ) %>% 
  split(interaction(.[["age"]], .[["sex"]])) |>
  map(
    \(x) glm.nb(
      events ~ period + ethnicity + imd + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
irr_age_sex_all <- names(irr_age_sex)
irr_age_sex <- irr_age_sex |>
  map2(irr_age_sex_all, \(x, y) mutate(x, group = y, .before = 1)) |>
  bind_rows() |>
  separate(group, into = c("group", "sex"), sep = "[.]") |>
  dplyr::select(
    group, sex, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "Age", .before = 1)

age_trends <- bind_rows(irr_age_all, irr_age_sex)


## IMD ----------------

irr_imd_all <- mutual %>%
  split(.[["imd"]]) |>
  map(
    \(x) glm.nb(
      events ~ period + age + sex + ethnicity + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
names_irr_imd_all <- names(irr_imd_all)
irr_imd_all <- irr_imd_all |>
  map2(names_irr_imd_all, \(x, y) mutate(x, group = y, .before = 1)) |>
  bind_rows() |>
  dplyr::select(
    group, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "IMD", sex = "All", .before = 1)

irr_imd_sex <- mutual %>%
  split(interaction(.[["imd"]], .[["sex"]])) |>
  map(
    \(x) glm.nb(
      events ~ period + age + ethnicity + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
names_imd_sex <- names(irr_imd_sex)
irr_imd_sex <- irr_imd_sex |>
  map2(names_imd_sex, \(x, y) mutate(x, group = y, .before = 1)) |>
  bind_rows() |>
  separate(group, into = c("group", "sex"), sep = "[.]") |>
  dplyr::select(
    group, sex, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "IMD", .before = 1)

imd_trends <- bind_rows(irr_imd_all, irr_imd_sex)


## Ethnicity ----------

irr_eth_all <- mutual %>%
  split(.[["ethnicity"]]) |>
  map(
    \(x) glm.nb(
      events ~ period + age + sex + imd + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
names_irr_eth_all <- names(irr_eth_all)
irr_eth_all <- irr_eth_all |>
  map2(names_irr_eth_all, \(x, y) mutate(x, group = y, .before = 1)) |>
  bind_rows() |>
  dplyr::select(
    group, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "Ethnicity", sex = "All", .before = 1)

irr_eth_sex <- mutual %>%
  split(interaction(.[["ethnicity"]], .[["sex"]])) |>
  map(
    \(x) glm.nb(
      events ~ period + age + imd + offset(log(py)),
      data = x
    ) |>
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)  |>
      filter(str_detect(term, "period")) 
  )
names_eth_sex <- names(irr_eth_sex)
irr_eth_sex <- irr_eth_sex |>
  map2(names_eth_sex, \(x, y) mutate(x, group = y, .before = 1)) |>
  bind_rows() |>
  separate(group, into = c("group", "sex"), sep = "[.]") |>
  dplyr::select(
    group, sex, period = term, estimate, lci = conf.low, uci = conf.high, 
    std.error, p = p.value
  ) |>
  mutate(strata = "Ethnicity", .before = 1)

ethnicity_trends <- bind_rows(irr_eth_all, irr_eth_sex)

bind_rows(
  overall_trends,
  age_trends,
  imd_trends,
  ethnicity_trends
) |>
  write_csv(here::here("figs", "incidence_rate_ratios", "incidence_rate_ratios_pre_post_30kgm2_14aug2025.csv"))


# Attendance frequency sensitivity analysis -------------------------------

num_182  <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_numerator_comb1_30kgm2_182_std.rds") %>% 
  pivot_longer(cols = !c(months, years, sex, age_grp_5yr), names_to = "combined_1lvl", values_to = "events") %>%
  mutate(events = if_else(is.na(events), 0, as.integer(events)), age = age_grp_5yr, years = as.character(years)) %>%
  dplyr::select(-age_grp_5yr)
num_365  <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_numerator_comb1_30kgm2_365_std.rds") %>% 
  pivot_longer(cols = !c(months, years, sex, age_grp_5yr), names_to = "combined_1lvl", values_to = "events") %>%
  mutate(events = if_else(is.na(events), 0, as.integer(events)), age = age_grp_5yr, years = as.character(years)) %>%
  dplyr::select(-age_grp_5yr)
num_730  <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_numerator_comb1_30kgm2_730_std.rds") %>% 
  pivot_longer(cols = !c(months, years, sex, age_grp_5yr), names_to = "combined_1lvl", values_to = "events") %>%
  mutate(events = if_else(is.na(events), 0, as.integer(events)), age = age_grp_5yr, years = as.character(years)) %>%
  dplyr::select(-age_grp_5yr)
num_1095 <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_numerator_comb1_30kgm2_1095_std.rds") %>% 
  pivot_longer(cols = !c(months, years, sex, age_grp_5yr), names_to = "combined_1lvl", values_to = "events") %>%
  mutate(events = if_else(is.na(events), 0, as.integer(events)), age = age_grp_5yr, years = as.character(years)) %>%
  dplyr::select(-age_grp_5yr)
num_1460 <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_numerator_comb1_30kgm2_1460_std.rds")%>% 
  pivot_longer(cols = !c(months, years, sex, age_grp_5yr), names_to = "combined_1lvl", values_to = "events") %>%
  mutate(events = if_else(is.na(events), 0, as.integer(events)), age = age_grp_5yr, years = as.character(years)) %>%
  dplyr::select(-age_grp_5yr)
num_1825 <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_numerator_comb1_30kgm2_1825_std.rds")%>% 
  pivot_longer(cols = !c(months, years, sex, age_grp_5yr), names_to = "combined_1lvl", values_to = "events") %>%
  mutate(events = if_else(is.na(events), 0, as.integer(events)), age = age_grp_5yr, years = as.character(years)) %>%
  dplyr::select(-age_grp_5yr)

den_182  <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_denominator_comb1_30kgm2_182_std.rds") %>%
  mutate(years = as.character(years))
den_365  <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_denominator_comb1_30kgm2_365_std.rds") %>%
  mutate(years = as.character(years))
den_730  <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_denominator_comb1_30kgm2_730_std.rds") %>%
  mutate(years = as.character(years))
den_1095 <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_denominator_comb1_30kgm2_1095_std.rds") %>%
  mutate(years = as.character(years))
den_1460 <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_denominator_comb1_30kgm2_1460_std.rds") %>%
  mutate(years = as.character(years))
den_1825 <- readRDS("~/collab/CCU096_01/data/ccu096_01_20191101_tmp_incidence_denominator_comb1_30kgm2_1825_std.rds") %>%
  mutate(years = as.character(years))

prep <- function(denominator, numerator, freq, filename) {
  denominator %>%
    left_join(numerator) %>%
    filter(!is.na(events)) %>%
    mutate(
      date = as.Date(paste0(years, "-", months, "-01"), format = "%Y-%b-%d")
    ) |>
    arrange(date) |>
    filter(date >= ymd(start_date)) |> 
    left_join(
      european_standard_popn_2013 |>
        mutate(
          age = str_remove(AgeGroup, " years"),
          age = str_replace(age, "plus", "+"),
          sex = str_extract(Sex, "^[M|F]")
        ),
      by = join_by(age, sex)
    ) %>%
    mutate(incidence = events / py * 1000) %>%
    dplyr::select(-events) %>%
    mutate(attendance_frequency = freq) %>%
    write_csv(here::here("figs", "visit_frequency_incidence_sensitivity_analysis", filename))
}

total_182  <- prep(den_182,  num_182,  "6 months", "01_obesity_incidence_05years.csv")
total_365  <- prep(den_365,  num_365,  "1 year",   "02_obesity_incidence_1year.csv")
total_730  <- prep(den_730,  num_730,  "2 years",  "03_obesity_incidence_2years.csv")
total_1095 <- prep(den_1095, num_1095, "3 years",  "04_obesity_incidence_3years.csv") 
total_1460 <- prep(den_1460, num_1460, "4 years",  "05_obesity_incidence_4years.csv") 
total_1825 <- prep(den_1825, num_1825, "5 years",  "06_obesity_incidence_5years.csv") 




# Save --------------------------------------------------------------------

write_chunk <- function(data, file_name) {
  
  # if (nrow(data) > 1000 && nrow(data) <= 2000) {
  # 
  #   data_1 <- data |> slice(1:1000)
  #   data_2 <- data |> slice(1001:2000)
  #     
  #   write_csv(
  #     data_1, 
  #     here::here("figs", "incidence_tables", paste0(file_name, "_chunk_01.csv"))
  #   )
  #    
  #   write_csv(
  #     data_2, 
  #     here::here("figs", "incidence_tables", paste0(file_name, "_chunk_02.csv"))
  #   )
  #   
  # } else if (nrow(data) > 2000 && nrow(data) <= 3000) {
  #   
  #   data_1 <- data |> slice(1:1000)
  #   data_2 <- data |> slice(1001:2000)
  #   data_3 <- data |> slice(2001:3000)
  #   
  #   write_csv(
  #     data_1, 
  #     here::here("figs", "incidence_tables", paste0(file_name, "_chunk_01.csv"))
  #   )
  #   
  #   write_csv(
  #     data_2, 
  #     here::here("figs", "incidence_tables", paste0(file_name, "_chunk_02.csv"))
  #   )
  #   
  #   write_csv(
  #     data_3, 
  #     here::here("figs", "incidence_tables", paste0(file_name, "_chunk_03.csv"))
  #   )
  #     
  # } else if (nrow(data) <= 1000) {
  #   
    write_csv(
      data, here::here("figs", "incidence_tables", paste0(file_name, ".csv"))
    )
    
  # } else {
  #   
  #   message("You need more segments.")
  #   
  # }
  
}

# Standardised trends
write_chunk(std_all,   "inc_all_std")
write_chunk(std_age,   "inc_age_std")
write_chunk(std_eth5,  "inc_eth5_std")
write_chunk(std_eth18, "inc_eth18_std")
write_chunk(std_imd5,  "inc_imd5_std")
write_chunk(std_imd10, "inc_imd10_std")
write_chunk(std_reg,   "inc_reg_std")
write_chunk(std_lacl,  "inc_lacl_std")


# Standardised (30kg/m2) OVERALL
write_chunk(std_all_30kgm2_overall,   "inc_all_30kgm2_std_overall_apr2025")
write_chunk(std_age_30kgm2_overall,   "inc_age_30kgm2_std_overall_apr2025")
write_chunk(c,  "inc_eth5_30kgm2_std_overall_apr2025")
write_chunk(std_imd5_30kgm2_overall,  "inc_imd5_30kgm2_std_overall_apr2025")

# Standardised (30kg/m2) BY SEX
write_chunk(std_all_30kgm2,   "inc_all_30kgm2_std_apr2025")
write_chunk(std_age_30kgm2,   "inc_age_30kgm2_std_apr2025")
write_chunk(std_eth5_30kgm2,  "inc_eth5_30kgm2_std_apr2025")
write_chunk(std_imd5_30kgm2,  "inc_imd5_30kgm2_std_apr2025")


write_irr <- function(data, file) { 
  
  write_csv(
    data, here::here("figs", "incidence_rate_ratios", paste0(file, ".csv"))
  )
  
}

# Incidence rate ratios (by winter)
write_irr(all_irr_winter,   "all_irr_winter")
write_irr(age_irr_winter,   "age_irr_winter")
write_irr(eth5_irr_winter,  "eth5_irr_winter")
write_irr(eth18_irr_winter, "eth18_irr_winter")
write_irr(imd5_irr_winter,  "imd5_irr_winter")
write_irr(imd10_irr_winter, "imd10_irr_winter")
write_irr(reg_irr_winter,   "reg_irr_winter")
write_irr(lacl_irr_winter,  "lacl_irr_winter")

# Incidence rate ratios (by pre/during/post pandemic)
write_irr(all_irr_pre_post,   "all_irr_pre_post")
write_irr(age_irr_pre_post,   "age_irr_pre_post")
write_irr(eth5_irr_pre_post,  "eth5_irr_pre_post")
write_irr(eth18_irr_pre_post, "eth18_irr_pre_post")
write_irr(imd5_irr_pre_post,  "imd5_irr_pre_post")
write_irr(imd10_irr_pre_post, "imd10_irr_pre_post")
write_irr(reg_irr_pre_post,   "reg_irr_pre_post")
write_irr(lacl_irr_pre_post,  "lacl_irr_pre_post")


# Incidence by 18 ethnicity groups ----------------------------------------

eth18_numerator <- tmp_numerator_comb1_18grp_30kgm2_std |> 
  dplyr::select(-row) |>
  pivot_longer(
    cols = !c(months, years, sex, age_grp_5yr, date),
    names_to = "combined_1lvl_18grp", values_to = "n"
  ) |>
  mutate(
    age = age_grp_5yr,
    ethnicity = str_remove(combined_1lvl_18grp, "_[1-5]$"),
    imd       = str_extract(combined_1lvl_18grp, "[1-5]$"),
    n = if_else(is.na(n), 0, as.double(n))
  ) |>
  filter(date < ymd("2025-05-01")) |>
  dplyr::select(-c(combined_1lvl_18grp, age_grp_5yr))

eth18_denominator <- tmp_denominator_comb1_18grp_30kgm2_std |>
  mutate(
    ethnicity = str_remove(combined_1lvl_18grp, "_[1-5]$"),
    imd       = str_extract(combined_1lvl_18grp, "[1-5]$")
  ) |>
  filter(date < ymd("2025-05-01")) |>
  dplyr::select(-c(row, combined_1lvl_18grp))

eth18_combined <- eth18_numerator |>
  left_join(
    eth18_denominator, 
    by = c("months", "years", "sex", "date", "age", "ethnicity", "imd")
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
  mutate(py = if_else(is.na(py), 0, as.double(py)))


eth18_overall <- eth18_combined |> 
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
      str_detect(age, "^9") ~ "90+"
    ),
    incidence = n / py * 1000,
    weighted_incidence = incidence * EuropeanStandardPopulation
  ) |>
  group_by(age, ethnicity, imd) |>
  summarise(
    incidence_std = sum(weighted_incidence) / sum(EuropeanStandardPopulation),
    incidence_crude = sum(n) / sum(py) * 1000,
    events = rnd(sum(n) / 5, 0) * 5, 
    py = sum(py)
  ) |>
  ungroup() |>
  mutate(events = if_else(events < 10, "REDACTED", as.character(events)))

eth18_by_sex <- eth18_combined |> 
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
      str_detect(age, "^9") ~ "90+"
    )
  ) |>
  group_by(age, sex, ethnicity, imd) |> 
  summarise(events = sum(n), py = sum(py), incidence = sum(n) / sum(py) * 1000) |>
  ungroup() |>
  mutate(events = if_else(events < 10, "REDACTED", as.character(events)))


# Export data for pregnancy sensitivity analysis --------------------------

# Normal
main_age_inc <- comb1_30kgm2_std_2 |> 
  filter(sex == "F") |>
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
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
      str_detect(age, "^9") ~ "90+"
    )
  ) |>
  group_by(age, sex, ethnicity, imd) |>
  summarise(n = sum(n), py = sum(py), incidence = n / py * 100000) |>
  ungroup() |>
  mutate(n = rnd(n / 5, 0) * 5) |>
  mutate(analysis = "main", .before = 1)

# Pregnancy derivation data
comb1_preg_std_2 |>
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
  group_by(age, sex, imd, ethnicity) |>
  summarise(
    n = sum(n), py = sum(py)
  ) |>
  ungroup() |>
  mutate(incidence = n / py * 100000) |>
  dplyr::select(-n) |>
  write_csv(here::here("figs", "pregnancy_sensitivity_analysis", "pregnancy_sensitivity_analysis_derivation_data.csv"))

# Pregnancy 
pregnancy_age_inc <- comb1_preg_std_2 |>
  separate(combined_1lvl, into = c("ethnicity", "imd"), sep = "_") |>
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
      str_detect(age, "^9") ~ "90+"
    )
  ) |>
  group_by(age, sex, ethnicity, imd) |>
  summarise(n = sum(n), py = sum(py), incidence = n / py * 100000) |>
  ungroup() |>
  mutate(n = rnd(n / 5, 0) * 5) |>
  mutate(analysis = "pregnancy", .before = 1)

pregnancy_sensitivity_analysis_trends |>
  write_csv(here::here("figs", "pregnancy_sensitivity_analysis", "01_pregnancy_sensitivity_analysis_trends.csv"))

bind_rows(main_age_inc, pregnancy_age_inc) |>
  mutate(n = if_else(n == 5, "REDACTED", as.character(n))) |>
  write_csv(here::here("figs", "pregnancy_sensitivity_analysis", "02_pregnancy_sensitivity_analysis_by_age.csv"))
