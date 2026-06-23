#*******************************************************************************
#
# Project: CCU096_01
# Title:   Whole-population trends in obesity and its treatment
# Date:    17-Feb-2025
# Author:  Robert Fletcher
# Purpose: Install libraries
#
#*******************************************************************************


# Install libraries -------------------------------------------------------

libs <- c("here", "tidyverse", "emmeans")
repo <- "https://packages.sde.digital.nhs.uk/repository/cran-mirror/"

for (l in libs) install.packages(l, repos = repo, dependencies = TRUE)
