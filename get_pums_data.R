library(tidyverse)
library(ipumsr)
library(duckplyr)
library(srvyr)
library(readxl)
library(flexiblas)
library(googledrive)

acs_samples <- get_sample_info("usa") |>
  filter(str_detect(name, pattern = "20\\d+a$")) |>
  pull(name)

drive_id <- as_id("10PZfhFsAZ-dl_w7bmohXdWUbJcxpRPlG")

ipums_id <- "00102"
ipums_dir <- "data/ipums_raw/"
ipums_ddi <- paste0("usa_", ipums_id, ".xml")
ipums_data <- paste0("usa_", ipums_id, ".dat.gz")

local_files <- list.files(
  ipums_dir,
  paste0(ipums_id, ".(dat.gz|xml)$"),
  full.names = TRUE
)


if (
  any(str_detect(local_files, ipums_ddi)) &
    any(str_detect(local_files, ipums_data))
) {
  acs_ddi <- read_ipums_ddi(paste0(ipums_dir, ipums_ddi))
  acs_data <- acs_ddi |>
    read_ipums_micro_list()
} else {
  drive_files <- drive_ls(path = drive_id)
  if (
    any(str_detect(drive_files$name, ipums_ddi)) &
      any(str_detect(drive_files$name, ipums_data))
  ) {
    drive_files %>%
      filter(str_detect(name, paste0(ipums_id, ".(dat.gz|xml)$"))) %>%
      select(name, id) |>
      pmap(
        \(name, id) {
          drive_download(id, path = paste0(ipums_dir, name), overwrite = TRUE)
        }
      )
  } else {
    acs_data <- define_extract_micro(
      collection = "usa",
      description = "ACS 1 year samples in Colorado of vacancy variables",
      samples = acs_samples,
      variables = list(
        var_spec("STATEFIP", case_selections = "08"),
        "COUNTYFIP",
        "PUMA",
        "CPUMA0010",
        "GQ",
        "OWNERSHP",
        "VACANCY",
        "REPWT",
        "AGE",
        "REPWTP"
      ),
      data_structure = "hierarchical",
    ) |>
      submit_extract() |>
      wait_for_extract() |>
      download_extract(download_dir = ipums_dir) |>
      read_ipums_micro_list()

    acs_ddi <- read_ipums_ddi(paste0(ipums_dir, ipums_ddi))

    list.files(
      ipums_dir,
      paste0(ipums_id, ".(dat.gz|xml)$"),
      full.names = TRUE
    ) |>
      map(~ drive_put(.x, path = drive_id))
  }
}

# create acs_info from the DDI which includes variable information
acs_info <- ipums_var_info(acs_ddi)

# generate table with variable information
# table will include columns {VAR}_val, and {VAR}_lbl with the value and
# label for each variable. This is setup so tab_spanner_delim can be used to
# separate out columns in a GT table
acs_var_tbl <- acs_info |>
  # only include necessary variables
  filter(var_name %in% c('OWNERSHP', 'VACANCY', "GQ")) |>
  # remove other columns
  select(var_name, val_labels) |>
  # unnest labels from nested columns
  unnest(val_labels) |>
  # group by var_name for generating row numbers to unpivot without lists
  group_by(var_name) |>
  # add row_numbers
  mutate(row = row_number()) |>
  ungroup() |>
  # pivot wider to generate table
  pivot_wider(
    names_from = var_name,
    values_from = c(val, lbl),
    # specify order of variable then val/lbl column
    names_glue = "{var_name}_{.value}",
    # use slowest so variables are grouped together
    names_vary = "slowest"
  ) |>
  # remove row index as no longer necessary
  select(-row) |>
  # move GQ to end
  relocate(starts_with("GQ"), .after = last_col())


# load in pums_region_xwalk
xwalk_puma <- read_csv("data/xwalk_puma_region.csv", col_types = "ccc")


# create survey design object
acs_hh_data <- acs_data$HOUSEHOLD |>
  filter(YEAR >= 2005, GQ %in% c("0", "1", "2", "5")) |>
  select(
    YEAR,
    PUMA,
    GQ,
    OWNERSHP,
    OWNERSHPD,
    VACANCY,
    HHWT,
    matches("REPWT[0-9]+")
  ) |>
  mutate(
    YEAR = as.integer(YEAR),
    PUMA = as.character(PUMA),
    PUMA_YEAR = case_when(
      YEAR >= 2022 ~ "PUMA20",
      YEAR >= 2012 ~ "PUMA10",
      YEAR >= 2005 ~ "PUMA00"
    ),
    across(c(GQ:VACANCY), as_factor),
    across(c(HHWT, matches("REPWT[0-9]+")), as.integer)
  ) |>
  mutate(
    across(where(is.factor), as.character),
    OWNERSHP = if_else(OWNERSHP == "N/A", "Vacant", OWNERSHP),
    VACANCY = if_else(VACANCY == "N/A", "Occupied", VACANCY)
  ) |>
  left_join(
    xwalk_puma,
    by = c("PUMA", "PUMA_YEAR"),
    relationship = "many-to-one"
  )

# create function for generating survey design object for a given year and processing it
process_vacancy <- function(year) {
  acs_hh_srvy <- acs_hh_data |>
    filter(YEAR == year) |>
    as_survey_design(
      weight = HHWT,
      repweights = matches("REPWT[0-9]+"),
      type = "ACS",
      mse = TRUE
    )

  acs_vacancy <- acs_hh_srvy |>
    group_by(YEAR, REGION, OWNERSHP, VACANCY) |>
    survey_count(VACANCY, name = "vac") |>
    mutate(
      vac_moe = vac_se * 1.645
    ) |>
    ungroup() |>
    arrange(REGION, VACANCY, OWNERSHP)

  return(acs_vacancy)
}
# vacancy_23 <- process_vacancy(2023)

year_vec <- c(2005:2024)

vacancy_df <- year_vec |>
  map(process_vacancy) |>
  set_names(year_vec) |>
  list_rbind(names_to = "YEAR")

xwalk_puma |>
  filter(REGION == "4", YEAR) |>
  distinct(PUMA_YEAR, PUMA) |>
  print(n = Inf)


write_excel_csv(vacancy_df, "acs_vacancy_region.csv")
