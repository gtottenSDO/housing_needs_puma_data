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

ipums_id <- "00101"
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
  acs_ddi <- read_ipums_ddi(file_loc)
  acs_00_23 <- acs_ddi |>
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
    acs_00_23 <- define_extract_micro(
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
      read_ipums_micro()

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
xwalk <- excel_sheets("data/puma_xwalks.xlsx") |>
  set_names() |>
  map(read_excel, path = "data/puma_xwalks.xlsx")

xwalk_puma <- xwalk$xwalk_20_region |>
  left_join(
    xwalk$xwalk10_20 |>
      select(GEOID20, PUMA10, GEOID10)
  ) |>
  left_join(
    xwalk$xwalk_00_10 |>
      select(GEOID10, PUMA00)
  ) |>
  select(Region, PUMA20, PUMA10, PUMA00) |>
  pivot_longer(
    cols = starts_with("PUMA"),
    names_to = "PUMA_year",
    values_to = "PUMA"
  ) |>
  distinct() |>
  left_join(
    tibble(YEAR = c(2005:2032)) |>
      mutate(
        PUMA_year = case_when(
          YEAR >= 2022 ~ "PUMA20",
          YEAR >= 2012 ~ "PUMA10",
          YEAR >= 2005 ~ "PUMA00"
        )
      )
  ) |>
  select(YEAR, PUMA, REGION = Region) |>
  mutate(
    PUMA = as.character(as.integer(PUMA))
  )


# create survey design object
acs_hh_data <- acs_00_23$HOUSEHOLD |>
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
    across(c(GQ:VACANCY), as_factor),
    across(c(HHWT, matches("REPWT[0-9]+")), as.integer)
  ) |>
  mutate(
    across(where(is.factor), as.character),
    OWNERSHP = if_else(OWNERSHP == "N/A", "Vacant", OWNERSHP),
    VACANCY = if_else(VACANCY == "N/A", "Occupied", VACANCY)
  ) |>
  left_join(
    xwalk_puma
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
vacancy_23 <- process_vacancy(2023)

year_vec <- c(2005:2023)

vacancy_df <- year_vec |>
  map(process_vacancy) |>
  set_names(year_vec) |>
  list_rbind(names_to = "YEAR")


write_excel_csv(vacancy_df, "acs_vacancy_region.csv")

# calculate occupied
acs_occupancy <- acs_hh_srvy |>
  group_by(OWNERSHP, YEAR) |>
  survey_count(OWNERSHP, name = "occ") |>
  mutate(
    occ_moe = occ_se * 1.645,
    tenure = case_when(
      OWNERSHP == 1 ~ "own",
      OWNERSHP == 2 ~ "rent"
    )
  ) |>
  filter(!is.na(tenure)) |>
  ungroup()

# calculate total occupied and total vacant
acs_combined_ov <- acs_vacancy |>
  select(YEAR, tenure, vac, vac_moe) |>
  left_join(
    acs_occupancy |>
      select(YEAR, tenure, occ, occ_moe),
    by = c("YEAR", "tenure")
  ) %>% # old style pipe needs to be used to pipe . into bind_rows
  bind_rows(
    summarize(
      .,
      tenure = "total",
      across(c(vac:occ_moe), sum),
      .by = YEAR
    )
  ) |>
  mutate(
    total = vac + occ,
    vac_rate = vac / total
  ) |>
  group_by(tenure) |>
  mutate(avg_vac_rate = mean(vac_rate)) |>
  ungroup() |>
  mutate(
    diff_from_avg = avg_vac_rate - vac_rate,
    shortfall = diff_from_avg * total
  )
