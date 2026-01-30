my_apdf <- read_csv(here("data", "my_apdf_2025_JanJun_RAW.csv"))


set.seed(42)
airports_of_interest <- c("EDDB","EDDF","EDDC","EDDN")
slot_core <- my_apdf |>
  rename(
    ADEP       = ADEP_ICAO,
    ADES       = ADES_ICAO,
    SLOT_TIME  = SCHED_TIME_UTC,
    BLOCK_TIME = BLOCK_TIME_UTC,
    PHASE      = SRC_PHASE
  ) |>
  mutate(
    SLOT_TIME  = ymd_hms(SLOT_TIME, tz = "UTC"),
    BLOCK_TIME = ymd_hms(BLOCK_TIME, tz = "UTC")
  ) |>
  # --- restrict to one-week window ---------------------------------
filter(
  SLOT_TIME >= ymd_hms("2025-02-03 00:00:00", tz = "UTC"),
  SLOT_TIME <= ymd_hms("2025-02-07 23:59:59", tz = "UTC")
) |>
  # --- keep only airports of interest -------------------------------
filter(
  (PHASE == "DEP" & ADEP %in% airports_of_interest) |
    (PHASE == "ARR" & ADES %in% airports_of_interest)
) |>
  # --- anonymise other airports ------------------------------------
mutate(
  ADEP = if_else(
    ADEP %in% airports_of_interest,
    ADEP,
    sprintf("%04d", sample(0:9999, n(), replace = TRUE))
  ),
  ADES = if_else(
    ADES %in% airports_of_interest,
    ADES,
    sprintf("%04d", sample(0:9999, n(), replace = TRUE))
  ),
  AIRPORT = if_else(PHASE == "DEP", ADEP, ADES),
  DEV_MIN = as.numeric(difftime(BLOCK_TIME, SLOT_TIME, units = "mins")),
  ABS_DEV = abs(DEV_MIN),
  CDM     = if_else(AIRPORT %in% c("EDDF", "EDDB"), "CDM", "non-CDM")
) |>
  filter(!is.na(DEV_MIN)) |>
  # --- trim extreme deviations (distribution-based) -----------------
mutate(
  dev_lower = quantile(DEV_MIN, 0.005),
  dev_upper = quantile(DEV_MIN, 0.995)
) |>
  filter(DEV_MIN >= dev_lower & DEV_MIN <= dev_upper) |>
  select(-dev_lower, -dev_upper)

