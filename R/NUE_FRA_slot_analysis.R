###############################################################################
# PROJECT: Slot Deviation Analysis & Buffer-Based Adherence Indicator
# AUTHOR: Tetiana Alifirenko
# DATE: 08.12.2025
###############################################################################

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(purrr)
})

# 0) Helpers + folders
source(here::here("R", "helpers.R"))
init_project_folders()
write_session_info()

# 1) Data loading
slot_core <- readr::read_csv(here::here("data_public", "slot_core.csv"),
                             show_col_types = FALSE)

required_cols <- c("DEV_MIN","AIRPORT","PHASE","CDM","SLOT_TIME")
assert_required_columns(slot_core, required_cols)
slot_core <- prep_slot_core(slot_core, tz = "UTC")

# 2) Tables
desc_by_ap <- summarise_deviations(slot_core)
save_table_csv(desc_by_ap, "desc_by_ap.csv")

coverage_grid <- seq(0.50, 0.95, by = 0.05)
buffer_sensitivity <- purrr::map_df(coverage_grid, ~compute_adherence_for_coverage(slot_core, .x))
save_table_csv(buffer_sensitivity, "buffer_sensitivity.csv")

delta_table <- compute_delta_table(buffer_sensitivity, cov_lo = 0.70, cov_hi = 0.80)
save_table_csv(delta_table, "delta_table_cov70_80.csv")

hybrid_bounds_both <- compute_hybrid_bounds_two_tailed(slot_core)
save_table_csv(hybrid_bounds_both, "hybrid_bounds_both.csv")

slot_with_hybrid_both <- compute_slot_with_hybrid(slot_core, hybrid_bounds_both)
adherence_cdm_phase <- summarise_hybrid_adherence(slot_with_hybrid_both)
save_table_csv(adherence_cdm_phase, "adherence_cdm_phase.csv")

# Weibull fits by group (for plotting + saving)
weib_fits <- fit_weibull_by_tail(slot_core)
save_table_csv(weib_fits$weib_late,  "weibull_late_by_airport_phase.csv")
save_table_csv(weib_fits$weib_early, "weibull_early_by_airport_phase.csv")

# 3) Plots
p1 <- plot_deviations_off_on(slot_core)
save_plot(p1, "PPT_DEV_OFF_ON_BLOCK.png")

p2 <- plot_adherence_vs_coverage(buffer_sensitivity)
save_plot(p2, "PPT_Buffer_Coverage_vs_Adherence.png")

p3 <- plot_width_vs_adherence(buffer_sensitivity)
save_plot(p3, "PPT_Buffer_Width_vs_Adherence.png")

p4 <- plot_weibull_global(slot_core)
save_plot(p4, "PPT_WEIBULL_BUFFER.png")

p5 <- plot_weibull_negative_compare(slot_core, weib_fits$weib_early)
save_plot(p5, "PPT_Weibull_Negative_Compare.png")

p6 <- plot_weibull_positive_compare(slot_core, weib_fits$weib_late)
save_plot(p6, "PPT_Weibull_Positive_Compare.png")

p7 <- plot_hybrid_adherence_cdm_phase(adherence_cdm_phase)
save_plot(p7, "PPT_Hybrid_Adherence_CDM_Phase.png")

p8 <- plot_hybrid_bounds(hybrid_bounds_both)
save_plot(p8, "PPT_Hybrid_Buffer_Bounds.png")

p9 <- plot_queueing_load(slot_with_hybrid_both, fixed_buffer = 15)
save_plot(p9, "PPT_Queuing_Load_Comparison.png")

# Optional: print key outputs in console
print(desc_by_ap)
print(delta_table)
print(adherence_cdm_phase)
