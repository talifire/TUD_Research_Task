# ============================================================
# helpers.R  — Slot Deviation Analysis helpers
# Author: Tetiana Alifirenko
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(here)
  library(scales)
})

# ------------------------------------------------------
# 1) Project utilities
# ------------------------------------------------------

init_project_folders <- function() {
  dir.create(here("plots"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(here("tables"), showWarnings = FALSE, recursive = TRUE)
  invisible(TRUE)
}

save_plot <- function(plot, filename, w = 9, h = 6, dpi = 300) {
  ggsave(
    filename = here("plots", filename),
    plot     = plot,
    width    = w,
    height   = h,
    dpi      = dpi
  )
}

save_table_csv <- function(df, filename) {
  readr::write_csv(df, here("tables", filename))
}

write_session_info <- function(filename = "sessionInfo.txt") {
  writeLines(capture.output(sessionInfo()), here("tables", filename))
}

# ------------------------------------------------------
# 2) Data validation & prep
# ------------------------------------------------------

assert_required_columns <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

prep_slot_core <- function(slot_core, tz = "UTC") {
  # Make types consistent (safe for publishing / reproducibility)
  slot_core |>
    mutate(
      DEV_MIN = as.numeric(DEV_MIN),
      PHASE   = as.character(PHASE),
      AIRPORT = as.character(AIRPORT),
      CDM     = as.character(CDM),
      SLOT_TIME = as.POSIXct(SLOT_TIME, tz = tz)  # adjust if your data uses local time
    ) |>
    mutate(
      PHASE = factor(PHASE, levels = c("ARR", "DEP"))
    )
}

# ------------------------------------------------------
# 3) Weibull fitting helpers
# ------------------------------------------------------

trim_upper <- function(x, upper_trim = 0.995) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(x)
  thr <- quantile(x, upper_trim, na.rm = TRUE)
  x[x <= thr]
}
fit_weibull_shape_scale <- function(x_pos, upper_trim = 0.995, min_n = 30) {
  x_pos <- x_pos[is.finite(x_pos) & x_pos > 0]
  n_raw <- length(x_pos)
  
  if (n_raw < min_n) {
    return(list(shape = NA_real_, scale = NA_real_, n_raw = n_raw, converged = FALSE))
  }
  
  thr    <- quantile(x_pos, upper_trim, na.rm = TRUE)
  x_trim <- x_pos[x_pos <= thr]
  if (length(x_trim) < min_n) {
    return(list(shape = NA_real_, scale = NA_real_, n_raw = n_raw, converged = FALSE))
  }
  
  weibull_nll <- function(par, x) {
    if (any(!is.finite(par))) return(1e12)
    
    shape <- exp(par[1])
    scale <- exp(par[2])
    
    if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) return(1e12)
    
    ll <- suppressWarnings(dweibull(x, shape = shape, scale = scale, log = TRUE))
    
    if (any(!is.finite(ll))) return(1e12)
    
    -sum(ll)
  }
  
  init_scale <- max(stats::median(x_trim), 1e-6)
  init <- c(log(1), log(init_scale))
  
  fit <- optim(
    par = init,
    fn  = weibull_nll,
    x   = x_trim,
    method  = "BFGS",
    control = list(maxit = 2000)
  )
  
  converged <- is.list(fit) && is.finite(fit$value) && fit$convergence == 0
  if (!converged) {
    return(list(shape = NA_real_, scale = NA_real_, n_raw = n_raw, converged = FALSE))
  }
  
  list(
    shape = exp(fit$par[1]),
    scale = exp(fit$par[2]),
    n_raw = n_raw,
    converged = TRUE
  )
}


weibull_q <- function(p, shape, scale) {
  if (is.na(shape) || is.na(scale)) return(NA_real_)
  qweibull(p, shape = shape, scale = scale)
}

fit_tail_q <- function(x_mag, p = 0.90, upper_trim = 0.995, min_n = 30) {
  # Returns a list with fit + requested quantile for positive magnitudes
  fit <- fit_weibull_shape_scale(x_mag, upper_trim = upper_trim, min_n = min_n)
  q   <- if (fit$converged) weibull_q(p, fit$shape, fit$scale) else NA_real_
  c(fit, list(q = q))
}

# ------------------------------------------------------
# 4) Buffer / adherence helpers
# ------------------------------------------------------

compute_empirical_bounds <- function(df, lower_p = 0.10, upper_p = 0.90) {
  df |>
    group_by(AIRPORT, PHASE) |>
    summarise(
      lower = quantile(DEV_MIN, lower_p, na.rm = TRUE),
      upper = quantile(DEV_MIN, upper_p, na.rm = TRUE),
      .groups = "drop"
    )
}

add_bounds_and_flag <- function(df, bounds_df,
                                lower_col = "lower", upper_col = "upper",
                                flag_name = "ADHERENT") {
  df |>
    left_join(bounds_df, by = c("AIRPORT", "PHASE")) |>
    mutate(
      "{flag_name}" := DEV_MIN >= .data[[lower_col]] & DEV_MIN <= .data[[upper_col]]
    )
}

compute_adherence_for_coverage <- function(slot_core, cov) {
  lower_p <- (1 - cov) / 2
  upper_p <- 1 - lower_p
  
  bounds <- compute_empirical_bounds(slot_core, lower_p = lower_p, upper_p = upper_p)
  
  slot_with_buf <- slot_core |>
    left_join(bounds, by = c("AIRPORT", "PHASE")) |>
    mutate(
      ADHERENT     = DEV_MIN >= lower & DEV_MIN <= upper,
      buffer_width = upper - lower
    )
  
  slot_with_buf |>
    group_by(CDM, AIRPORT, PHASE) |>
    summarise(
      coverage_target  = cov,
      n_flights        = n(),
      adherence_rate   = mean(ADHERENT, na.rm = TRUE),
      mean_buf_width   = mean(buffer_width, na.rm = TRUE),
      median_buf_width = median(buffer_width, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_hybrid_bounds_two_tailed <- function(slot_core,
                                             emp_lower_p = 0.10,
                                             emp_upper_p = 0.90,
                                             weib_p = 0.90,
                                             upper_trim = 0.995,
                                             min_n = 30) {
  emp_bounds <- slot_core |>
    group_by(AIRPORT, PHASE) |>
    summarise(
      q10_emp = quantile(DEV_MIN, emp_lower_p, na.rm = TRUE),
      q90_emp = quantile(DEV_MIN, emp_upper_p, na.rm = TRUE),
      .groups = "drop"
    )
  
  weib_tails <- slot_core |>
    group_by(AIRPORT, PHASE) |>
    group_modify(~{
      dev <- .x$DEV_MIN
      early_mag <- -dev[dev < 0]
      late_mag  <-  dev[dev > 0]
      
      early_fit <- fit_tail_q(early_mag, p = weib_p, upper_trim = upper_trim, min_n = min_n)
      late_fit  <- fit_tail_q(late_mag,  p = weib_p, upper_trim = upper_trim, min_n = min_n)
      
      tibble(
        n_early     = early_fit$n_raw,
        shape_early = early_fit$shape,
        scale_early = early_fit$scale,
        q_weib_early_mag = early_fit$q,
        lower_weib  = ifelse(is.na(early_fit$q), NA_real_, -early_fit$q),
        
        n_late      = late_fit$n_raw,
        shape_late  = late_fit$shape,
        scale_late  = late_fit$scale,
        q_weib_late = late_fit$q,
        upper_weib  = late_fit$q
      )
    }) |>
    ungroup()
  
  emp_bounds |>
    left_join(weib_tails, by = c("AIRPORT", "PHASE")) |>
    mutate(
      lower_hybrid_both = pmin(q10_emp, lower_weib, na.rm = TRUE),
      upper_hybrid_both = pmax(q90_emp, upper_weib, na.rm = TRUE)
    )
}

# ---------------------------
# 5) Plot helpers
# ---------------------------

plot_deviations_off_on <- function(slot_core) {
  slot_core_plot <- slot_core |>
    mutate(
      PHASE_LABEL = dplyr::case_when(
        PHASE == "DEP" ~ "Departure: off-block vs slot",
        PHASE == "ARR" ~ "Arrival: on-block vs slot",
        TRUE ~ as.character(PHASE)
      )
    )
  
  ggplot(slot_core_plot, aes(x = DEV_MIN, fill = AIRPORT)) +
    geom_histogram(
      position = "identity",
      alpha = 0.45,
      bins = 80,
      colour = "grey30"
    ) +
    geom_vline(
      xintercept = 0,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    facet_grid(PHASE_LABEL ~ AIRPORT, scales = "free_y") +
    coord_cartesian(xlim = c(-60, 60)) +
    labs(
      title = "Deviations between scheduled and actual off-/on-block times",
      subtitle = "Negative = earlier than slot, Positive = later than slot",
      x = "Deviation (BLOCK − SLOT) [min]",
      y = "Flight count",
      fill = "Airport"
    ) +
    theme_minimal(base_size = 14)
}

summarise_deviations <- function(slot_core) {
  slot_core |>
    group_by(AIRPORT, PHASE) |>
    summarise(
      n          = n(),
      mean_dev   = mean(DEV_MIN, na.rm = TRUE),
      median_dev = median(DEV_MIN, na.rm = TRUE),
      sd_dev     = sd(DEV_MIN, na.rm = TRUE),
      p10        = quantile(DEV_MIN, 0.10, na.rm = TRUE),
      p90        = quantile(DEV_MIN, 0.90, na.rm = TRUE),
      .groups    = "drop"
    )
}

plot_weibull_global <- function(slot_core, upper_trim = 0.995, q = 0.90,
                                x_max = 200, bins = 40) {
  dev_pos <- slot_core$DEV_MIN
  dev_pos <- dev_pos[is.finite(dev_pos) & dev_pos > 0]
  dev_pos <- trim_upper(dev_pos, upper_trim = upper_trim)
  
  fit <- fit_weibull_shape_scale(dev_pos, upper_trim = 1.0) # already trimmed
  shape_hat <- fit$shape
  scale_hat <- fit$scale
  buffer_q  <- weibull_q(q, shape_hat, scale_hat)
  
  x_grid <- seq(0, x_max, length.out = 300)
  dens   <- dweibull(x_grid, shape = shape_hat, scale = scale_hat)
  
  ggplot(data.frame(del = dev_pos), aes(x = del)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = bins,
      fill = "grey80",
      color = "grey40"
    ) +
    geom_line(
      data = data.frame(x = x_grid, y = dens),
      aes(x = x, y = y),
      linewidth = 1.1
    ) +
    geom_vline(xintercept = buffer_q, linetype = "dashed", linewidth = 1) +
    coord_cartesian(xlim = c(0, x_max)) +
    labs(
      title = "Weibull-Based Operational Deviation Model",
      subtitle = paste0(
        "Shape = ", round(shape_hat, 2),
        ", Scale = ", round(scale_hat, 2),
        ", Buffer (q", round(q * 100), ") ≈ ", round(buffer_q, 1), " min"
      ),
      x = "Deviation (minutes)",
      y = "Density"
    ) +
    theme_minimal(base_size = 16)
}

plot_adherence_vs_coverage <- function(buffer_sensitivity) {
  ggplot(
    buffer_sensitivity,
    aes(
      x = coverage_target,
      y = adherence_rate,
      color = CDM,
      linetype = PHASE,
      group = interaction(CDM, AIRPORT, PHASE)
    )
  ) +
    geom_line(alpha = 0.6, linewidth = 1) +
    geom_point(alpha = 0.7, size = 2) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title    = "Adherence rate as a function of buffer coverage",
      x        = "Central buffer coverage (quantile band)",
      y        = "Observed adherence rate",
      linetype = "Phase"
    ) +
    theme_minimal(base_size = 14)
}


plot_width_vs_adherence <- function(buffer_sensitivity) {
  ggplot(
    buffer_sensitivity,
    aes(
      x = mean_buf_width,
      y = adherence_rate,
      color = CDM,
      linetype = PHASE
    )
  ) +
    geom_line(aes(group = interaction(CDM, AIRPORT, PHASE)), alpha = 0.5) +
    geom_point(alpha = 0.7) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title    = "Trade-off between buffer width and adherence",
      x        = "Mean buffer width [min] (upper - lower quantile)",
      y        = "Adherence rate",
      linetype = "Phase"
    ) +
    theme_minimal(base_size = 14)
}

compute_delta_table <- function(buffer_sensitivity, cov_lo = 0.70, cov_hi = 0.80) {
  buffer_sensitivity |>
    filter(coverage_target %in% c(cov_lo, cov_hi)) |>
    select(CDM, AIRPORT, PHASE, coverage_target, adherence_rate, mean_buf_width) |>
    pivot_wider(
      names_from  = coverage_target,
      values_from = c(adherence_rate, mean_buf_width),
      names_glue  = "{.value}_cov{coverage_target*100}"
    ) |>
    mutate(
      delta_adherence = .data[[paste0("adherence_rate_cov", cov_hi * 100)]] -
        .data[[paste0("adherence_rate_cov", cov_lo * 100)]],
      delta_width = .data[[paste0("mean_buf_width_cov", cov_hi * 100)]] -
        .data[[paste0("mean_buf_width_cov", cov_lo * 100)]]
    )
}

fit_weibull_by_tail <- function(slot_core) {
  weib_late <- slot_core |>
    group_by(AIRPORT, PHASE) |>
    group_modify(~{
      late_mag <- .x$DEV_MIN[.x$DEV_MIN > 0]
      fit <- fit_weibull_shape_scale(late_mag)
      tibble(n_raw = fit$n_raw, shape = fit$shape, scale = fit$scale)
    }) |>
    ungroup()
  
  weib_early <- slot_core |>
    group_by(AIRPORT, PHASE) |>
    group_modify(~{
      early_mag <- - .x$DEV_MIN[.x$DEV_MIN < 0]
      fit <- fit_weibull_shape_scale(early_mag)
      tibble(n_raw = fit$n_raw, shape = fit$shape, scale = fit$scale)
    }) |>
    ungroup()
  
  list(weib_late = weib_late, weib_early = weib_early)
}

plot_weibull_positive_compare <- function(slot_core, weib_late, x_max = 100) {
  x_grid <- seq(0, 200, length.out = 400)
  
  dens_late_df <- tidyr::expand_grid(
    AIRPORT = unique(slot_core$AIRPORT),
    PHASE   = unique(slot_core$PHASE),
    x       = x_grid
  ) |>
    left_join(weib_late, by = c("AIRPORT", "PHASE")) |>
    mutate(
      density = ifelse(
        is.na(shape) | is.na(scale),
        NA_real_,
        dweibull(x, shape = shape, scale = scale)
      )
    )
  
  ggplot(
    dens_late_df,
    aes(
      x = x, y = density,
      color = interaction(AIRPORT, PHASE),
      linetype = PHASE
    )
  ) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    labs(
      title    = "Weibull fits for positive slot deviations",
      subtitle = "Late operations: comparison by airport and phase",
      x        = "Positive deviation (minutes)",
      y        = "Density",
      color    = "Airport–Phase",
      linetype = "Phase"
    ) +
    coord_cartesian(xlim = c(0, x_max)) +
    theme_minimal(base_size = 14)
}

plot_weibull_negative_compare <- function(slot_core, weib_early, x_min = -60) {
  x_grid <- seq(0, 200, length.out = 400)
  
  dens_early_df <- tidyr::expand_grid(
    AIRPORT = unique(slot_core$AIRPORT),
    PHASE   = unique(slot_core$PHASE),
    x_mag   = x_grid
  ) |>
    left_join(weib_early, by = c("AIRPORT", "PHASE")) |>
    mutate(
      x       = -x_mag,
      density = ifelse(
        is.na(shape) | is.na(scale),
        NA_real_,
        dweibull(x_mag, shape = shape, scale = scale)
      )
    )
  
  ggplot(
    dens_early_df,
    aes(
      x = x, y = density,
      color = interaction(AIRPORT, PHASE),
      linetype = PHASE
    )
  ) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    labs(
      title    = "Weibull fits for negative slot deviations",
      subtitle = "Early operations (mirrored): comparison by airport and phase",
      x        = "Negative deviation (minutes)",
      y        = "Density",
      color    = "Airport–Phase",
      linetype = "Phase"
    ) +
    coord_cartesian(xlim = c(x_min, 0)) +
    theme_minimal(base_size = 14)
}

plot_hybrid_bounds <- function(hybrid_bounds_both) {
  ggplot(
    hybrid_bounds_both,
    aes(
      x = PHASE,
      ymin = lower_hybrid_both,
      ymax = upper_hybrid_both,
      color = AIRPORT
    )
  ) +
    geom_errorbar(width = 0.25, linewidth = 1.1) +
    labs(
      title = "Hybrid Buffer Bounds (Empirical + Weibull)",
      y     = "Deviation range [min]",
      x     = "Phase",
      color = "Airport"
    ) +
    theme_minimal(base_size = 14)
}

compute_slot_with_hybrid <- function(slot_core, hybrid_bounds_both) {
  slot_core |>
    left_join(
      hybrid_bounds_both |>
        select(AIRPORT, PHASE, lower_hybrid_both, upper_hybrid_both),
      by = c("AIRPORT", "PHASE")
    ) |>
    mutate(
      ADHERENT_HYB_BOTH = if_else(
        DEV_MIN >= lower_hybrid_both & DEV_MIN <= upper_hybrid_both,
        1L, 0L
      )
    )
}

summarise_hybrid_adherence <- function(slot_with_hybrid_both) {
  slot_with_hybrid_both |>
    group_by(CDM, PHASE) |>
    summarise(
      n_flights   = n(),
      adherent_n  = sum(ADHERENT_HYB_BOTH, na.rm = TRUE),
      adherence_h = adherent_n / n_flights,
      .groups = "drop"
    )
}

plot_hybrid_adherence_cdm_phase <- function(adherence_cdm_phase) {
  ggplot(adherence_cdm_phase, aes(x = PHASE, y = adherence_h, fill = CDM)) +
    geom_col(position = position_dodge(width = 0.6)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                               limits = c(0, 1)) +
    labs(
      title    = "Hybrid buffer adherence by CDM environment and phase",
      subtitle = "Two-tailed hybrid buffer (empirical + Weibull, early and late deviations)",
      x        = "Phase (ARR / DEP)",
      y        = "Adherence rate",
      fill     = "Environment"
    ) +
    theme_minimal(base_size = 14)
}

plot_queueing_load <- function(slot_with_hybrid_both, fixed_buffer = 15) {
  sim <- slot_with_hybrid_both |>
    arrange(SLOT_TIME) |>
    mutate(
      delay_pos       = pmax(DEV_MIN, 0),
      workload_fixed  = cumsum(pmax(delay_pos - fixed_buffer, 0)),
      workload_hybrid = cumsum(pmax(delay_pos - upper_hybrid_both, 0))
    )
  
  ggplot(sim, aes(x = SLOT_TIME)) +
    geom_line(aes(y = workload_fixed,  color = paste0("Fixed ±", fixed_buffer), group = 1),
              linewidth = 1.0, na.rm = TRUE) +
    geom_line(aes(y = workload_hybrid, color = "Hybrid buffer", group = 1),
              linewidth = 1.0, na.rm = TRUE) +
    facet_wrap(~ CDM, scales = "free_y") +
    labs(
      title = "Queueing Load Before vs After Hybrid Buffer",
      x     = "Time",
      y     = "Accumulated delay (proxy workload)",
      color = "Method"
    ) +
    theme_minimal(base_size = 14)
}

