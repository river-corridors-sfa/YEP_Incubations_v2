# =============================================================================
# 01_model_DO.R -- best-fit kL and Vmax for every microcosm, across Km values
# -----------------------------------------------------------------------------
# Model (one combined model for every vial):
#
#   dDO/dt = -kL * DO  -  Vmax * DO / (Km + DO)
#
#   kL   : first-order (DO-dependent) rate constant, h-1
#   Vmax : DO-independent rate at saturating DO, mg O2 L-1 h-1
#   Km   : half-saturation constant, mg L-1. FIXED in each fit (not estimated)
#   DO0  : DO at time zero, mg L-1. Estimated together with kL and Vmax
#
# Km is varied from 0.05 to 10 mg/L:
#   - 0.13-1.22 mg/L is the Sihi et al. (2020) posterior range (main figure)
#   - 10 mg/L is the value used in the submitted manuscript (SI figure)
#   - 0.2 mg/L is the base value used for the treatment statistics
#
# Per-mass parameters use the liters of water per kg dry sediment in each vial:
#   kL_kg   = kL   * L_water_per_kg   (L kg-1 h-1)
#   Vmax_kg = Vmax * L_water_per_kg   (mg O2 kg-1 h-1)
#
# Inputs : data_files$raw_do, data_files$water_mass
# Outputs: outputs/modeling/YEP_DataPrep_Log.csv
#          outputs/modeling/YEP_DO_Series_Prepared.csv   (data actually fitted)
#          outputs/modeling/YEP_BestFit_AllKm.csv        (all vials x all Km)
#          outputs/modeling/YEP_BestFit_DOabove2_Km0.2.csv (sensitivity refit)
#          outputs/figures/FigS2_QC_Model_Fits_Km0.2.pdf/.png
# Run time: about 5-10 min.
# =============================================================================

rm(list = ls(all = TRUE))
library(deSolve)      # ode()
library(minpack.lm)   # nls.lm()
library(dplyr)
library(ggplot2)
source("src/00_data_paths.R")

# ---- Settings ----------------------------------------------------------------
Km_values <- c(0.05, 0.1, 0.13, 0.2, 0.3, 0.51, 0.75, 1, 1.22, 3.2, 5.28, 10)  # mg/L
Km_base   <- 0.2                                                                  # mg/L

trim_start_min   <- 2     # minutes removed at the start of each run
trim_end_min     <- 2     # minutes removed at the end of each run
anoxic_threshold <- 0.1   # mg/L; series is cut before the first reading <= this
DO_min_sensitivity <- 2   # mg/L; sensitivity refit uses only readings >= this

# Parameter bounds for the fit (kL and Vmax may go to zero)
lower <- c(kL = 0,  Vmax = 0,   DO0 = 0)
upper <- c(kL = 50, Vmax = 500, DO0 = 15)

out_dir <- file.path("outputs", "modeling")
fig_dir <- file.path("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Read DO data and label each vial ------------------------------------
do_raw <- read.csv(data_files$raw_do, stringsAsFactors = FALSE)

do_data <- do_raw %>%
  filter(!grepl("^YEP_INC-W", Sample_Name)) %>%       # water-only controls: see 07_water_controls.R
  mutate(Elapsed_Seconds = as.numeric(Elapsed_Seconds),
         DO_mg_per_L     = as.numeric(DO_mg_per_L)) %>%
  filter(!is.na(Elapsed_Seconds), !is.na(DO_mg_per_L))

# Sample names look like YEP1A_INC-H2:
#   YEP1 = dry, YEP2 = wet; A/B/C = incubation day 1/2/3; S/U/H = treatment
vial_info <- do_data %>%
  distinct(Sample_Name) %>%
  mutate(Moisture  = ifelse(substr(Sample_Name, 4, 4) == "1", "Dry", "Wet"),
         Day       = recode(substr(Sample_Name, 5, 5), A = "1", B = "2", C = "3"),
         Treatment = recode(sub(".*_INC-([SUH]).*", "\\1", Sample_Name),
                            S = "Control", U = "Unburned + DOC", H = "High Burn + DOC")) %>%
  arrange(Sample_Name)

# ---- 2. Liters of water per kg dry sediment in each vial --------------------
mass <- read.csv(data_files$water_mass, stringsAsFactors = FALSE) %>%
  filter(Sample_Name %in% vial_info$Sample_Name) %>%
  mutate(L_water_per_kg = as.numeric(Water_Mass_g) / as.numeric(Dry_Sediment_Mass_g)) %>%
  select(Sample_Name, Dry_Sediment_Mass_g, Water_Mass_g, L_water_per_kg)

stopifnot(nrow(mass) == nrow(vial_info))   # every vial needs a water mass
vial_info <- left_join(vial_info, mass, by = "Sample_Name")

# ---- 3. Prepare the DO series of each vial ----------------------------------
# For each vial: trim 2 min at each end, cut at the DO minimum (removes the
# end-of-run rise), and cut before the first reading <= 0.1 mg/L (anoxic
# plateau). Time is in hours from the first retained reading.
series   <- list()
prep_log <- data.frame()

for (s in vial_info$Sample_Name) {
  d <- do_data %>% filter(Sample_Name == s) %>% arrange(Elapsed_Seconds)
  t_h <- (d$Elapsed_Seconds - min(d$Elapsed_Seconds)) / 3600
  DO  <- d$DO_mg_per_L
  n_raw <- length(DO)

  keep <- t_h >= trim_start_min / 60 & t_h <= max(t_h) - trim_end_min / 60
  t_h <- t_h[keep]; DO <- DO[keep]

  i_min <- which.min(DO)                       # cut after the DO minimum
  t_h <- t_h[1:i_min]; DO <- DO[1:i_min]

  i_anox <- which(DO <= anoxic_threshold)[1]   # cut before anoxia
  reached_anoxia <- !is.na(i_anox)
  if (reached_anoxia) { t_h <- t_h[1:(i_anox - 1)]; DO <- DO[1:(i_anox - 1)] }

  t_h <- t_h - t_h[1]
  series[[s]] <- data.frame(time = t_h, DO = DO)
  prep_log <- rbind(prep_log, data.frame(
    Sample_Name = s, n_raw = n_raw, n_fit = length(DO),
    duration_fit_min = round(60 * max(t_h), 1), reached_anoxia = reached_anoxia,
    DO_first_obs = DO[1], DO_last_obs = DO[length(DO)]))
}
write.csv(left_join(vial_info, prep_log, by = "Sample_Name"),
          file.path(out_dir, "YEP_DataPrep_Log.csv"), row.names = FALSE)
# Save the prepared series so 02_model_mcmc.R fits exactly the same data
write.csv(bind_rows(series, .id = "Sample_Name"),
          file.path(out_dir, "YEP_DO_Series_Prepared.csv"), row.names = FALSE)

# ---- 4. Fitting function ----------------------------------------------------
# One function because the same fit is repeated for 30 vials x 12 Km values
# (and again for the DO >= 2 sensitivity fit in section 6). It solves the ODE,
# minimizes the squared residuals with nls.lm from three starting guesses, and
# keeps the best one.
#   - DO0 is the starting DO of the fitted curve, at the first fitted reading
#     (2 min into the incubation, after the start trim). The model needs it to
#     draw the curve, but it is not a result we interpret.
#   - If the ODE solver or nls.lm fails, the failure is recorded (n_starts_failed,
#     fit_ok) instead of stopping the whole run.
#   - nls_info is the nls.lm exit code; 1-4 mean it converged normally.
#   - n_starts_agree: how many of the 3 starts reached the same answer (diagnostic only).
model_DO <- function(kL, Vmax, DO0, Km, times) {
  rhs <- function(t, y, p) list(-kL * max(y, 0) - Vmax * max(y, 0) / (Km + max(y, 0)))
  out <- tryCatch(ode(y = c(DO = DO0), times = times, func = rhs, parms = NULL,
                      method = "lsoda", rtol = 1e-8, atol = 1e-10),
                  error = function(e) NULL)
  if (is.null(out) || nrow(out) != length(times)) return(rep(NA, length(times)))
  out[, "DO"]
}

fit_vial <- function(dat, Km) {
  residuals_fn <- function(p) {
    r <- dat$DO - model_DO(p[["kL"]], p[["Vmax"]], p[["DO0"]], Km, dat$time)
    r[!is.finite(r)] <- 1e3                      # solver failure is penalized, not hidden
    r
  }
  run_nls <- function(par) {                     # returns NULL if nls.lm stops with an error
    tryCatch(nls.lm(par, fn = residuals_fn, lower = lower, upper = upper,
                    control = nls.lm.control(maxiter = 300, ftol = 1e-12, ptol = 1e-12)),
             error = function(e) NULL)
  }
  slope <- abs(coef(lm(DO ~ time, data = dat))[[2]])   # mg/L/h, for starting guesses
  starts <- list(c(kL = 0.01,                       Vmax = slope,     DO0 = dat$DO[1]),
                 c(kL = slope / mean(dat$DO),       Vmax = 0.1,       DO0 = dat$DO[1]),
                 c(kL = slope / mean(dat$DO) / 2,   Vmax = slope / 2, DO0 = dat$DO[1]))
  fits <- list()
  for (i in seq_along(starts)) {
    s  <- pmin(pmax(starts[[i]], lower + 1e-6), upper - 1e-6)
    f1 <- run_nls(s)
    if (is.null(f1)) next                        # this start failed; try the next one
    # Polish step: run nls.lm again from its own answer. This removes a small
    # stopping error when a parameter sits on its bound (e.g. kL = 0 in dry vials).
    f2 <- run_nls(f1$par)
    fits[[length(fits) + 1]] <- if (!is.null(f2) && sum(f2$fvec^2) <= sum(f1$fvec^2)) f2 else f1
  }
  
  # Every start failed: return a row of NA so the run continues and the vial is flagged
  if (length(fits) == 0) {
    return(data.frame(kL_per_h = NA, Vmax_mg_L_h = NA, DO0_est = NA, n_points = nrow(dat),
                      RMSE = NA, R2 = NA, Resid_Lag1_Autocorr = NA, n_starts_agree = 0,
                      n_starts_failed = length(starts), nls_info = NA, fit_ok = FALSE,
                      kL_at_lower = NA, Vmax_at_lower = NA))
  }
  
  ssr  <- sapply(fits, function(f) sum(f$fvec^2))
  best_fit <- fits[[which.min(ssr)]]
  best <- best_fit$par
  agree <- sapply(fits, function(f)
    abs(sum(f$fvec^2) - min(ssr)) <= 0.01 * min(ssr) &&
      all(abs(f$par[1:2] - best[1:2]) <= pmax(0.02 * abs(best[1:2]), 0.05)))
  
  res <- dat$DO - model_DO(best[["kL"]], best[["Vmax"]], best[["DO0"]], Km, dat$time)
  data.frame(kL_per_h = best[["kL"]], Vmax_mg_L_h = best[["Vmax"]], DO0_est = best[["DO0"]],
             n_points = length(res), RMSE = sqrt(mean(res^2)),
             R2 = 1 - sum(res^2) / sum((dat$DO - mean(dat$DO))^2),
             Resid_Lag1_Autocorr = cor(res[-1], res[-length(res)]),
             n_starts_agree = sum(agree),
             n_starts_failed = length(starts) - length(fits),
             nls_info = best_fit$info,
             fit_ok = best_fit$info %in% 1:4,
             kL_at_lower = best[["kL"]] < 1e-3, Vmax_at_lower = best[["Vmax"]] < 1e-2)
}
# ---- 5. Fit every vial at every Km ------------------------------------------
best_fits <- data.frame()
for (Km in Km_values) {
  message("Km = ", Km)
  for (s in vial_info$Sample_Name) {
    best_fits <- rbind(best_fits, data.frame(Sample_Name = s, Km_mg_L = Km, fit_vial(series[[s]], Km)))
  }
}

# Per-mass parameters (see header)
best_fits <- best_fits %>%
  left_join(vial_info, by = "Sample_Name") %>%
  mutate(kL_L_kg_h    = kL_per_h    * L_water_per_kg,
         Vmax_mg_kg_h = Vmax_mg_L_h * L_water_per_kg) %>%
  relocate(Moisture, Treatment, Day, .after = Sample_Name)
write.csv(best_fits, file.path(out_dir, "YEP_BestFit_AllKm.csv"), row.names = FALSE)

# ---- 6. Sensitivity refit: only readings with DO >= 2 mg/L, at Km = 0.2 -----
# Answers R1 (Methods 5): some vials went close to anoxia, others did not.
sens <- data.frame()
for (s in vial_info$Sample_Name) {
  dat <- series[[s]] %>% filter(DO >= DO_min_sensitivity)
  sens <- rbind(sens, data.frame(Sample_Name = s, Km_mg_L = Km_base, fit_vial(dat, Km_base)))
}
sens <- sens %>%
  left_join(vial_info, by = "Sample_Name") %>%
  mutate(kL_L_kg_h = kL_per_h * L_water_per_kg, Vmax_mg_kg_h = Vmax_mg_L_h * L_water_per_kg)
write.csv(sens, file.path(out_dir, "YEP_BestFit_DOabove2_Km0.2.csv"), row.names = FALSE)

# ---- 7. QC figure: observed vs fitted DO for every vial at Km = 0.2 ---------
qc <- data.frame()
for (s in vial_info$Sample_Name) {
  p <- best_fits %>% filter(Sample_Name == s, Km_mg_L == Km_base)
  qc <- rbind(qc, data.frame(Sample_Name = s, time_min = 60 * series[[s]]$time,
                             Observed = series[[s]]$DO,
                             Fitted = model_DO(p$kL_per_h, p$Vmax_mg_L_h, p$DO0_est, Km_base, series[[s]]$time)))
}
qc$Panel <- sub("_INC-", " ", qc$Sample_Name)   # shorter panel labels, e.g. "YEP1A S1"
qc_plot <- ggplot(qc, aes(time_min)) +
  geom_point(aes(y = Observed), size = 0.3, color = "gray55") +
  geom_line(aes(y = Fitted), color = "black", linewidth = 0.5) +
  facet_wrap(~ Panel, scales = "free_x", ncol = 6) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +   # fewer ticks so labels do not overlap
  labs(x = "Time (min)", y = expression(DO~(mg~L^{-1})),
       caption = "Gray: observed DO (5-s readings); black: model fit with Km = 0.2 mg/L.") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 9), axis.text = element_text(size = 7.5),
        panel.spacing = unit(0.15, "lines"))
# Saved at the size it is printed in the SI (landscape page), so the text stays readable
ggsave(file.path(fig_dir, "FigS2_QC_Model_Fits_Km0.2.pdf"), qc_plot, width = 9, height = 6)
ggsave(file.path(fig_dir, "FigS2_QC_Model_Fits_Km0.2.png"), qc_plot, width = 9, height = 6, dpi = 300)

# ---- 8. Console check -------------------------------------------------------
print(best_fits %>% filter(Km_mg_L == Km_base) %>% group_by(Moisture) %>%
        summarize(kL_L_kg_h = mean(kL_L_kg_h), Vmax_mg_kg_h = mean(Vmax_mg_kg_h),
                  starts_agree = min(n_starts_agree), RMSE_max = max(RMSE)))
