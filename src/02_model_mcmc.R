# =============================================================================
# 02_model_mcmc.R -- Bayesian (MCMC) estimates of kL and Vmax for every vial
# -----------------------------------------------------------------------------
# Same model as 01_model_DO.R:  dDO/dt = -kL * DO - Vmax * DO / (Km + DO)
#
# What is estimated:  kL, Vmax and DO0 (plus the residual variance sigma2).
#   - Km is FIXED at each value of the Sihi et al. (2020) posterior range.
#     Each Km gives a separate posterior, conditional on that Km (a
#     sensitivity analysis, not one posterior averaged over Km).
#   - DO0 is sampled (not fixed at its best fit), so its uncertainty carries
#     into kL and Vmax. Its prior is normal, centered on the first fitted DO
#     reading with SD = DO0_prior_sd (see mcmc_sampler.R for why).
#   - rho (AR(1) residual autocorrelation) is FIXED at the lag-1
#     autocorrelation of the best-fit residuals from 01_model_DO.R.
#     Sampling rho as well did not converge in tests; 02b_mcmc_checks.R shows
#     how much the credible intervals depend on rho.
#   The sampler itself (priors, likelihood, steps) is in src/mcmc_sampler.R.
#
# Chains: 4 per vial, started from different balances of kL and Vmax
#   (see dispersed_starts() in mcmc_sampler.R) to test for other solutions.
#
# Convergence (posterior package): rank-normalized split R-hat, bulk and tail
#   effective sample size (ESS), for kL, Vmax, DO0 and sigma2, computed from
#   all post-burn-in draws. A vial-Km fit is flagged converged when every
#   R-hat < 1.01 and every bulk and tail ESS > 400.
#
# Inputs : outputs/modeling/YEP_DO_Series_Prepared.csv, YEP_BestFit_AllKm.csv
# Outputs: outputs/modeling/mcmc/YEP_MCMC_Summary_Km<Km>.csv (one per Km)
#          outputs/modeling/mcmc/YEP_MCMC_Draws_Km<Km>.rds  (every 10th draw,
#                                                            incl. burn-in)
#          outputs/modeling/YEP_MCMC_Summary.csv               (all Km combined)
#          outputs/figures/FigS_MCMC_*.pdf/.png
# Packages: deSolve, posterior, dplyr, tidyr, ggplot2
#   (install.packages("posterior") if it is missing)
# Run time: about 20-40 min per Km, 2-4.5 h for all 7 Km depending on the computer. To split the work,
#   give Km values on the command line, e.g.  Rscript src/02_model_mcmc.R 0.13 0.2
#   Then run  Rscript src/02_model_mcmc.R combine  to combine and plot only.
# =============================================================================

rm(list = ls(all = TRUE))
library(deSolve)
library(posterior)  # rhat(), ess_bulk(), ess_tail()
library(dplyr)
library(tidyr)
library(ggplot2)
source("src/mcmc_sampler.R")

# ---- Settings ----------------------------------------------------------------
Km_values <- c(0.13, 0.2, 0.3, 0.51, 0.75, 1, 1.22)   # Sihi et al. (2020) range, mg/L
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) Km_values <- suppressWarnings(as.numeric(args))   # optional: a subset of Km
Km_values <- Km_values[!is.na(Km_values)]   # "combine" on the command line skips sampling

n_iter  <- 10000      # iterations per chain
n_burn  <- 3000       # burn-in (proposal tuning; discarded)
thin_saved <- 10      # every 10th draw is saved to the .rds file (summaries use all draws)
set.seed(42)

# Prior bounds (same as the best-fit bounds in 01_model_DO.R)
lower <- c(kL = 0,  Vmax = 0,   DO0 = 0)
upper <- c(kL = 50, Vmax = 500, DO0 = 15)
DO0_prior_sd <- 0.2   # mg/L; SD of the normal prior on DO0 around the first fitted reading

# Convergence thresholds
rhat_max <- 1.01
ess_min  <- 400

in_dir   <- file.path("outputs", "modeling")
mcmc_dir <- file.path(in_dir, "mcmc")
fig_dir  <- file.path("outputs", "figures")
dir.create(mcmc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

# ---- 1. Inputs from 01_model_DO.R and checks ---------------------------------
series    <- read.csv(file.path(in_dir, "YEP_DO_Series_Prepared.csv"))
best_fits <- read.csv(file.path(in_dir, "YEP_BestFit_AllKm.csv"))
vials     <- sort(unique(series$Sample_Name))

# The AR(1) likelihood assumes evenly spaced readings: check the spacing (s)
spacing_s <- unlist(tapply(series$time, series$Sample_Name, function(t) diff(t) * 3600))
cat("Spacing between readings (s): min", round(min(spacing_s), 2), "max", round(max(spacing_s), 2), "\n")
stopifnot(all(spacing_s > 4 & spacing_s < 6))   # all close to 5 s

# ---- 2. Run the sampler for every Km and every vial -------------------------
for (Km in Km_values) {
  message("Km = ", Km)
  summary_rows <- data.frame()
  saved_draws  <- data.frame()

  for (s in vials) {
    dat <- series %>% filter(Sample_Name == s) %>% arrange(time)
    bf  <- best_fits %>% filter(Sample_Name == s, Km_mg_L == Km)
    stopifnot(nrow(bf) == 1, nrow(dat) >= 20, all(diff(dat$time) > 0), bf$L_water_per_kg > 0)
    rho <- max(bf$Resid_Lag1_Autocorr, 0)       # AR(1) coefficient, fixed

    # Four chains from different starting points
    starts <- dispersed_starts(bf$kL_per_h, bf$Vmax_mg_L_h, bf$DO0_est)
    chains <- list()
    for (ch in 1:4) {
      chains[[ch]] <- run_chain(dat, Km, rho, starts[[ch]], lower, upper, n_iter, n_burn, DO0_prior_sd)
    }

    # Post-burn-in draws as a posterior draws_array: iterations x chains x parameters
    post_array <- array(NA, dim = c(n_iter - n_burn, 4, 4),
                        dimnames = list(NULL, NULL, c("kL", "Vmax", "DO0", "sigma2")))
    for (ch in 1:4) post_array[, ch, ] <- chains[[ch]][(n_burn + 1):n_iter, ]
    post_array <- as_draws_array(post_array)

    # Convergence for every sampled parameter
    rhat_all <- sapply(c("kL", "Vmax", "DO0", "sigma2"), function(p) rhat(extract_variable_matrix(post_array, p)))
    bulk_all <- sapply(c("kL", "Vmax", "DO0", "sigma2"), function(p) ess_bulk(extract_variable_matrix(post_array, p)))
    tail_all <- sapply(c("kL", "Vmax", "DO0", "sigma2"), function(p) ess_tail(extract_variable_matrix(post_array, p)))

    # All post-burn-in draws pooled, and per-kg versions
    post    <- do.call(rbind, lapply(chains, function(m) m[(n_burn + 1):n_iter, ]))
    kL_kg   <- post[, "kL"]   * bf$L_water_per_kg
    Vmax_kg <- post[, "Vmax"] * bf$L_water_per_kg

    summary_rows <- rbind(summary_rows, data.frame(
      Sample_Name = s, Km_mg_L = Km,
      Moisture = bf$Moisture, Treatment = bf$Treatment, Day = bf$Day,
      L_water_per_kg = bf$L_water_per_kg, rho_AR1 = rho,
      kL_bestfit = bf$kL_per_h, Vmax_bestfit = bf$Vmax_mg_L_h, DO0_bestfit = bf$DO0_est,
      DO_first_obs = dat$DO[1],
      kL_median    = median(post[, "kL"]),
      kL_lower95   = quantile(post[, "kL"], 0.025),  kL_upper95 = quantile(post[, "kL"], 0.975),
      kL_one_sided_upper95 = quantile(post[, "kL"], 0.95),
      Vmax_median  = median(post[, "Vmax"]),
      Vmax_lower95 = quantile(post[, "Vmax"], 0.025), Vmax_upper95 = quantile(post[, "Vmax"], 0.975),
      DO0_median   = median(post[, "DO0"]),
      DO0_lower95  = quantile(post[, "DO0"], 0.025),  DO0_upper95 = quantile(post[, "DO0"], 0.975),
      kL_kg_median    = median(kL_kg),
      kL_kg_lower95   = quantile(kL_kg, 0.025),   kL_kg_upper95   = quantile(kL_kg, 0.975),
      Vmax_kg_median  = median(Vmax_kg),
      Vmax_kg_lower95 = quantile(Vmax_kg, 0.025), Vmax_kg_upper95 = quantile(Vmax_kg, 0.975),
      r_kL_Vmax  = cor(post[, "kL"],   post[, "Vmax"]),
      r_kL_DO0   = cor(post[, "kL"],   post[, "DO0"]),
      r_Vmax_DO0 = cor(post[, "Vmax"], post[, "DO0"]),
      Rhat_kL = rhat_all[["kL"]], Rhat_Vmax = rhat_all[["Vmax"]],
      Rhat_DO0 = rhat_all[["DO0"]], Rhat_sigma2 = rhat_all[["sigma2"]],
      ESS_bulk_min = min(bulk_all), ESS_tail_min = min(tail_all),
      acceptance_rate = mean(sapply(chains, attr, "acceptance")),
      converged = all(rhat_all < rhat_max) && all(bulk_all > ess_min) && all(tail_all > ess_min),
      row.names = NULL))

    # Save every 10th draw (with burn-in, for trace plots)
    keep <- seq(thin_saved, n_iter, by = thin_saved)
    for (ch in 1:4) {
      saved_draws <- rbind(saved_draws, data.frame(
        Sample_Name = s, Moisture = bf$Moisture, Treatment = bf$Treatment, chain = ch,
        iteration = keep, burn_in = keep <= n_burn, chains[[ch]][keep, ],
        kL_kg = chains[[ch]][keep, "kL"] * bf$L_water_per_kg,
        Vmax_kg = chains[[ch]][keep, "Vmax"] * bf$L_water_per_kg))
    }
  }

  write.csv(summary_rows, file.path(mcmc_dir, paste0("YEP_MCMC_Summary_Km", Km, ".csv")), row.names = FALSE)
  saveRDS(saved_draws, file.path(mcmc_dir, paste0("YEP_MCMC_Draws_Km", Km, ".rds")))
}

# ---- 3. Combine all Km files ------------------------------------------------
summary_files <- list.files(mcmc_dir, pattern = "^YEP_MCMC_Summary_Km.*\\.csv$", full.names = TRUE)
mcmc_summary  <- bind_rows(lapply(summary_files, read.csv)) %>% arrange(Km_mg_L, Sample_Name)
write.csv(mcmc_summary, file.path(in_dir, "YEP_MCMC_Summary.csv"), row.names = FALSE)

cat("\nKm values found:", sort(unique(mcmc_summary$Km_mg_L)), "\n")
cat("Fits flagged converged:", sum(mcmc_summary$converged), "of", nrow(mcmc_summary), "\n")
cat("Max R-hat (kL, Vmax, DO0, sigma2):",
    round(max(mcmc_summary[, c("Rhat_kL", "Rhat_Vmax", "Rhat_DO0", "Rhat_sigma2")]), 3),
    "; min bulk ESS:", round(min(mcmc_summary$ESS_bulk_min)),
    "; min tail ESS:", round(min(mcmc_summary$ESS_tail_min)), "\n")
if (any(!mcmc_summary$converged)) {
  cat("Not converged:\n")
  print(mcmc_summary %>% filter(!converged) %>%
          select(Sample_Name, Km_mg_L, Rhat_kL, Rhat_Vmax, Rhat_DO0, ESS_bulk_min, ESS_tail_min))
}

# ---- 4. SI figures (Km = 0.2) -----------------------------------------------
Km_fig <- 0.2
draws_file <- file.path(mcmc_dir, paste0("YEP_MCMC_Draws_Km", Km_fig, ".rds"))
if (!file.exists(draws_file)) stop("Run Km = 0.2 to make the SI figures")
draws <- readRDS(draws_file)

# Same treatment names, colors and shapes as every other figure. The saved
# tables above keep the treatment codes; only the figures show the new names.
source("src/00_plot_style.R")
draws$Treatment        <- relabel_treatment(draws$Treatment)
mcmc_summary$Treatment <- relabel_treatment(mcmc_summary$Treatment)

# (a) Trace plots (all 4 chains, with burn-in) for a typical dry vial, a typical
#     wet vial, and the vial with the widest credible intervals (YEP2A_INC-U1)
trace_vials <- c("YEP1C_INC-S4", "YEP2C_INC-S4", "YEP2A_INC-U1")
trace_data <- draws %>%
  filter(Sample_Name %in% trace_vials) %>%
  pivot_longer(c(kL, Vmax, DO0), names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = factor(Parameter, levels = c("kL", "Vmax", "DO0"),
                            labels = c("kL (h-1)", "Vmax (mg O2 L-1 h-1)", "DO0 (mg L-1)")))
p_trace <- ggplot(trace_data, aes(iteration, Value, color = factor(chain))) +
  annotate("rect", xmin = 0, xmax = n_burn, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_line(linewidth = 0.3) +
  facet_grid(Parameter ~ Sample_Name, scales = "free_y") +
  labs(x = "Iteration", y = NULL, color = "Chain",
       caption = "Shaded: burn-in (proposal tuning; discarded). Every 10th draw shown. Km = 0.2 mg/L.") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "FigS_MCMC_Traces_Km0.2.pdf"), p_trace, width = 10, height = 7)
ggsave(file.path(fig_dir, "FigS_MCMC_Traces_Km0.2.png"), p_trace, width = 10, height = 7, dpi = 300)

# (b) Joint posterior of kL and Vmax (per kg) for every vial: shows how much
#     the two parameters trade off (equifinality) within each vial
post_draws <- draws %>% filter(!burn_in) %>% mutate(Panel = sub("_INC-", " ", Sample_Name))
r_labels <- mcmc_summary %>% filter(Km_mg_L == Km_fig) %>%
  mutate(label = sprintf("r = %.2f", r_kL_Vmax), Panel = sub("_INC-", " ", Sample_Name)) %>%
  select(Panel, label)
p_joint <- ggplot(post_draws, aes(Vmax_kg, kL_kg, color = Treatment, shape = Moisture)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_text(data = r_labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.4, size = 2.8) +
  facet_wrap(~ Panel, scales = "free", ncol = 6) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = moisture_shapes) +
  labs(x = expression(V[max]~(mg~O[2]~kg^{-1}~h^{-1})), y = expression(k[L]~(L~kg^{-1}~h^{-1})),
       color = "DOM treatment", shape = "Inundation history",
       caption = "Posterior draws (after burn-in, every 10th) at Km = 0.2 mg/L; r = posterior correlation of kL and Vmax.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9),
        axis.text = element_text(size = 7), panel.spacing = unit(0.15, "lines"),
        plot.caption = element_text(size = 7)) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
# Saved at the size it is printed in the SI (landscape page), so the text stays readable
ggsave(file.path(fig_dir, "FigS_MCMC_Joint_Posteriors_Km0.2.pdf"), p_joint, width = 9, height = 6.2)
ggsave(file.path(fig_dir, "FigS_MCMC_Joint_Posteriors_Km0.2.png"), p_joint, width = 9, height = 6.2, dpi = 300)

# (c) Posterior median and 95% credible interval for every vial, at each Km in
#     the Sihi range, grouped by sediment and treatment
ci_long <- bind_rows(
  mcmc_summary %>% transmute(Sample_Name, Km_mg_L, Moisture, Treatment, Parameter = "Vmax",
                             median = Vmax_kg_median, lower = Vmax_kg_lower95, upper = Vmax_kg_upper95),
  mcmc_summary %>% transmute(Sample_Name, Km_mg_L, Moisture, Treatment, Parameter = "kL",
                             median = kL_kg_median, lower = kL_kg_lower95, upper = kL_kg_upper95)) %>%
  mutate(Parameter = factor(Parameter, levels = c("Vmax", "kL"),
                            labels = c("V[max]~(mg~O[2]~kg^{-1}~h^{-1})", "k[L]~(L~kg^{-1}~h^{-1})")),
         # order vials by sediment, then treatment, so groups sit together on the y axis
         Sample_Name = factor(Sample_Name, levels = unique(Sample_Name[order(Moisture, Treatment, Sample_Name)])),
         Panel = factor(sub("_INC-", " ", Sample_Name), levels = sub("_INC-", " ", levels(Sample_Name))))   # shorter labels
ci_long$Km_label <- factor(paste0("K[m]==", ci_long$Km_mg_L),
                           levels = paste0("K[m]==", sort(unique(ci_long$Km_mg_L))))
# One page per parameter, so the 30 vial labels on the y axis stay readable
ci_parts <- list(a_Vmax = levels(ci_long$Parameter)[1], b_kL = levels(ci_long$Parameter)[2])
for (part in names(ci_parts)) {
  p_ci <- ggplot(filter(ci_long, Parameter == ci_parts[[part]]),
                 aes(x = median, y = Panel, color = Treatment, shape = Moisture)) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, orientation = "y", linewidth = 0.5) +
    geom_point(size = 1.6) +
    facet_wrap(~ Km_label, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_color_manual(values = treatment_colors) +
    scale_shape_manual(values = moisture_shapes) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    labs(y = NULL, tag = substr(part, 1, 1), color = "DOM treatment", shape = "Inundation history",
         x = parse(text = paste0("Posterior~median~and~95*'%'~credible~interval~of~", ci_parts[[part]]))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 8.5),
          axis.text.x = element_text(size = 8.5), strip.text = element_text(size = 10),
          panel.spacing = unit(0.3, "lines"))
  # Saved at the size it is printed in the SI (landscape page), so the text stays readable
  ggsave(file.path(fig_dir, paste0("FigS3_MCMC_Credible_Intervals_", part, ".pdf")), p_ci, width = 9, height = 6)
  ggsave(file.path(fig_dir, paste0("FigS3_MCMC_Credible_Intervals_", part, ".png")), p_ci, width = 9, height = 6, dpi = 300)
}
