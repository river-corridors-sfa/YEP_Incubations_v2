# =============================================================================
# mcmc_sampler.R -- MCMC sampler shared by 02_model_mcmc.R and
#                   02b_mcmc_checks.R (not run on its own)
# -----------------------------------------------------------------------------
# Model:  dDO/dt = -kL * DO - Vmax * DO / (Km + DO),  Km fixed
#
# Sampled parameters: kL (h-1), Vmax (mg O2 L-1 h-1), DO0 (mg L-1), sigma2.
#   Priors: kL, Vmax : flat (uniform) on [lower, upper]
#           DO0      : normal, centered on the first fitted DO reading with
#                      SD = DO0_prior_sd (sensor uncertainty), and kept within
#                      [lower, upper]. Without this prior, DO0 drifted to
#                      values the probe never recorded (e.g. 12 mg/L when the
#                      first reading was 8) in vials with rho close to 1,
#                      because the AR(1) likelihood then mostly sees the slope
#                      of the curve, not its level.
#           sigma2   : Jeffreys prior (1/sigma2).
#
# Likelihood: the 5-s readings are autocorrelated, so residuals are AR(1):
#     e[t] = rho * e[t-1] + u[t],   u ~ Normal(0, sigma2)
#   rho is FIXED (given as an argument). Readings are 5 s apart (checked in
#   02_model_mcmc.R), so one rho applies to every pair of neighbors.
#
# Sampler, each iteration:
#   (a) Metropolis step for (kL, Vmax, DO0) together, with a multivariate
#       normal proposal. A proposal outside the prior bounds is rejected; the
#       acceptance ratio includes the DO0 prior.
#   (b) Gibbs step for sigma2 (inverse-gamma draw).
#   During burn-in, every 200 iterations, the proposal covariance is re-tuned
#   to the draws so far (adaptive Metropolis, Haario et al. 2001). After
#   burn-in the proposal is fixed, so the retained draws are a valid sample.
# =============================================================================

# Model prediction: DO at each time (NA if the ODE solver fails)
model_DO <- function(kL, Vmax, DO0, Km, times) {
  rhs <- function(t, y, p) list(-kL * max(y, 0) - Vmax * max(y, 0) / (Km + max(y, 0)))
  out <- tryCatch(deSolve::ode(y = c(DO = DO0), times = times, func = rhs, parms = NULL,
                               method = "lsoda", rtol = 1e-6, atol = 1e-8),
                  error = function(e) NULL)
  if (is.null(out) || nrow(out) != length(times)) return(rep(NA, length(times)))
  out[, "DO"]
}

# Run one chain.
#   dat   : data frame with time (h) and DO (mg/L), sorted by time
#   start : named vector c(kL, Vmax, DO0)
#   DO0_prior_sd : SD (mg/L) of the normal prior on DO0, centered on dat$DO[1]
# Returns a matrix (one row per iteration: kL, Vmax, DO0, sigma2) with the
# acceptance rate after burn-in as an attribute.
run_chain <- function(dat, Km, rho, start, lower, upper, n_iter, n_burn, DO0_prior_sd) {
  n <- nrow(dat)

  # Log prior density of DO0 (normal around the first fitted reading)
  log_prior_DO0 <- function(DO0) dnorm(DO0, mean = dat$DO[1], sd = DO0_prior_sd, log = TRUE)

  # Sum of squared AR(1) innovations for one parameter set (Inf if the solver fails)
  ss_innov <- function(p) {
    e <- dat$DO - model_DO(p[["kL"]], p[["Vmax"]], p[["DO0"]], Km, dat$time)
    u <- e[-1] - rho * e[-n]
    if (any(!is.finite(u))) return(Inf)
    sum(u^2)
  }

  cur    <- pmin(pmax(start, lower + 1e-6), upper - 1e-6)   # start inside the prior
  cur_ss <- ss_innov(cur)
  if (!is.finite(cur_ss) || cur_ss <= 0) stop("Invalid starting point for chain")
  sigma2 <- cur_ss / (n - 1)
  prop_cov <- diag((0.02 * abs(cur) + c(0.01, 0.05, 0.01))^2)

  out <- matrix(NA, n_iter, 4, dimnames = list(NULL, c("kL", "Vmax", "DO0", "sigma2")))
  n_acc <- 0
  for (i in 1:n_iter) {
    # (a) Metropolis step for kL, Vmax and DO0
    prop <- cur + as.vector(t(chol(prop_cov)) %*% rnorm(3))
    names(prop) <- names(cur)
    if (all(prop >= lower & prop <= upper)) {
      prop_ss <- ss_innov(prop)
      log_ratio <- (cur_ss - prop_ss) / (2 * sigma2) +                         # likelihood
                   log_prior_DO0(prop[["DO0"]]) - log_prior_DO0(cur[["DO0"]])   # DO0 prior
      if (log(runif(1)) < log_ratio) {
        cur <- prop; cur_ss <- prop_ss
        if (i > n_burn) n_acc <- n_acc + 1
      }
    }
    # (b) Gibbs step for sigma2
    sigma2 <- 1 / rgamma(1, shape = (n - 1) / 2, rate = cur_ss / 2)
    out[i, ] <- c(cur, sigma2)

    # (c) Burn-in only: re-tune the proposal to the recent draws
    if (i <= n_burn && i %% 200 == 0 && i >= 400) {
      prop_cov <- (2.38^2 / 3) * cov(out[(i / 2):i, 1:3]) + diag(1e-10, 3)
    }
  }
  attr(out, "acceptance") <- n_acc / (n_iter - n_burn)
  out
}

# Four starting points spread around the best fit, with different balances of
# kL and Vmax (kL-heavy, Vmax-heavy), so the chains test for other solutions.
dispersed_starts <- function(kL, Vmax, DO0) {
  list(c(kL = kL,           Vmax = Vmax,       DO0 = DO0),
       c(kL = 2 * kL + 0.5, Vmax = 0.6 * Vmax, DO0 = DO0),
       c(kL = 0.3 * kL,     Vmax = 1.3 * Vmax, DO0 = DO0),
       c(kL = 1.5 * kL + 1, Vmax = 0.8 * Vmax, DO0 = DO0 + 0.1))
}
