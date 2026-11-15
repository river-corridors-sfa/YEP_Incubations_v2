# =============================================================================
# 03_stats.R -- statistics for the figures and text
# -----------------------------------------------------------------------------
# All tests are non-parametric (rank-based), so a single extreme microcosm
# cannot drive a result.
#
# Part A. Dry vs wet at every Km (stars in Fig. 1 and Fig. S Km full range)
#   Mann-Whitney (Wilcoxon rank-sum) test of dry vs wet, treatments pooled
#   (15 vs 15 vials), for kL and Vmax per kg (leas-square estimates 01_model_DO.R).
#
# Part B. DOM treatment effects (letters in Fig. 1 and Fig. S Km full range)
#   Within each sediment history (dry, wet), n = 5 vials per treatment:
#     1. Kruskal-Wallis test among Control, Unburned DOM and Burned DOM
#     2. Dunn's test for each pair of treatments, with no p-value adjustment.
#        Saved for every cell (Dunn_Treatment_byKm.csv, with the Kruskal-Wallis
#        p), but a pair is called significant, and letters are drawn, ONLY
#        if Kruskal-Wallis p < 0.05
#   Responses: Vmax and kL per kg at every Km, and the measured O2 consumption
#   rate per kg (does not depend on Km).
#
# Part C. Two more checks, with the same tests
#   C1. Measured O2 consumption rate, dry vs wet (Mann-Whitney, treatments
#       pooled). Part A covers only kL and Vmax.
#   C2. Do low-DO readings drive the results? kL and Vmax refit at Km = 0.2
#       using only readings with DO >= 2 mg/L (01_model_DO.R), tested again:
#       dry vs wet (Mann-Whitney) and DOM treatments within dry and within wet
#       (Kruskal-Wallis, then Dunn if significant). Each row is shown next to
#       the same test on the fits that use all readings.
#
# Part D. H2: is the burned - control difference in Vmax larger in previously
#   dry than in wet sediments? Done at every Km. (kL is not tested because it
#   is 0 in nearly every dry vial, so the dry burned - control difference is 0.)
#     Burned - control difference = median Vmax of the 5 burned vials
#                                   - median Vmax of the 5 control vials
#     Interaction = (dry burned - dry control) - (wet burned - wet control)
#     > 0: difference larger in dry (supports H2); < 0: larger in wet.
#   95% percentile bootstrap CI: 10,000 resamples, each drawing 5 vials with
#   replacement within each of the 4 groups. A CI that includes 0 means no
#   evidence that the difference differs between histories.
#
# The linear-model version of Part B (reviewed in November 2026) is in
# Archive/2026-11_treatment_linear_model_reviewed/.
#
# Inputs : data_files$rates, data_files$npoc_tn, data_files$ions,
#          outputs/modeling/YEP_BestFit_AllKm.csv, YEP_BestFit_DOabove2_Km0.2.csv
# Outputs: outputs/stats/MannWhitney_DryWet_byKm.csv
#          outputs/stats/microcosm_responses.csv
#          outputs/stats/KruskalWallis_Treatment_byKm.csv
#          outputs/stats/Dunn_Treatment_byKm.csv
#          outputs/stats/Treatment_Letters_byKm.csv
#          outputs/stats/MannWhitney_DryWet_TotalO2.csv
#          outputs/stats/DOabove2_Refit_Tests.csv
#          outputs/stats/DOabove2_Refit_Dunn.csv (only if a Kruskal-Wallis test is significant)
#          outputs/stats/Bootstrap_Interaction_Vmax_byKm.csv
# =============================================================================

rm(list = ls(all = TRUE))
library(dplyr)
library(tidyr)     # pivot_longer()
library(rstatix)   # dunn_test()
source("src/00_plot_style.R")   # treatment names shown in tables and figures
source("src/00_data_paths.R")

Km_base <- 0.2
out_dir <- file.path("outputs", "stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

best_fits <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_AllKm.csv"))

# =============================================================================
# Part A. Mann-Whitney dry vs wet at every Km
# =============================================================================
mw <- data.frame()
for (Km in sort(unique(best_fits$Km_mg_L))) {
  for (param in c("kL_L_kg_h", "Vmax_mg_kg_h")) {
    d   <- best_fits %>% filter(Km_mg_L == Km)
    dry <- d[[param]][d$Moisture == "Dry"]
    wet <- d[[param]][d$Moisture == "Wet"]
    test <- wilcox.test(dry, wet, exact = FALSE)   # exact = FALSE: dry kL has ties at 0
    mw <- rbind(mw, data.frame(Km_mg_L = Km, Parameter = param,
                               n_dry = length(dry), n_wet = length(wet),
                               median_dry = median(dry), median_wet = median(wet),
                               W = unname(test$statistic), p = test$p.value))
  }
}
mw$stars <- as.character(cut(mw$p, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns")))
mw$stars[is.na(mw$p)] <- ""   # no test possible: e.g. kL = 0 in every vial at Km = 10
write.csv(mw, file.path(out_dir, "MannWhitney_DryWet_byKm.csv"), row.names = FALSE)

# =============================================================================
# Part B. DOM treatment effects
# =============================================================================

# ---- B1. One row per microcosm ----------------------------------------------
# A microcosm has different sample names per analysis (INC = O2, SOC = NPOC,
# SIN = ions). The shared ID drops that code: YEP1A_INC-H2 -> YEP1A_H2.
make_id <- function(x) sub("^(YEP[12][ABC])_[A-Z]+-([HSU][0-9])$", "\\1_\\2", x)
is_microcosm <- function(x) grepl("^YEP[12][ABC]_[A-Z]+-[HSU][0-9]$", x)

# Flagged cells look like "TN_Above_...|3.3_ppm_Raw_Not_Corrected|23.09_ppm_Final_Corrected";
# use the Final_Corrected number, otherwise the cell itself.
flag_value <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("Final_Corrected", x), sub(".*\\|([0-9.]+)_ppm_Final_Corrected.*", "\\1", x), x)
  suppressWarnings(as.numeric(x))
}

rates <- read_data(data_files$rates) %>%
  filter(is_microcosm(Sample_Name)) %>%
  transmute(ID = make_id(Sample_Name),
            O2_rate_mg_kg_h = -as.numeric(Respiration_Rate_mg_DO_per_kg_per_H))   # consumption as positive

params <- best_fits %>%
  filter(Km_mg_L == Km_base) %>%
  transmute(ID = make_id(Sample_Name), Vmax_mg_kg_h, kL_L_kg_h)

params_DO2 <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_DOabove2_Km0.2.csv")) %>%
  transmute(ID = make_id(Sample_Name), Vmax_DO2 = Vmax_mg_kg_h, kL_DO2 = kL_L_kg_h)

doc <- read_data(data_files$npoc_tn) %>%
  filter(is_microcosm(Sample_Name)) %>%
  transmute(ID = make_id(Sample_Name), DOC_mg_kg = flag_value(Extractable_NPOC_mg_per_kg))

ions <- read_data(data_files$ions) %>%
  filter(is_microcosm(Sample_Name)) %>%
  transmute(ID = make_id(Sample_Name), NO3_mg_N_L = flag_value(`00618_NO3_mg_per_L_as_N`))

dat <- rates %>%
  full_join(params, by = "ID") %>%
  full_join(params_DO2, by = "ID") %>%
  full_join(doc, by = "ID") %>%
  full_join(ions, by = "ID") %>%
  mutate(Moisture  = factor(ifelse(substr(ID, 4, 4) == "1", "Dry", "Wet"), levels = c("Dry", "Wet")),
         Treatment = factor(dplyr::recode(sub(".*_([HSU])[0-9]$", "\\1", ID),
                                          S = "Control", U = "Unburned + DOC", H = "High Burn + DOC"),
                            levels = c("Control", "Unburned + DOC", "High Burn + DOC")),
         Day       = factor(substr(ID, 5, 5), levels = c("A", "B", "C"), labels = c("1", "2", "3")))

stopifnot(nrow(dat) == 30, all(table(dat$Moisture, dat$Treatment) == 5))
write.csv(dat, file.path(out_dir, "microcosm_responses.csv"), row.names = FALSE)

# ---- B2. Long table: one row per microcosm x response x Km -------------------
# Measured O2 consumption has no Km (Km_mg_L = NA); Vmax and kL have one value
# per Km.
o2_long <- dat %>%
  transmute(ID, Moisture, Treatment, Km_mg_L = NA_real_,
            Response = "O2_rate_mg_kg_h", Value = O2_rate_mg_kg_h)

param_long <- best_fits %>%
  transmute(ID = make_id(Sample_Name), Km_mg_L, Vmax_mg_kg_h, kL_L_kg_h) %>%
  left_join(dat %>% select(ID, Moisture, Treatment), by = "ID") %>%
  pivot_longer(c(Vmax_mg_kg_h, kL_L_kg_h), names_to = "Response", values_to = "Value")

long <- bind_rows(o2_long, param_long) %>%
  mutate(Treatment = relabel_treatment(Treatment))   # names used in the manuscript
treatments <- levels(long$Treatment)

# ---- B3. Letters that summarize the Dunn comparisons -------------------------
# Treatments that share a letter are NOT significantly different.
# groups: treatment names in order; sig_pairs: list of c(treatment1, treatment2)
letters_from_pairs <- function(groups, sig_pairs) {
  sets <- list(groups)                       # start: everything shares one letter
  for (pair in sig_pairs) {
    new_sets <- list()
    for (s in sets) {
      if (all(pair %in% s)) {                # this letter joins a significant pair: split it
        new_sets <- c(new_sets, list(setdiff(s, pair[1])), list(setdiff(s, pair[2])))
      } else {
        new_sets <- c(new_sets, list(s))
      }
    }
    # drop letters that are contained in another letter (or are duplicates)
    keep <- rep(TRUE, length(new_sets))
    for (i in seq_along(new_sets)) for (j in seq_along(new_sets)) {
      if (i != j && all(new_sets[[i]] %in% new_sets[[j]]) &&
          (length(new_sets[[i]]) < length(new_sets[[j]]) || j < i)) keep[i] <- FALSE
    }
    sets <- new_sets[keep]
  }
  sets <- sets[order(sapply(sets, function(s) min(match(s, groups))))]   # "a" starts at the control
  sapply(groups, function(g) paste(letters[which(sapply(sets, function(s) g %in% s))], collapse = ""))
}

# ---- B4. Kruskal-Wallis, then Dunn if significant ----------------------------
alpha <- 0.05
kw_all <- data.frame(); dunn_all <- data.frame(); letters_all <- data.frame()

tests <- long %>% distinct(Response, Km_mg_L, Moisture) %>% arrange(Response, Km_mg_L, Moisture)
for (k in seq_len(nrow(tests))) {
  resp <- tests$Response[k]; Km <- tests$Km_mg_L[k]; mo <- tests$Moisture[k]
  d <- long %>% filter(Response == resp, Moisture == mo,
                       if (is.na(Km)) is.na(Km_mg_L) else Km_mg_L == Km)
  stopifnot(nrow(d) == 15, all(table(d$Treatment) == 5))
  medians <- tapply(d$Value, d$Treatment, median)

  if (length(unique(d$Value)) == 1) {
    # every vial has the same value (kL = 0 in all dry vials at some Km): nothing to test
    H <- NA; p <- NA; note <- "all values identical, not tested"
  } else {
    kw <- kruskal.test(Value ~ Treatment, data = d)
    H <- unname(kw$statistic); p <- kw$p.value; note <- ""
  }
  significant <- !is.na(p) && p < alpha

  kw_all <- rbind(kw_all, data.frame(
    Response = resp, Km_mg_L = Km, Moisture = mo, n_per_treatment = 5,
    median_Control = medians[[treatments[1]]], median_Unburned = medians[[treatments[2]]],
    median_Burned = medians[[treatments[3]]],
    H = H, df = 2, p = p, significant = significant, note = note))

  # Dunn's test is run and saved for every testable cell, so every pairwise
  # p value cited in the text can be checked. It is INTERPRETED (letters,
  # "significant" column) only when Kruskal-Wallis p < 0.05.
  if (!is.na(p)) {
    dn <- dunn_test(d, Value ~ Treatment, p.adjust.method = "none")   # two-sided, no adjustment
    dunn_all <- rbind(dunn_all, data.frame(
      Response = resp, Km_mg_L = Km, Moisture = mo,
      Treatment_1 = dn$group1, Treatment_2 = dn$group2, Z = dn$statistic, p = dn$p,
      KW_p = p, KW_significant = significant,
      significant = significant & dn$p < alpha))   # pair counts only if Kruskal-Wallis is significant
  }
  if (significant) {
    sig <- dn[dn$p < alpha, ]
    sig_pairs <- lapply(seq_len(nrow(sig)), function(i) c(sig$group1[i], sig$group2[i]))
    letters_all <- rbind(letters_all, data.frame(
      Response = resp, Km_mg_L = Km, Moisture = mo, Treatment = treatments,
      Letter = unname(letters_from_pairs(treatments, sig_pairs))))
  }
}

write.csv(kw_all,      file.path(out_dir, "KruskalWallis_Treatment_byKm.csv"), row.names = FALSE)
write.csv(dunn_all,    file.path(out_dir, "Dunn_Treatment_byKm.csv"), row.names = FALSE)
write.csv(letters_all, file.path(out_dir, "Treatment_Letters_byKm.csv"), row.names = FALSE)

# =============================================================================
# Part C. Two more checks
# =============================================================================

# ---- C1. Measured O2 consumption rate: dry vs wet ----------------------------
o2_dry <- dat$O2_rate_mg_kg_h[dat$Moisture == "Dry"]
o2_wet <- dat$O2_rate_mg_kg_h[dat$Moisture == "Wet"]
o2_test <- wilcox.test(o2_dry, o2_wet, exact = FALSE)
mw_o2 <- data.frame(Response = "O2_rate_mg_kg_h", n_dry = length(o2_dry), n_wet = length(o2_wet),
                    median_dry = median(o2_dry), median_wet = median(o2_wet),
                    W = unname(o2_test$statistic), p = o2_test$p.value)
write.csv(mw_o2, file.path(out_dir, "MannWhitney_DryWet_TotalO2.csv"), row.names = FALSE)

# ---- C2. Refits with DO >= 2 mg/L only (Km = 0.2) ----------------------------
# One row per test. "all_readings" = the standard fit at Km = 0.2;
# "DO_above_2" = the refit that leaves out readings below 2 mg/L.
refit_columns <- list(Vmax_mg_kg_h = c(all_readings = "Vmax_mg_kg_h", DO_above_2 = "Vmax_DO2"),
                      kL_L_kg_h    = c(all_readings = "kL_L_kg_h",    DO_above_2 = "kL_DO2"))
dat_c <- dat %>% mutate(Treatment = relabel_treatment(Treatment))

# p value of a test, or NA when every value is identical (nothing to test)
safe_p <- function(values, run_test) if (length(unique(values)) == 1) NA else run_test()$p.value

refit_tests <- data.frame(); refit_dunn <- data.frame()
for (param in names(refit_columns)) {
  cols <- refit_columns[[param]]

  # Dry vs wet, treatments pooled
  p_mw <- sapply(cols, function(col) safe_p(dat_c[[col]], function()
    wilcox.test(dat_c[[col]][dat_c$Moisture == "Dry"], dat_c[[col]][dat_c$Moisture == "Wet"], exact = FALSE)))
  refit_tests <- rbind(refit_tests, data.frame(
    Parameter = param, Test = "Mann-Whitney, dry vs wet (treatments pooled)", Moisture = "Dry vs Wet",
    p_all_readings = p_mw[["all_readings"]], p_DO_above_2 = p_mw[["DO_above_2"]]))

  # DOM treatments within dry and within wet
  for (mo in c("Dry", "Wet")) {
    d <- dat_c[dat_c$Moisture == mo, ]
    p_kw <- sapply(cols, function(col) safe_p(d[[col]], function() kruskal.test(d[[col]] ~ d$Treatment)))
    refit_tests <- rbind(refit_tests, data.frame(
      Parameter = param, Test = "Kruskal-Wallis, DOM treatments", Moisture = mo,
      p_all_readings = p_kw[["all_readings"]], p_DO_above_2 = p_kw[["DO_above_2"]]))

    if (!is.na(p_kw[["DO_above_2"]]) && p_kw[["DO_above_2"]] < alpha) {
      d$Value <- d[[cols[["DO_above_2"]]]]
      dn <- dunn_test(d, Value ~ Treatment, p.adjust.method = "none")
      refit_dunn <- rbind(refit_dunn, data.frame(
        Parameter = param, Moisture = mo, Treatment_1 = dn$group1, Treatment_2 = dn$group2,
        Z = dn$statistic, p = dn$p, significant = dn$p < alpha))
    }
  }
}
refit_tests$changed_at_0.05 <- (refit_tests$p_all_readings < alpha) != (refit_tests$p_DO_above_2 < alpha)
write.csv(refit_tests, file.path(out_dir, "DOabove2_Refit_Tests.csv"), row.names = FALSE)
if (nrow(refit_dunn) > 0) write.csv(refit_dunn, file.path(out_dir, "DOabove2_Refit_Dunn.csv"), row.names = FALSE)

# =============================================================================
# Part D. H2: burned - control difference in Vmax, dry vs wet (bootstrap), every Km
# =============================================================================
n_boot <- 10000

# Burned - control differences and the interaction, from the medians of the
# four groups of 5 vials
differences <- function(dc, db, wc, wb) {
  m <- sapply(list(dc, db, wc, wb), median)
  c(dry = m[2] - m[1], wet = m[4] - m[3], interaction = (m[2] - m[1]) - (m[4] - m[3]))
}

h2 <- data.frame()
for (Km in sort(unique(best_fits$Km_mg_L))) {
  set.seed(2026)   # same resamples every run, and at every Km
  d <- best_fits %>% filter(Km_mg_L == Km)
  group_vmax <- function(mo, tr) d$Vmax_mg_kg_h[d$Moisture == mo & d$Treatment == tr]
  DC <- group_vmax("Dry", "Control"); DB <- group_vmax("Dry", "High Burn + DOC")
  WC <- group_vmax("Wet", "Control"); WB <- group_vmax("Wet", "High Burn + DOC")
  stopifnot(all(lengths(list(DC, DB, WC, WB)) == 5))

  obs  <- differences(DC, DB, WC, WB)
  boot <- t(replicate(n_boot, differences(sample(DC, replace = TRUE), sample(DB, replace = TRUE),
                                          sample(WC, replace = TRUE), sample(WB, replace = TRUE))))
  ci <- apply(boot, 2, quantile, c(0.025, 0.975))

  h2 <- rbind(h2, data.frame(
    Response = "Vmax_mg_kg_h", Km_mg_L = Km, Comparison = "Burned DOM - Control",
    median_dry_control = median(DC), median_dry_burned = median(DB),
    median_wet_control = median(WC), median_wet_burned = median(WB),
    difference_dry = obs[["dry"]], difference_wet = obs[["wet"]], interaction = obs[["interaction"]],
    CI_low = ci[1, "interaction"], CI_high = ci[2, "interaction"], n_boot = n_boot))
}
h2$CI_excludes_0 <- h2$CI_low > 0 | h2$CI_high < 0
write.csv(h2, file.path(out_dir, "Bootstrap_Interaction_Vmax_byKm.csv"), row.names = FALSE)

# ---- Console summary --------------------------------------------------------
options(width = 150)
sihi <- function(x) is.na(x) | (x >= 0.13 & x <= 1.22)   # Km values reported in the manuscript
cat("\nMann-Whitney dry vs wet, Km in the Sihi range and Km = 10:\n")
print(mw %>% filter(Km_mg_L >= 0.13 & Km_mg_L <= 1.22 | Km_mg_L == 10) %>%
        mutate(across(c(median_dry, median_wet), ~ round(., 2)), p = signif(p, 2)), row.names = FALSE)
cat("\nKruskal-Wallis among DOM treatments (measured O2 rate, and Km in the Sihi range):\n")
print(kw_all %>% filter(sihi(Km_mg_L)) %>%
        mutate(across(starts_with("median"), ~ round(., 2)), H = round(H, 2), p = round(p, 3)) %>%
        select(-n_per_treatment, -df), row.names = FALSE)
cat("\nDunn's test where Kruskal-Wallis p < 0.05 (no p-value adjustment):\n")
print(dunn_all %>% filter(sihi(Km_mg_L), KW_significant) %>% select(-KW_p, -KW_significant) %>%
        mutate(Z = round(Z, 2), p = round(p, 3)), row.names = FALSE)
cat("\nDunn's test where Kruskal-Wallis p >= 0.05 (saved, not interpreted): smallest pairwise p by response and history\n")
print(dunn_all %>% filter(sihi(Km_mg_L), !KW_significant) %>% group_by(Response, Moisture) %>%
        summarize(min_pair_p = round(min(p), 3), max_pair_p = round(max(p), 3), .groups = "drop"), row.names = FALSE)
cat("\nMeasured O2 consumption rate, dry vs wet (Mann-Whitney):\n")
print(mw_o2 %>% mutate(across(starts_with("median"), ~ round(., 1)), p = signif(p, 2)), row.names = FALSE)
cat("\nRefits with DO >= 2 mg/L only (Km = 0.2) compared with all readings:\n")
print(refit_tests %>% mutate(across(starts_with("p_"), ~ signif(., 2))), row.names = FALSE)
cat("\nH2: burned - control difference in Vmax (mg O2 kg-1 h-1), previously dry vs wet (bootstrap 95% CI):\n")
print(h2 %>% select(Km_mg_L, difference_dry, difference_wet, interaction, CI_low, CI_high, CI_excludes_0) %>%
        mutate(across(where(is.numeric), ~ round(., 2))), row.names = FALSE)
