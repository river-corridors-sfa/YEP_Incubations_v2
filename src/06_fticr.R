# =============================================================================
# 06_fticr.R -- FTICR-MS: NMDS, Wd* tests, dispersion, AI_mod and NOSC (Fig. 3)
# -----------------------------------------------------------------------------
# Post-incubation aqueous samples (SIR). All analyses are intensity weighted:
# peak intensities are converted to relative intensities within each sample
# (each sample sums to 1).
#   - Bray-Curtis dissimilarity on Hellinger-transformed relative intensities
#     (square root of relative intensity)
#   - NMDS on Bray-Curtis dissimilarity (Fig. 3)
#   - group differences in composition (location/centroid): distance-based
#     multivariate Welch ANOVA, Wd* test (Hamidi et al. 2019, Microbiome 7:51),
#     which is robust to unequal dispersion and unequal group sizes. Used
#     instead of PERMANOVA because wet-sediment treatments differ in dispersion.
#       * inundation history (all samples)
#       * DOM treatment within previously dry and within wet sediments
#       permutations restricted within incubation day (strata)
#     followed by pairwise distance-based Welch tests (Tw2; Alekseyenko 2016,
#     Bioinformatics 32:3552-3558), unadjusted and interpreted only when the
#     Wd* test for that inundation history is significant (same approach as
#     Dunn's tests after Kruskal-Wallis)
#   - dispersion (betadisper / PERMDISP), reported alongside the Wd* tests
#   - mean Bray-Curtis dissimilarity of each wet treatment to the previously dry
#     samples (used in the Results text)
#   - intensity-weighted mean AI_mod and NOSC per sample (table)
#   - Figure 3: (a) AI_mod and (b) NOSC by DOM treatment within each
#     history (Kruskal-Wallis, Dunn letters), (c) NMDS
# Note: YEP2A_SIR-S2 and YEP2B_SIR-S3 were reacquired at 0.2 s ion
#   accumulation (low signal at 0.075 s); the processed data package contains
#   only the 0.2 s spectra for these two samples, so no special handling here.
# Formula filters (applied before all analyses; see "Filter Formulas"):
#   1) blank removal: formulas detected in > 50% of the methanol SPE process
#      blanks are removed (data package readme, note 11);
#   2) mass window: only formulas with 200 <= m/z <= 900 are kept (the
#      processed data package was filtered at 100-1000 m/z).
#
# Inputs : data_files$fticr_data, data_files$fticr_mol, data_files$sample_metadata
# Outputs: outputs/figures/Fig3_FTICR.pdf/.png  (manuscript Figure 3)
#          outputs/fticr/fticr_formula_filter_log.csv, fticr_formula_filter_per_sample.csv,
#          fticr_sample_counts.csv, fticr_nmds_site_scores.csv, fticr_nmds_stress.csv,
#          fticr_wds_tests.csv, fticr_tw2_pairwise.csv, fticr_betadisper.csv,
#          fticr_mean_distance_to_dry.csv, fticr_chemistry_per_sample.csv,
#          fticr_chemistry_by_group.csv, fticr_chemistry_tests.csv
# =============================================================================

library(tidyverse)
library(vegan)
library(WdStar)   # Wd* and Tw2 tests; install with remotes::install_github("alekseyenko/WdStar")

# Setup ----
set.seed(123)

# File paths (run from the repo root)
source("src/00_data_paths.R")
source("src/00_plot_style.R")   # treatment names, colors, shapes used in every figure
metadata_file <- data_files$sample_metadata
fticr_file    <- data_files$fticr_data
mol_file      <- data_files$fticr_mol

# Output folders
out_dir <- file.path("outputs", "fticr")
fig_dir <- file.path("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Treatment codes used in the data and in the output tables. The figure shows
# the manuscript names from src/00_plot_style.R (see "Create NMDS Plot").
treatment_codes <- names(treatment_labels)

# Target samples: post-incubation aqueous samples
target_sample_type <- "Aqueous sample post-incubation (treatment solution and sediment)"


# Load and Clean Data ----

## Read metadata and standardize treatment/moisture names
metadata <- read_csv(metadata_file, show_col_types = FALSE, trim_ws = TRUE) %>%
  mutate(
    # Simplify treatment names
    Treatment = case_when(
      str_detect(Treatment, "Columbia River") ~ "Control",
      str_detect(Treatment, "Unburned") ~ "Unburned + DOC",
      str_detect(Treatment, "Burned") ~ "High Burn + DOC"
    ),
    # Standardize moisture categories
    Moisture = str_to_title(str_squish(Sediment_Field_Moisture_Conditions)),
    # Convert to factors with defined order
    Treatment = factor(Treatment, levels = treatment_codes),
    Moisture = factor(Moisture, levels = names(moisture_shapes)),
    Incubation_Day = factor(Incubation_Day, levels = c("Day_1", "Day_2", "Day_3"))
  )

## Read FTICR data (rows = molecules, columns = samples)
fticr_raw <- read_csv(fticr_file, show_col_types = FALSE, trim_ws = TRUE)

## Read molecular formula information (rows match the Data file)
mol <- read_csv(mol_file, show_col_types = FALSE)
stopifnot(nrow(mol) == nrow(fticr_raw), all(abs(mol$Calibrated_Mass - fticr_raw$Calibrated_Mass) < 1e-6))


# Filter Formulas ----
# 1) Blank removal: remove formulas detected in more than half of the methanol
#    SPE process blanks (>= 3 of 4), as recommended in the data package readme
#    (note 11). No blank subtraction was done in the processed data package.
# 2) Mass window: keep formulas with 200 <= m/z <= 900, the reliable m/z
#    window. The processed data package was filtered at 100-1000 m/z.

## Blank columns and number of blanks in which each formula was detected
blank_cols <- grep("MeOH_blank", names(fticr_raw), value = TRUE)
stopifnot(length(blank_cols) == 4)
n_blank_detected <- rowSums(fticr_raw[blank_cols] > 0, na.rm = TRUE)
in_blanks <- n_blank_detected > length(blank_cols) / 2

## Mass window
in_window <- fticr_raw$Calibrated_Mass >= 200 & fticr_raw$Calibrated_Mass <= 900

## Formulas kept: not in blanks and inside the mass window
keep <- !in_blanks & in_window

## Log of formulas removed by each filter (whole dataset)
filter_log <- tibble(
  Step = c("Formulas in processed data package",
           "Detected in > 50% of methanol blanks",
           "Outside 200-900 m/z",
           "Removed by either filter",
           "Formulas kept"),
  n_formulas = c(nrow(fticr_raw), sum(in_blanks), sum(!in_window),
                 sum(!keep), sum(keep)))
print(filter_log)
write_csv(filter_log, file.path(out_dir, "fticr_formula_filter_log.csv"))

## Formulas per sample before and after filtering
sample_cols_all <- setdiff(names(fticr_raw), c("Calibrated_Mass", blank_cols))
filter_per_sample <- tibble(
  Sample_Name = sample_cols_all,
  n_formulas_before = colSums(fticr_raw[sample_cols_all] > 0, na.rm = TRUE),
  n_removed_blanks = colSums(fticr_raw[in_blanks, sample_cols_all] > 0, na.rm = TRUE),
  n_removed_mass = colSums(fticr_raw[!in_window, sample_cols_all] > 0, na.rm = TRUE),
  n_formulas_after = colSums(fticr_raw[keep, sample_cols_all] > 0, na.rm = TRUE)) %>%
  mutate(pct_removed = round(100 * (1 - n_formulas_after / n_formulas_before), 1))
print(filter_per_sample, n = Inf)
write_csv(filter_per_sample, file.path(out_dir, "fticr_formula_filter_per_sample.csv"))

## Apply the filter to both files (keeps their rows aligned) and drop the blanks
fticr_raw <- fticr_raw[keep, ] %>% select(-all_of(blank_cols))
mol <- mol[keep, ]
stopifnot(nrow(mol) == nrow(fticr_raw), all(abs(mol$Calibrated_Mass - fticr_raw$Calibrated_Mass) < 1e-6))


# Create Intensity Matrix ----

## Extract sample columns (exclude mass column)
sample_cols <- setdiff(names(fticr_raw), "Calibrated_Mass")

## Peak intensities (0 = not detected); rows = samples, columns = formulas
int_matrix <- fticr_raw %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  as.matrix() %>%
  t()

rownames(int_matrix) <- sample_cols
colnames(int_matrix) <- paste0("mol_", seq_len(nrow(fticr_raw)))


# Filter Samples ----

## Keep only target samples with complete metadata
sample_info <- tibble(Sample_Name = rownames(int_matrix)) %>%
  left_join(metadata, by = "Sample_Name") %>%
  filter(
    Sample_Material == target_sample_type,
    !is.na(Treatment), !is.na(Moisture), !is.na(Incubation_Day)
  )

## Filter intensity matrix to match
int_matrix <- int_matrix[sample_info$Sample_Name, , drop = FALSE]

## Keep formula columns aligned with mol (formulas absent from every target
## sample have zero relative intensity and do not affect any result)
stopifnot(all(rowSums(int_matrix) > 0))
stopifnot(identical(sample_info$Sample_Name, rownames(int_matrix)))

## Relative intensity within each sample (each row sums to 1)
rel_matrix <- int_matrix / rowSums(int_matrix)

## Hellinger transformation (square root of relative intensity), which
## down-weights the most intense peaks before computing Bray-Curtis
hel_matrix <- sqrt(rel_matrix)

## Number of samples per group (Fig. 3 caption)
sample_counts <- sample_info %>% count(Moisture, Treatment, name = "n_samples")
print(sample_counts)
write_csv(sample_counts, file.path(out_dir, "fticr_sample_counts.csv"))


# Run NMDS Ordination ----

## Bray-Curtis dissimilarity on Hellinger-transformed relative intensities
dist_bc <- vegdist(hel_matrix, method = "bray")

## Run NMDS to create 2D representation
set.seed(123)
nmds <- metaMDS(dist_bc, k = 2, trymax = 200, trace = FALSE)

## Extract sample coordinates and add metadata
site_scores <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_info, by = "Sample_Name")

write_csv(site_scores, file.path(out_dir, "fticr_nmds_site_scores.csv"))

## NMDS stress and run details (the stress printed on Fig. 3c and in its caption)
nmds_stress <- tibble(
  Distance = "Bray-Curtis on Hellinger-transformed relative intensities",
  k = nmds$ndim, n_samples = nrow(hel_matrix), n_formulas = ncol(hel_matrix),
  Stress = nmds$stress, Stress_rounded = round(nmds$stress, 2),
  Converged = nmds$converged, Tries = nmds$tries, Seed = 123)
print(nmds_stress)
write_csv(nmds_stress, file.path(out_dir, "fticr_nmds_stress.csv"))


# Group differences in composition: Wd* tests ----
# PERMANOVA's pseudo-F can be significant when groups differ only in spread
# (dispersion), and the wet-sediment treatments differ in spread (see
# betadisper below). The Wd* test compares group centroids while accounting for
# unequal dispersion and unequal group sizes. Permutations are restricted
# within incubation day (strata) so that day differences do not drive the
# result. Tests:
#   1) inundation history, all samples
#   2) DOM treatment within previously dry sediments
#   3) DOM treatment within wet sediments
n_perm <- 9999
bc_mat <- as.matrix(dist_bc)

run_wds <- function(samples, group_var, label) {
  s_info <- sample_info %>% filter(Sample_Name %in% samples)
  dm <- as.dist(bc_mat[s_info$Sample_Name, s_info$Sample_Name])
  set.seed(123)
  w <- WdS.test(dm, group = droplevels(s_info[[group_var]]), nrep = n_perm,
                strata = droplevels(s_info$Incubation_Day))
  tibble(Test = label, Grouping = group_var, n_samples = nrow(s_info),
         n_groups = nlevels(droplevels(s_info[[group_var]])),
         WdS = unname(w$statistic), omega_squared = unname(w$estimate),
         p = w$p.value, n_perm = n_perm, strata = "Incubation_Day")
}

dry_samples <- sample_info$Sample_Name[sample_info$Moisture == "Dry"]
wet_samples <- sample_info$Sample_Name[sample_info$Moisture == "Wet"]

wds_tests <- bind_rows(
  run_wds(sample_info$Sample_Name, "Moisture", "Inundation history (all samples)"),
  run_wds(dry_samples, "Treatment", "DOM treatment, previously dry"),
  run_wds(wet_samples, "Treatment", "DOM treatment, wet"))
print(wds_tests)
write_csv(wds_tests, file.path(out_dir, "fticr_wds_tests.csv"))

# Pairwise comparisons among DOM treatments: distance-based Welch test (Tw2)
# for every pair of treatments within each inundation history. As for Dunn's
# tests elsewhere in the paper, p values are not adjusted and pairs are
# interpreted only when the Wd* test for that history is significant
# (p < 0.05); with three groups this protected approach controls the overall
# error rate. Note: Tw2 permutations are not stratified by day (the function
# does not support strata). With 4-5 samples per group, the smallest possible
# Tw2 p value is about 0.008.
tw2_pairwise <- tibble()
for (mo in c("Dry", "Wet")) {
  s_info <- sample_info %>% filter(Moisture == mo)
  trt <- droplevels(s_info$Treatment)
  wds_p <- wds_tests$p[wds_tests$Test == paste0("DOM treatment, ",
                         ifelse(mo == "Dry", "previously dry", "wet"))]
  for (pair in combn(levels(trt), 2, simplify = FALSE)) {
    keep_s <- trt %in% pair
    dm <- as.dist(bc_mat[s_info$Sample_Name[keep_s], s_info$Sample_Name[keep_s]])
    set.seed(123)
    tw <- Tw2.test(dm, droplevels(trt[keep_s]), nrep = n_perm)
    tw2_pairwise <- bind_rows(tw2_pairwise, tibble(
      Moisture = mo, Treatment_1 = pair[1], Treatment_2 = pair[2],
      n_1 = sum(trt == pair[1]), n_2 = sum(trt == pair[2]),
      Tw2 = unname(tw$statistic), p = tw$p.value, n_perm = n_perm,
      WdS_p = wds_p, significant = wds_p < 0.05 & tw$p.value < 0.05))
  }
}
print(tw2_pairwise)
write_csv(tw2_pairwise, file.path(out_dir, "fticr_tw2_pairwise.csv"))


# Dispersion test (betadisper / PERMDISP) ----
# Tests whether the spread (distance to group centroid) differs among groups.
# Reported with the Wd* tests: Wd* tests differences in location (centroids);
# betadisper tests differences in spread.
#   - between inundation histories (all samples)
#   - among DOM treatments within previously dry and within wet sediments
dispersion <- tibble()
bd <- betadisper(dist_bc, sample_info$Moisture)
set.seed(123)
bd_test <- permutest(bd, permutations = 9999)
dispersion <- bind_rows(dispersion, tibble(
  Comparison = "Inundation history (all samples)", Group = names(bd$group.distances),
  Mean_distance_to_centroid = as.numeric(bd$group.distances),
  F = bd_test$tab$F[1], p = bd_test$tab$`Pr(>F)`[1]))
for (mo in c("Dry", "Wet")) {
  samples_subset <- sample_info %>% filter(Moisture == mo)
  dist_subset <- as.dist(bc_mat[samples_subset$Sample_Name, samples_subset$Sample_Name])
  bd <- betadisper(dist_subset, droplevels(samples_subset$Treatment))
  set.seed(123)
  bd_test <- permutest(bd, permutations = 9999)
  dispersion <- bind_rows(dispersion, tibble(
    Comparison = paste("DOM treatment,", ifelse(mo == "Dry", "previously dry", "wet")),
    Group = names(bd$group.distances),
    Mean_distance_to_centroid = as.numeric(bd$group.distances),
    F = bd_test$tab$F[1], p = bd_test$tab$`Pr(>F)`[1]))
}
print(dispersion)
write_csv(dispersion, file.path(out_dir, "fticr_betadisper.csv"))


# Mean dissimilarity to previously dry samples ----
# Used in the Results to ask whether a DOM treatment shifted wet-sediment DOM
# toward the composition of previously dry sediments. For each group, the mean
# Bray-Curtis dissimilarity between its samples and all previously dry samples
# (for the previously dry group itself: the mean among previously dry samples).
distance_to_dry <- sample_info %>%
  group_by(Moisture, Treatment) %>%
  summarize(Samples = list(Sample_Name), .groups = "drop") %>%
  bind_rows(tibble(Moisture = factor("Dry", levels = levels(sample_info$Moisture)),
                   Treatment = NA, Samples = list(dry_samples))) %>%
  rowwise() %>%
  mutate(Mean_BC_to_dry = {
    d <- bc_mat[Samples, dry_samples, drop = FALSE]
    mean(d[outer(Samples, dry_samples, "!=")])
  }) %>%
  ungroup() %>%
  mutate(Group = ifelse(is.na(Treatment), "All previously dry", paste(Moisture, Treatment))) %>%
  select(Group, Moisture, Treatment, Mean_BC_to_dry)
print(distance_to_dry)
write_csv(distance_to_dry, file.path(out_dir, "fticr_mean_distance_to_dry.csv"))


# Create Plot Labels ----

## Wd* results for the treatment tests (shown on the NMDS)
wds_label <- wds_tests %>%
  filter(Grouping == "Treatment") %>%
  mutate(p_text = case_when(p < 0.001 ~ "< 0.001", TRUE ~ paste0("= ", sprintf("%.3f", p))))

## Create text labels for plot
plot_labels <- tibble(
  x = min(site_scores$NMDS1) + 0.03 * diff(range(site_scores$NMDS1)),
  y = max(site_scores$NMDS2) - c(0.06, 0.12, 0.18) * diff(range(site_scores$NMDS2)),
  label = c(
    sprintf("NMDS stress = %.2f", nmds$stress),
    sprintf("Previously dry, Wd* test: p %s", wds_label$p_text[1]),
    sprintf("Wet, Wd* test: p %s", wds_label$p_text[2])
  )
)


# Create NMDS Plot ----

nmds_plot <- ggplot(site_scores %>% mutate(Treatment = relabel_treatment(Treatment)),
                    aes(x = NMDS1, y = NMDS2, 
                        color = Treatment, shape = Moisture)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text(data = plot_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.4, color = "black") +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  scale_shape_manual(values = moisture_shapes, drop = FALSE) +
  labs(x = "NMDS1", y = "NMDS2", color = "DOM treatment", shape = "Inundation history") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right", aspect.ratio = 1)

print(nmds_plot)

# Save plot
# Saved below as panel (c) of the combined Figure 3 (see "Figure 3")


# Intensity-weighted AI_mod and NOSC per sample ----
# For each sample, AI_mod and NOSC of the formulas are averaged using their
# relative intensities as weights (sum of relative intensity x value).
# AI_Mod and NOSC come from the Mol file (rows match the Data file) and follow
# the standard formulas:
#   AI_mod = (1 + C - 0.5*O - S - 0.5*(N + P + H)) / (C - 0.5*O - S - N - P),
#            set to 0 when the numerator or denominator is <= 0 (Koch & Dittmar 2006, 2016)
#   NOSC   = 4 - (4*C + H - 3*N - 2*O + 5*P - 2*S) / C (Riedel et al. 2012)
# mol was read and filtered together with fticr_raw (see "Filter Formulas").
stopifnot(nrow(mol) == ncol(rel_matrix))

chemistry <- tibble(
  Sample_Name = rownames(rel_matrix),
  n_formulas  = rowSums(int_matrix > 0),
  AI_mod_wmean = as.numeric(rel_matrix %*% mol$AI_Mod),
  NOSC_wmean   = as.numeric(rel_matrix %*% mol$NOSC)) %>%
  left_join(sample_info %>% select(Sample_Name, Moisture, Treatment, Incubation_Day),
            by = "Sample_Name")
write_csv(chemistry, file.path(out_dir, "fticr_chemistry_per_sample.csv"))

chemistry_table <- chemistry %>%
  group_by(Moisture, Treatment) %>%
  summarize(n = n(),
            across(c(n_formulas, AI_mod_wmean, NOSC_wmean),
                   list(median = median, q25 = ~quantile(., 0.25), q75 = ~quantile(., 0.75))),
            .groups = "drop")
print(chemistry_table)
write_csv(chemistry_table, file.path(out_dir, "fticr_chemistry_by_group.csv"))


# Figure 3: (a) AI_mod, (b) NOSC, (c) NMDS ----
# (a, b): per-sample intensity-weighted AI_mod and NOSC by DOM treatment, within each
#   inundation history. Same style as Fig_Total_O2_by_Treatment: large symbol =
#   median, bar = interquartile range, small symbols = single samples.
#   Kruskal-Wallis among the three treatments in each panel; letters only where
#   it is significant (Dunn's test, no p adjustment; treatments sharing a
#   letter do not differ). Tests are saved to fticr_chemistry_tests.csv.
# (c): the NMDS above.
library(patchwork)   # combine panels
history_titles <- c(Dry = "Previously dry", Wet = "Wet")

chem_long <- chemistry %>%
  mutate(Treatment = relabel_treatment(Treatment)) %>%
  select(Sample_Name, Moisture, Treatment, AI_mod_wmean, NOSC_wmean) %>%
  pivot_longer(c(AI_mod_wmean, NOSC_wmean), names_to = "Response", values_to = "Value")

# Kruskal-Wallis, then Dunn for every pair (interpreted only if Kruskal-Wallis p < 0.05)
chem_tests <- tibble()
for (resp in unique(chem_long$Response)) for (mo in levels(chem_long$Moisture)) {
  d  <- chem_long %>% filter(Response == resp, Moisture == mo)
  kw <- kruskal.test(Value ~ Treatment, data = d)
  dn <- rstatix::dunn_test(d, Value ~ Treatment, p.adjust.method = "none")
  chem_tests <- bind_rows(chem_tests, tibble(
    Response = resp, Moisture = mo, KW_H = unname(kw$statistic), KW_p = kw$p.value,
    Treatment_1 = dn$group1, Treatment_2 = dn$group2, Dunn_Z = dn$statistic, Dunn_p = dn$p,
    significant = kw$p.value < 0.05 & dn$p < 0.05))
}
write_csv(chem_tests, file.path(out_dir, "fticr_chemistry_tests.csv"))

# Letters (compact letter display): where Kruskal-Wallis is significant,
# treatments are ordered by median and each letter marks a run of adjacent
# treatments in which no pair differs (Dunn p >= 0.05). Treatments that share a
# letter do not differ. Example: if only the lowest and highest treatments
# differ, they get "a" and "b" and the middle treatment gets "ab".
chem_letters <- function(resp) {
  out <- tibble()
  for (mo in levels(chem_long$Moisture)) {
    t <- chem_tests %>% filter(Response == resp, Moisture == mo)
    if (t$KW_p[1] >= 0.05) next
    meds <- chem_long %>% filter(Response == resp, Moisture == mo) %>%
      group_by(Treatment) %>% summarize(m = median(Value), .groups = "drop") %>% arrange(m)
    trts <- as.character(meds$Treatment)
    differs <- function(x, y) any(t$significant & ((t$Treatment_1 == x & t$Treatment_2 == y) |
                                                   (t$Treatment_1 == y & t$Treatment_2 == x)))
    # runs of adjacent treatments with no significant pair
    runs <- list()
    for (i in seq_along(trts)) {
      j <- i
      while (j < length(trts) &&
             !any(sapply(i:j, function(k) differs(trts[k], trts[j + 1])))) j <- j + 1
      runs[[length(runs) + 1]] <- i:j
    }
    # drop runs contained in an earlier, longer run
    runs <- runs[!sapply(seq_along(runs), function(r)
      any(sapply(seq_along(runs)[-r], function(q) all(runs[[r]] %in% runs[[q]]) &&
                   length(runs[[q]]) > length(runs[[r]]))))]
    runs <- unique(runs)
    lets <- sapply(seq_along(trts), function(k)
      paste(letters[which(sapply(runs, function(r) k %in% r))], collapse = ""))
    out <- bind_rows(out, tibble(Moisture = mo,
                                 Treatment = factor(trts, levels = levels(chem_long$Treatment)),
                                 Letter = lets))
  }
  out
}

chem_panel <- function(resp, ylab) {
  d <- chem_long %>% filter(Response == resp)
  g <- d %>% group_by(Moisture, Treatment) %>%
    summarize(median = median(Value), q25 = quantile(Value, 0.25), q75 = quantile(Value, 0.75),
              top = max(Value), .groups = "drop")
  span <- diff(range(d$Value))
  kw_lab <- chem_tests %>% filter(Response == resp) %>% distinct(Moisture, KW_p) %>%
    mutate(Moisture = factor(Moisture, levels = levels(d$Moisture)),
           label = sprintf("Kruskal-Wallis p = %.3f", KW_p), y = max(d$Value) + 0.16 * span)
  lets <- chem_letters(resp)
  if (nrow(lets) > 0) lets <- lets %>% mutate(Moisture = factor(Moisture, levels = levels(d$Moisture))) %>%
    left_join(g %>% select(Moisture, Treatment, top), by = c("Moisture", "Treatment")) %>%
    mutate(y = top + 0.07 * span)
  ggplot(g, aes(Treatment, median, color = Treatment, shape = Moisture)) +
    geom_jitter(data = d, aes(y = Value), width = 0.12, height = 0, size = 1.8, alpha = 0.45) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, linewidth = 0.8) +
    geom_point(size = 3.8) +
    geom_text(data = kw_lab, aes(x = 0.55, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3.6) +
    { if (nrow(lets) > 0) geom_text(data = lets, aes(x = Treatment, y = y, label = Letter),
                                    inherit.aes = FALSE, size = 4.5, fontface = "bold") } +
    facet_wrap(~ Moisture, labeller = labeller(Moisture = history_titles)) +
    scale_color_manual(values = treatment_colors, guide = "none") +
    scale_shape_manual(values = moisture_shapes, guide = "none") +
    scale_x_discrete(labels = c("Control", "Unburned\nDOM", "Burned\nDOM")) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
}

set.seed(1)   # same horizontal jitter of the small points every run
fig2 <- (chem_panel("AI_mod_wmean", expression(Intensity*"-"*weighted~AI[mod])) /
         chem_panel("NOSC_wmean", "Intensity-weighted NOSC")) |
        (nmds_plot + theme_bw(base_size = 12) + theme(aspect.ratio = 1, panel.grid.minor = element_blank()) +
           scale_shape_manual(values = moisture_shapes, labels = history_titles, drop = FALSE))
fig2 <- fig2 + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(file.path(fig_dir, "Fig3_FTICR.png"), fig2, width = 13, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig3_FTICR.pdf"), fig2, width = 13, height = 7)
cat(sprintf("\nFigure 3 saved. NMDS stress = %.3f (for the caption)\n", nmds$stress))