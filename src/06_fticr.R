# =============================================================================
# 06_fticr.R -- FTICR-MS: NMDS, PERMANOVA, dispersion, aromaticity (Fig. 2)
# -----------------------------------------------------------------------------
# Post-incubation aqueous samples (SIR), presence/absence of molecular formulas.
#   - NMDS on Jaccard distance (Fig. 2)
#   - PERMANOVA across all samples (inundation history, treatment, day, and
#     the history x treatment interaction)
#   - PERMANOVA (Treatment + Incubation_Day) within dry and within wet
#   - betadisper: are group dispersions equal? (PERMANOVA assumes they are)
#   - aromaticity and oxidation state per sample (table)
#   - Figure 2: (a) mean AI_mod and (b) mean NOSC by DOM treatment within each
#     history (Kruskal-Wallis, Dunn letters), (c) NMDS
# Note: YEP2A_SIR-S2 and YEP2B_SIR-S3 were reacquired at 0.2 s ion
#   accumulation (low signal at 0.075 s); the processed data package contains
#   only the 0.2 s spectra for these two samples, so no special handling here.
#
# Inputs : data_files$fticr_data, data_files$fticr_mol, data_files$sample_metadata
# Outputs: outputs/figures/Fig2_FTICR.pdf/.png  (replaces Fig2_FTICR_NMDS)
#          outputs/fticr/*.csv  (incl. fticr_chemistry_tests.csv)
# =============================================================================

library(tidyverse)
library(vegan)

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


# Create Presence/Absence Matrix ----

## Extract sample columns (exclude mass column)
sample_cols <- setdiff(names(fticr_raw), "Calibrated_Mass")

## Convert to presence/absence (1 if detected, 0 if not)
pa_matrix <- fticr_raw %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), ~as.numeric(. > 0 & !is.na(.)))) %>%
  as.matrix() %>%
  t()  # Transpose: rows = samples, columns = molecules

rownames(pa_matrix) <- sample_cols
colnames(pa_matrix) <- paste0("mol_", seq_len(nrow(fticr_raw)))


# Filter Samples ----

## Keep only target samples with complete metadata
sample_info <- tibble(Sample_Name = rownames(pa_matrix)) %>%
  left_join(metadata, by = "Sample_Name") %>%
  filter(
    Sample_Material == target_sample_type,
    !is.na(Treatment), !is.na(Moisture), !is.na(Incubation_Day)
  )

## Filter presence/absence matrix to match
pa_matrix <- pa_matrix[sample_info$Sample_Name, , drop = FALSE]

## Remove molecules never detected & samples with no detections
pa_matrix <- pa_matrix[, colSums(pa_matrix) > 0, drop = FALSE]
pa_matrix <- pa_matrix[rowSums(pa_matrix) > 0, , drop = FALSE]

## Update sample info to match filtered matrix
sample_info <- sample_info %>%
  filter(Sample_Name %in% rownames(pa_matrix)) %>%
  arrange(match(Sample_Name, rownames(pa_matrix)))

stopifnot(identical(sample_info$Sample_Name, rownames(pa_matrix)))

## Number of samples per group (Fig. 2 caption)
sample_counts <- sample_info %>% count(Moisture, Treatment, name = "n_samples")
print(sample_counts)
write_csv(sample_counts, file.path(out_dir, "fticr_sample_counts.csv"))


# Run NMDS Ordination ----

## Calculate Jaccard distance (presence/absence similarity)
dist_jaccard <- vegdist(pa_matrix, method = "jaccard", binary = TRUE)

## Run NMDS to create 2D representation
nmds <- metaMDS(pa_matrix, distance = "jaccard", binary = TRUE, k = 2, 
                trymax = 200, autotransform = FALSE, trace = FALSE)

## Extract sample coordinates and add metadata
site_scores <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_info, by = "Sample_Name")

write_csv(site_scores, file.path(out_dir, "fticr_pa_nmds_site_scores.csv"))


# Statistical Testing: PERMANOVA ----
# Test if treatments differ significantly, separately for dry and wet samples

run_permanova <- function(moisture_level) {
  # Subset data for this moisture level
  samples_subset <- sample_info %>% filter(Moisture == moisture_level)
  dist_subset <- as.dist(as.matrix(dist_jaccard)[samples_subset$Sample_Name, 
                                                 samples_subset$Sample_Name])
  
  # Run PERMANOVA
  set.seed(123)
  adonis2(dist_subset ~ Treatment + Incubation_Day, 
          data = samples_subset, permutations = 9999, by = "margin")
}

## Run for dry and wet separately
dry_permanova <- run_permanova("Dry")
wet_permanova <- run_permanova("Wet")

## Save full results
write_csv(as.data.frame(dry_permanova) %>% rownames_to_column("Term"),
          file.path(out_dir, "fticr_pa_permanova_dry_only.csv"))
write_csv(as.data.frame(wet_permanova) %>% rownames_to_column("Term"),
          file.path(out_dir, "fticr_pa_permanova_wet_only.csv"))


# Statistical Testing: PERMANOVA across all samples ----
# Does composition differ between dry and wet sediments, and among treatments,
# when all samples are analyzed together?
#   Model 1 (main effects): each of inundation history, treatment and day is
#     tested while accounting for the other two (marginal tests).
#   Model 2: adds the history x treatment interaction, i.e. whether the
#     treatment effect differs between dry and wet sediments. Only the
#     interaction row is kept from this model.
stopifnot(identical(labels(dist_jaccard), sample_info$Sample_Name))   # same sample order

set.seed(123)
overall_main <- adonis2(dist_jaccard ~ Moisture + Treatment + Incubation_Day,
                        data = sample_info, permutations = 9999, by = "margin")
set.seed(123)
overall_interaction <- adonis2(dist_jaccard ~ Moisture * Treatment + Incubation_Day,
                               data = sample_info, permutations = 9999, by = "margin")

overall_permanova <- bind_rows(
  as.data.frame(overall_main) %>% rownames_to_column("Term") %>%
    mutate(Model = "Moisture + Treatment + Incubation_Day"),
  as.data.frame(overall_interaction) %>% rownames_to_column("Term") %>%
    filter(Term == "Moisture:Treatment") %>%
    mutate(Model = "Moisture * Treatment + Incubation_Day")) %>%
  mutate(n_samples = nrow(sample_info)) %>%
  relocate(Model, n_samples)
print(overall_permanova)
write_csv(overall_permanova, file.path(out_dir, "fticr_pa_permanova_overall.csv"))


# Dispersion test (betadisper) ----
# PERMANOVA can report a group difference when groups differ only in spread.
# betadisper tests whether the spread (distance to group centroid) differs
# among treatments, within dry and within wet.
dispersion <- tibble()
for (mo in c("Dry", "Wet")) {
  samples_subset <- sample_info %>% filter(Moisture == mo)
  dist_subset <- as.dist(as.matrix(dist_jaccard)[samples_subset$Sample_Name,
                                                 samples_subset$Sample_Name])
  bd <- betadisper(dist_subset, samples_subset$Treatment)
  set.seed(123)
  bd_test <- permutest(bd, permutations = 9999)
  dispersion <- bind_rows(dispersion, tibble(
    Moisture = mo, Treatment = names(bd$group.distances),
    Mean_distance_to_centroid = as.numeric(bd$group.distances),
    F = bd_test$tab$F[1], p = bd_test$tab$`Pr(>F)`[1]))
}
print(dispersion)
write_csv(dispersion, file.path(out_dir, "fticr_pa_betadisper.csv"))


# Create Plot Labels ----

## Extract treatment effect statistics
extract_treatment_stats <- function(permanova_result, moisture_name) {
  as.data.frame(permanova_result) %>%
    rownames_to_column("Term") %>%
    filter(Term == "Treatment") %>%
    mutate(Moisture = moisture_name)
}

permanova_stats <- bind_rows(
  extract_treatment_stats(dry_permanova, "Dry"),
  extract_treatment_stats(wet_permanova, "Wet")
) %>%
  mutate(
    p_text = case_when(
      `Pr(>F)` < 0.001 ~ "< 0.001",
      TRUE ~ paste0("= ", sprintf("%.3f", `Pr(>F)`))
    )
  )

## Create text labels for plot
plot_labels <- tibble(
  x = min(site_scores$NMDS1) + 0.03 * diff(range(site_scores$NMDS1)),
  y = max(site_scores$NMDS2) - c(0.06, 0.12, 0.18) * diff(range(site_scores$NMDS2)),
  label = c(
    sprintf("NMDS stress = %.2f", nmds$stress),
    sprintf("Dry PERMANOVA: R\u00b2 = %.2f, p %s", 
            permanova_stats$R2[1], permanova_stats$p_text[1]),
    sprintf("Wet PERMANOVA: R\u00b2 = %.2f, p %s", 
            permanova_stats$R2[2], permanova_stats$p_text[2])
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
# Saved below as panel (c) of the combined Figure 2 (see "Figure 2")


# Aromaticity and oxidation state per sample ----
# Uses the formulas present in each sample (presence/absence, not weighted by
# intensity). AI_Mod and NOSC come from the Mol file (rows match the Data file)
# and follow the standard formulas:
#   AI_mod = (1 + C - 0.5*O - S - 0.5*(N + P + H)) / (C - 0.5*O - S - N - P),
#            set to 0 when the numerator or denominator is <= 0 (Koch & Dittmar 2006, 2016)
#   NOSC   = 4 - (4*C + H - 3*N - 2*O + 5*P - 2*S) / C (Riedel et al. 2012)
# Classes: condensed aromatic AI_mod >= 0.67; aromatic 0.5 < AI_mod < 0.67.
mol <- read_csv(mol_file, show_col_types = FALSE)
stopifnot(nrow(mol) == nrow(fticr_raw), all(abs(mol$Calibrated_Mass - fticr_raw$Calibrated_Mass) < 1e-6))

aromaticity <- tibble()
for (s in sample_info$Sample_Name) {
  present <- !is.na(fticr_raw[[s]]) & fticr_raw[[s]] > 0
  aromaticity <- bind_rows(aromaticity, tibble(
    Sample_Name = s,
    n_formulas = sum(present),
    AI_mod_mean = mean(mol$AI_Mod[present]),
    NOSC_mean = mean(mol$NOSC[present]),
    Pct_condensed_aromatic = 100 * mean(mol$AI_Mod[present] >= 0.67),
    Pct_aromatic = 100 * mean(mol$AI_Mod[present] > 0.5 & mol$AI_Mod[present] < 0.67)))
}
aromaticity <- aromaticity %>% left_join(sample_info %>% select(Sample_Name, Moisture, Treatment, Incubation_Day),
                                         by = "Sample_Name")
write_csv(aromaticity, file.path(out_dir, "fticr_aromaticity_per_sample.csv"))

aromaticity_table <- aromaticity %>%
  group_by(Moisture, Treatment) %>%
  summarize(n = n(),
            across(c(n_formulas, AI_mod_mean, NOSC_mean, Pct_condensed_aromatic, Pct_aromatic),
                   list(mean = mean, sd = sd)), .groups = "drop")
print(aromaticity_table)
write_csv(aromaticity_table, file.path(out_dir, "fticr_aromaticity_by_group.csv"))


# Figure 2: (a) AI_mod, (b) NOSC, (c) NMDS ----
# (a, b): per-sample mean AI_mod and NOSC by DOM treatment, within each
#   inundation history. Same style as Fig_Total_O2_by_Treatment: large symbol =
#   median, bar = interquartile range, small symbols = single samples.
#   Kruskal-Wallis among the three treatments in each panel; letters only where
#   it is significant (Dunn's test, no p adjustment; treatments sharing a
#   letter do not differ). Tests are saved to fticr_chemistry_tests.csv.
# (c): the NMDS above.
library(patchwork)   # combine panels
history_titles <- c(Dry = "Previously dry", Wet = "Wet")

chem_long <- aromaticity %>%
  mutate(Treatment = relabel_treatment(Treatment)) %>%
  select(Sample_Name, Moisture, Treatment, AI_mod_mean, NOSC_mean, Pct_aromatic) %>%
  pivot_longer(c(AI_mod_mean, NOSC_mean, Pct_aromatic), names_to = "Response", values_to = "Value")

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

# Letters for three groups: where Kruskal-Wallis is significant, a treatment
# that differs from both others gets "b" and the other two get "a"; any other
# pattern stops the script so the letters can be checked by hand.
chem_letters <- function(resp) {
  out <- tibble()
  for (mo in levels(chem_long$Moisture)) {
    t <- chem_tests %>% filter(Response == resp, Moisture == mo)
    if (t$KW_p[1] >= 0.05) next
    trts <- levels(chem_long$Treatment)
    n_sig <- sapply(trts, function(x) sum(t$significant & (t$Treatment_1 == x | t$Treatment_2 == x)))
    stopifnot(sum(n_sig == 2) == 1, sum(t$significant) == 2)
    out <- bind_rows(out, tibble(Moisture = mo, Treatment = factor(trts, levels = trts),
                                 Letter = ifelse(n_sig == 2, "b", "a")))
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
fig2 <- (chem_panel("AI_mod_mean", expression(Mean~AI[mod])) /
         chem_panel("NOSC_mean", "Mean NOSC")) |
        (nmds_plot + theme_bw(base_size = 12) + theme(aspect.ratio = 1, panel.grid.minor = element_blank()) +
           scale_shape_manual(values = moisture_shapes, labels = history_titles, drop = FALSE))
fig2 <- fig2 + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(file.path(fig_dir, "Fig2_FTICR.png"), fig2, width = 13, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig2_FTICR.pdf"), fig2, width = 13, height = 7)
cat(sprintf("\nFigure 2 saved. NMDS stress = %.3f (for the caption)\n", nmds$stress))