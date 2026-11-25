# =============================================================================
# 08_tables.R -- manuscript and SI tables
# -----------------------------------------------------------------------------
# Tables are numbered in the order they are first cited in the manuscript
# (renumbered 23 Nov 2026; old numbers in brackets):
# Table S1 [old S3]: sediment water extracts at collection and after
#            pre-incubation. Ions now per kg dry sediment, like WEOC.
# Table S2 [old S5]: treatment solutions: DOC, TDN, ions, and FT-ICR MS
#            intensity-weighted AI_mod and NOSC of the solutions as added
# Table S3 [old S4]: carbon added with the treatment solution vs
#            water-extractable OC, and the share of post-incubation DOC above
#            the controls
# Table S4 [old S2]: model fit for every vial at Km = 0.2 mg/L (least squares + MCMC)
# Table S5 [old S1]: median (IQR) of each response by sediment history x
#            treatment (rank-based tests -> medians); nitrate per kg dry
#            sediment; asterisks for every dry vs wet Mann-Whitney test.
#            A mean (SD) version is saved for reference.
# Data S1  : model fits for every vial at the seven Sihi-range Km values
#            discussed in the paper (least squares + MCMC; Supplementary Data CSV)
#
# Inputs : outputs/stats/microcosm_responses.csv, MannWhitney_*.csv (03_stats.R)
#          outputs/modeling/YEP_DataPrep_Log.csv, YEP_BestFit_AllKm.csv (01)
#          outputs/modeling/YEP_MCMC_Summary.csv (02)
#          data_files$water_mass, data_files$npoc_tn, data_files$ions,
#          data_files$fticr_data, data_files$fticr_mol
# Outputs: outputs/tables/*.csv
# =============================================================================

rm(list = ls(all = TRUE))
library(dplyr)
library(tidyr)
source("src/00_data_paths.R")
source("src/00_plot_style.R")   # treatment names used in the figures and tables

Km_base <- 0.2
out_dir <- file.path("outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
trt_levels <- c("Control", "Unburned + DOC", "High Burn + DOC")

# "mean (SD)" and "median (Q1-Q3)" text for tables
mean_sd <- function(x, digits) sprintf(paste0("%.", digits, "f (%.", digits, "f)"),
                                       mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
med_iqr <- function(x, digits) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), q[1], q[2], q[3])
}

# Treatment names shown in the tables. The calculations below use the treatment
# codes from the data ("Unburned + DOC", "High Burn + DOC"); this function swaps
# them for the manuscript names just before each table is saved.
table_names <- function(d) {
  codes <- as.character(d$Treatment)
  d$Treatment <- ifelse(codes %in% names(treatment_labels), treatment_labels[codes], codes)
  d
}

# Flagged cells look like "..._Final_Corrected"; use the corrected number,
# otherwise the cell itself (below-LOD flags become NA)
flag_value <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("Final_Corrected", x), sub(".*\\|([0-9.]+)_ppm_Final_Corrected.*", "\\1", x), x)
  suppressWarnings(as.numeric(x))
}

npoc     <- read_data(data_files$npoc_tn)
ions     <- read_data(data_files$ions)
mass_all <- read.csv(data_files$water_mass)
responses <- read.csv(file.path("outputs", "stats", "microcosm_responses.csv"))
responses$Treatment <- factor(responses$Treatment, levels = trt_levels)

# ---- Table S1: sediment extracts at collection vs after pre-incubation -------
# WEOC: mean (SD) of 3 replicate extracts (mg C per kg dry sediment, from the
# data package). Ions: one composite extract per date (n = 1, no SD), measured
# in mg/L of extract and converted here to mg per kg dry sediment with the
# extract water per dry sediment (water incl. pore water / dry mass, mean of the
# 3 replicate extracts; YEP*_SED-1..3 rows of the mass file). The 1:2 extracts
# used 10 g field-moist sediment, so water per kg dry sediment differs between
# dry (~2.1 L/kg) and wet (~3.4 L/kg) sediments and mg/L values are not
# comparable between histories.
extract_water <- mass_all %>%
  filter(grepl("^YEP[12][ABC]?_SED-[123]$", Sample_Name)) %>%
  mutate(Sample = sub("_SED-[123]$", "", Sample_Name),
         L_per_kg = as.numeric(Water_Mass_g) / as.numeric(Dry_Sediment_Mass_g)) %>%
  group_by(Sample) %>% summarize(Extract_L_per_kg = mean(L_per_kg))
extract_ions <- ions %>%
  filter(grepl("^YEP[12][ABC]?_SIN-1$", Sample_Name)) %>%
  transmute(Sample = sub("_SIN-1$", "", Sample_Name),
            NO3_N_mg_L = flag_value(`00618_NO3_mg_per_L_as_N`),
            SO4_mg_L   = flag_value(`00945_SO4_mg_per_L_as_SO4`))
extract_weoc <- npoc %>%
  filter(grepl("^YEP[12][ABC]?_SOC-[123]$", Sample_Name)) %>%
  mutate(Sample = sub("_SOC-[123]$", "", Sample_Name),
         WEOC = flag_value(Extractable_NPOC_mg_per_kg),
         WEOC_mg_L = flag_value(Extractable_NPOC_mg_per_L)) %>%
  group_by(Sample) %>%
  summarize(WEOC_txt = mean_sd(WEOC, 1),
            ratio_check = mean(WEOC / WEOC_mg_L))   # L/kg implied by the package's WEOC normalization

tableS1 <- extract_weoc %>% left_join(extract_water, by = "Sample") %>% left_join(extract_ions, by = "Sample") %>%
  mutate(Sediment = ifelse(substr(Sample, 4, 4) == "1", "Dry", "Wet"),
         Stage = ifelse(nchar(Sample) == 4, "Field (collection)",
                        paste0("Pre-incubated, incubation day ", match(substr(Sample, 5, 5), LETTERS))))
# The package normalized WEOC with the same water/dry-mass ratio (check, within 2%)
stopifnot(all(abs(tableS1$ratio_check / tableS1$Extract_L_per_kg - 1) < 0.02))
tableS1 <- tableS1 %>%
  arrange(Sediment, Sample) %>%
  transmute(Sediment, Stage, Sample,
            `Extract water (L kg-1 dry sediment)` = round(Extract_L_per_kg, 2),
            `WEOC (mg C kg-1), mean (SD), n = 3` = WEOC_txt,
            `NO3-N (mg N kg-1), n = 1` = signif(NO3_N_mg_L * Extract_L_per_kg, 3),
            `SO4 (mg kg-1), n = 1`     = signif(SO4_mg_L * Extract_L_per_kg, 3),
            `NO3-N in extract (mg N L-1)` = NO3_N_mg_L,
            `SO4 in extract (mg L-1)`     = SO4_mg_L)
write.csv(tableS1, file.path(out_dir, "TableS1_Sediment_Extracts.csv"), row.names = FALSE)

# ---- Table S2: treatment solutions ------------------------------------------
# DOC/TDN: solutions as added (OCN-SK*, -UK*, -HK*). Ions: water-only incubations
# (ION-WS*, -WU*, -WH*), i.e. the same solutions without sediment.
# FT-ICR MS: solutions as added (ICR-UKA-C, ICR-HKA-C; one per incubation day).
# Same formula filters as 06_fticr.R (formulas in > 2 of 4 methanol blanks
# removed; 200-900 m/z), relative intensity within each sample, and
# intensity-weighted AI_mod and NOSC. No FT-ICR sample of the control solution
# (synthetic river water, no DOM) was analyzed with the treatment solutions.
ion_cols <- c("00940_Cl_mg_per_L", "00945_SO4_mg_per_L_as_SO4", "00618_NO3_mg_per_L_as_N",
              "00915_Ca_mg_per_L", "00925_Mg_mg_per_L", "00935_K_mg_per_L", "00930_Na_mg_per_L")
recode_trt <- function(code) dplyr::recode(code, S = "Control", U = "Unburned + DOC", H = "High Burn + DOC")

solution_doc <- npoc %>%
  filter(grepl("^YEP_OCN-[SUH]K", Sample_Name)) %>%
  mutate(Treatment = recode_trt(sub("^YEP_OCN-([SUH]).*", "\\1", Sample_Name)),
         DOC = flag_value(Extractable_NPOC_mg_per_L), TDN = flag_value(Extractable_TN_mg_per_L)) %>%
  group_by(Treatment) %>% summarize(Solution_DOC_mg_L = mean(DOC), TDN_mg_L = mean(TDN, na.rm = TRUE))

solution_ions <- ions %>%
  filter(grepl("^YEP_ION-W[SUH]", Sample_Name)) %>%
  mutate(Treatment = recode_trt(sub("^YEP_ION-W([SUH]).*", "\\1", Sample_Name))) %>%
  mutate(across(all_of(ion_cols), flag_value)) %>%
  group_by(Treatment) %>% summarize(across(all_of(ion_cols), mean))

fticr <- read.csv(data_files$fticr_data, check.names = FALSE)
mol   <- read.csv(data_files$fticr_mol, check.names = FALSE)
stopifnot(nrow(mol) == nrow(fticr), all(abs(mol$Calibrated_Mass - fticr$Calibrated_Mass) < 1e-6))
blank_cols <- grep("MeOH_blank", names(fticr), value = TRUE)
stopifnot(length(blank_cols) == 4)
keep <- rowSums(fticr[blank_cols] > 0, na.rm = TRUE) <= length(blank_cols) / 2 &
        fticr$Calibrated_Mass >= 200 & fticr$Calibrated_Mass <= 900
solution_icr <- grep("^YEP_ICR-[UH]K[ABC]$", names(fticr), value = TRUE)
stopifnot(length(solution_icr) == 6)
icr_int <- as.matrix(fticr[keep, solution_icr]); icr_int[is.na(icr_int)] <- 0
icr_rel <- sweep(icr_int, 2, colSums(icr_int), "/")
solution_chem <- data.frame(
  Sample_Name  = solution_icr,
  Treatment    = recode_trt(sub("^YEP_ICR-([UH])K.*", "\\1", solution_icr)),
  n_formulas   = colSums(icr_int > 0),
  AI_mod_wmean = colSums(icr_rel * mol$AI_Mod[keep]),
  NOSC_wmean   = colSums(icr_rel * mol$NOSC[keep]))
write.csv(solution_chem, file.path(out_dir, "TableS2_support_Solution_FTICR_per_sample.csv"), row.names = FALSE)
solution_chem_sum <- solution_chem %>% group_by(Treatment) %>%
  summarize(`FT-ICR formulas, mean (SD), n = 3` = mean_sd(n_formulas, 0),
            `AI_mod (intensity-weighted), mean (SD)` = mean_sd(AI_mod_wmean, 3),
            `NOSC (intensity-weighted), mean (SD)`   = mean_sd(NOSC_wmean, 3))

tableS2 <- solution_doc %>% left_join(solution_ions, by = "Treatment") %>%
  mutate(across(where(is.numeric), ~ signif(., 3))) %>%
  left_join(solution_chem_sum, by = "Treatment") %>%
  mutate(Treatment = factor(Treatment, levels = trt_levels)) %>% arrange(Treatment) %>%
  rename(`DOC (mg C L-1)` = Solution_DOC_mg_L, `TDN (mg N L-1)` = TDN_mg_L,
         `Cl (mg L-1)` = `00940_Cl_mg_per_L`, `SO4 (mg L-1)` = `00945_SO4_mg_per_L_as_SO4`,
         `NO3-N (mg N L-1)` = `00618_NO3_mg_per_L_as_N`, `Ca (mg L-1)` = `00915_Ca_mg_per_L`,
         `Mg (mg L-1)` = `00925_Mg_mg_per_L`, `K (mg L-1)` = `00935_K_mg_per_L`, `Na (mg L-1)` = `00930_Na_mg_per_L`) %>%
  mutate(across(starts_with(c("FT-ICR", "AI_mod", "NOSC")), ~ ifelse(is.na(.), "n.a.", .))) %>%
  table_names()
write.csv(tableS2, file.path(out_dir, "TableS2_Treatment_Solutions.csv"), row.names = FALSE)

# ---- Table S3: carbon added vs water-extractable OC -------------------------
# Added solution volume per vial = total water - pore water already in the
# sediment (field-moist mass - dry mass). Added C = solution DOC x volume / dry mass.
# Post-incubation DOC: median of the 5 vials. "Above control" = treatment median
# - control median of the same history; its share of the treatment median is
# the fraction of post-incubation DOC (and of the carbon analyzed by FT-ICR MS,
# which was normalized to 1.5 mg C/L) that can be attributed to the addition.
mass <- mass_all %>%
  filter(grepl("^YEP[12][ABC]_INC-[HSU]", Sample_Name)) %>%
  mutate(Moisture  = ifelse(substr(Sample_Name, 4, 4) == "1", "Dry", "Wet"),
         Treatment = recode_trt(sub(".*_INC-([HSU]).*", "\\1", Sample_Name)),
         Pore_water_mL      = as.numeric(Field_Moist_Sediment_Mass_g) - as.numeric(Dry_Sediment_Mass_g),
         Added_solution_mL  = as.numeric(Water_Mass_g) - Pore_water_mL)

weoc <- npoc %>%
  filter(grepl("^YEP[12][ABC]?_SOC-[123]$", Sample_Name)) %>%
  mutate(Moisture = ifelse(substr(Sample_Name, 4, 4) == "1", "Dry", "Wet"),
         Stage    = ifelse(grepl("^YEP[12]_", Sample_Name), "Field", "Pre-incubation"),
         WEOC     = flag_value(Extractable_NPOC_mg_per_kg)) %>%
  group_by(Moisture, Stage) %>% summarize(WEOC_mg_kg = mean(WEOC), .groups = "drop") %>%
  pivot_wider(names_from = Stage, values_from = WEOC_mg_kg, names_prefix = "WEOC_mg_kg_")

post_doc <- responses %>% group_by(Moisture, Treatment) %>%
  summarize(Post_DOC_median = median(DOC_mg_kg), .groups = "drop") %>%
  mutate(Treatment = as.character(Treatment)) %>%
  group_by(Moisture) %>%
  mutate(Post_DOC_above_control = Post_DOC_median - Post_DOC_median[Treatment == "Control"],
         Share_above_control_pct = 100 * Post_DOC_above_control / Post_DOC_median) %>% ungroup()

tableS3 <- mass %>%
  left_join(solution_doc %>% select(Treatment, Solution_DOC_mg_L), by = "Treatment") %>%
  mutate(Added_C_mg_kg = Solution_DOC_mg_L * (Added_solution_mL / 1000) / (as.numeric(Dry_Sediment_Mass_g) / 1000)) %>%
  group_by(Moisture, Treatment) %>%
  summarize(Solution_DOC_mg_L = mean(Solution_DOC_mg_L), Added_solution_mL = mean(Added_solution_mL),
            Added_C_mg_kg = mean(Added_C_mg_kg), .groups = "drop") %>%
  left_join(weoc, by = "Moisture") %>%
  left_join(post_doc, by = c("Moisture", "Treatment")) %>%
  mutate(Added_C_pct_of_field_WEOC          = 100 * Added_C_mg_kg / `WEOC_mg_kg_Field`,
         Added_C_pct_of_preincubation_WEOC  = 100 * Added_C_mg_kg / `WEOC_mg_kg_Pre-incubation`,
         Share_above_control_pct = ifelse(Treatment == "Control", NA, Share_above_control_pct),
         Post_DOC_above_control  = ifelse(Treatment == "Control", NA, Post_DOC_above_control)) %>%
  mutate(Treatment = factor(Treatment, levels = trt_levels)) %>%
  arrange(Moisture, Treatment) %>%
  mutate(across(where(is.numeric), ~ round(., 1))) %>%
  transmute(Moisture, Treatment,
            `Solution DOC (mg C L-1)` = Solution_DOC_mg_L,
            `Solution added (mL per vial)` = Added_solution_mL,
            `C added (mg C kg-1)` = Added_C_mg_kg,
            `WEOC at collection (mg C kg-1)` = `WEOC_mg_kg_Field`,
            `WEOC after pre-incubation (mg C kg-1)` = `WEOC_mg_kg_Pre-incubation`,
            `C added (% of WEOC at collection)` = Added_C_pct_of_field_WEOC,
            `C added (% of WEOC after pre-incubation)` = Added_C_pct_of_preincubation_WEOC,
            `Post-incubation DOC, median (mg C kg-1)` = Post_DOC_median,
            `Post-incubation DOC above control (mg C kg-1)` = Post_DOC_above_control,
            `DOC above control (% of post-incubation DOC)` = Share_above_control_pct) %>%
  table_names()
write.csv(tableS3, file.path(out_dir, "TableS3_Carbon_Added_vs_WEOC.csv"), row.names = FALSE, na = "")

# ---- Table S4: per-vial model fits at Km = 0.2 ----------------------------------
prep <- read.csv(file.path("outputs", "modeling", "YEP_DataPrep_Log.csv"))
best <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_AllKm.csv")) %>% filter(Km_mg_L == Km_base)
mcmc <- read.csv(file.path("outputs", "modeling", "YEP_MCMC_Summary.csv")) %>% filter(Km_mg_L == Km_base)

tableS4 <- prep %>%
  left_join(best %>% select(Sample_Name, n_points, RMSE, R2, kL_L_kg_h, Vmax_mg_kg_h), by = "Sample_Name") %>%
  left_join(mcmc %>% select(Sample_Name, kL_kg_median, kL_kg_lower95, kL_kg_upper95,
                            Vmax_kg_median, Vmax_kg_lower95, Vmax_kg_upper95,
                            r_kL_Vmax, Rhat_kL, Rhat_Vmax, ESS_bulk_min, ESS_tail_min, converged),
            by = "Sample_Name") %>%
  arrange(Moisture, factor(Treatment, levels = c("Control", "Unburned + DOC", "High Burn + DOC")), Sample_Name) %>%
  transmute(Sample = Sample_Name, Sediment = Moisture, Treatment, Day,
            `Points fitted` = n_fit,
            `Fit duration (min)` = duration_fit_min,
            `DO range fitted (mg/L)` = sprintf("%.2f-%.2f", DO_last_obs, DO_first_obs),
            `Water (L per kg sediment)` = round(L_water_per_kg, 2),
            `RMSE (mg/L)` = signif(RMSE, 2),
            `kL, least squares (L kg-1 h-1)` = round(kL_L_kg_h, 2),          # per kg, same units
            `Vmax, least squares (mg O2 kg-1 h-1)` = round(Vmax_mg_kg_h, 1), # and digits as the MCMC columns
            `kL (L kg-1 h-1), posterior median [95% CrI]` =
              sprintf("%.2f [%.2f-%.2f]", kL_kg_median, kL_kg_lower95, kL_kg_upper95),
            `Vmax (mg O2 kg-1 h-1), posterior median [95% CrI]` =
              sprintf("%.1f [%.1f-%.1f]", Vmax_kg_median, Vmax_kg_lower95, Vmax_kg_upper95),
            `r(kL, Vmax)` = round(r_kL_Vmax, 2),
            `R-hat kL` = round(Rhat_kL, 3), `R-hat Vmax` = round(Rhat_Vmax, 3),
            `Min bulk ESS` = round(ESS_bulk_min), `Min tail ESS` = round(ESS_tail_min),
            Converged = ifelse(converged, "yes", "no")) %>%
  table_names()
write.csv(tableS4, file.path(out_dir, "TableS4_Model_Fits_Km0.2.csv"), row.names = FALSE)

# ---- Data S1: per-vial model fits at every Km (Supplementary Data) -----------
# Same fits as Table S4, for every Km discussed in the paper. Plain numbers (no
# formatting), one row per vial x Km. 01_model_DO.R also fits Km = 0.05, 0.1,
# 3.2, 5.28 and 10 mg/L for the response to reviewers; those are left out here so
# Data S1 holds only the seven Km values in the Sihi et al. (2020) range
# (0.13-1.22 mg/L), which are also the Km values with MCMC.
Km_paper <- c(0.13, 0.2, 0.3, 0.51, 0.75, 1, 1.22)   # mg/L, as in 02_model_mcmc.R
best_all <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_AllKm.csv")) %>%
  filter(round(Km_mg_L, 2) %in% Km_paper)
mcmc_all <- read.csv(file.path("outputs", "modeling", "YEP_MCMC_Summary.csv"))

dataS1 <- best_all %>%
  select(Sample_Name, Moisture, Treatment, Day, Km_mg_L, n_points, RMSE, R2, kL_L_kg_h, Vmax_mg_kg_h) %>%
  left_join(mcmc_all %>% select(Sample_Name, Km_mg_L, kL_kg_median, kL_kg_lower95, kL_kg_upper95,
                                Vmax_kg_median, Vmax_kg_lower95, Vmax_kg_upper95,
                                r_kL_Vmax, Rhat_kL, Rhat_Vmax, ESS_bulk_min, ESS_tail_min, converged),
            by = c("Sample_Name", "Km_mg_L")) %>%
  arrange(Km_mg_L, Moisture, factor(Treatment, levels = c("Control", "Unburned + DOC", "High Burn + DOC")),
          Sample_Name) %>%
  transmute(Sample = Sample_Name, Sediment = Moisture, Treatment, Day,
            `Km (mg/L)` = Km_mg_L,
            `Points fitted` = n_points,
            `RMSE (mg/L)` = signif(RMSE, 3),
            R2 = round(R2, 4),
            `kL, least squares (L kg-1 h-1)` = signif(kL_L_kg_h, 4),
            `Vmax, least squares (mg O2 kg-1 h-1)` = signif(Vmax_mg_kg_h, 4),
            `kL, posterior median (L kg-1 h-1)` = signif(kL_kg_median, 4),
            `kL, 95% CrI lower` = signif(kL_kg_lower95, 4),
            `kL, 95% CrI upper` = signif(kL_kg_upper95, 4),
            `Vmax, posterior median (mg O2 kg-1 h-1)` = signif(Vmax_kg_median, 4),
            `Vmax, 95% CrI lower` = signif(Vmax_kg_lower95, 4),
            `Vmax, 95% CrI upper` = signif(Vmax_kg_upper95, 4),
            `r(kL, Vmax)` = round(r_kL_Vmax, 2),
            `R-hat kL` = round(Rhat_kL, 3), `R-hat Vmax` = round(Rhat_Vmax, 3),
            `Min bulk ESS` = round(ESS_bulk_min), `Min tail ESS` = round(ESS_tail_min),
            Converged = ifelse(is.na(converged), NA, ifelse(converged, "yes", "no"))) %>%
  table_names()
# 30 vials x 7 Km, and every row has an MCMC result
stopifnot(nrow(dataS1) == 30 * length(Km_paper),
          nrow(dataS1) == nrow(mcmc_all),
          all(!is.na(dataS1$`kL, posterior median (L kg-1 h-1)`)))
write.csv(dataS1, file.path(out_dir, "DataS1_Model_Fits_AllKm.csv"), row.names = FALSE, na = "")

# ---- Table S5: group summaries -----------------------------------------------
# Median (Q1-Q3) because every test is rank-based (03_stats.R). Nitrate per kg
# dry sediment (normalized like DOC). Asterisk (*) on the wet "All treatments"
# value where wet differs from previously dry, p < 0.05 (Mann-Whitney,
# treatments pooled, n = 15 vs 15; Vmax and kL at Km = 0.2).
s5_cols <- c(O2_rate_mg_kg_h = "O2 consumption (mg O2 kg-1 h-1)",
             Vmax_mg_kg_h    = "Vmax (mg O2 kg-1 h-1)",
             kL_L_kg_h       = "kL (L kg-1 h-1)",
             DOC_mg_kg       = "DOC (mg C kg-1)",
             NO3_mg_N_kg     = "NO3-N (mg N kg-1)")
s5_digits <- c(O2_rate_mg_kg_h = 1, Vmax_mg_kg_h = 1, kL_L_kg_h = 2, DOC_mg_kg = 1, NO3_mg_N_kg = 1)

summarize_groups <- function(d, f) {
  out <- d %>% summarize(n = n(), .groups = "drop")
  for (col in names(s5_cols)) {
    out[[s5_cols[[col]]]] <- (d %>% summarize(v = f(.data[[col]], s5_digits[[col]]), .groups = "drop"))$v
  }
  out
}
build_s5 <- function(f) {
  by_trt <- summarize_groups(responses %>% group_by(Moisture, Treatment), f) %>%
    mutate(Treatment = as.character(Treatment))
  pooled <- summarize_groups(responses %>% group_by(Moisture), f) %>% mutate(Treatment = "All treatments")
  bind_rows(by_trt, pooled) %>% select(Moisture, Treatment, everything())
}

stats_dir   <- file.path("outputs", "stats")
dry_wet_sig <- bind_rows(
  read.csv(file.path(stats_dir, "MannWhitney_DryWet_TotalO2.csv")) %>% transmute(Response, p),
  read.csv(file.path(stats_dir, "MannWhitney_DryWet_Chemistry.csv")) %>% transmute(Response, p),
  read.csv(file.path(stats_dir, "MannWhitney_DryWet_byKm.csv")) %>% filter(Km_mg_L == Km_base) %>%
    transmute(Response = Parameter, p))
stopifnot(all(names(s5_cols) %in% dry_wet_sig$Response))
add_stars <- function(tab) {
  for (resp in names(s5_cols)) {
    if (dry_wet_sig$p[dry_wet_sig$Response == resp] < 0.05) {
      rows <- tab$Moisture == "Wet" & tab$Treatment == "All treatments"
      tab[[s5_cols[[resp]]]][rows] <- paste0(tab[[s5_cols[[resp]]]][rows], "*")
    }
  }
  tab
}
tableS5      <- build_s5(med_iqr) %>% add_stars() %>% table_names()
tableS5_mean <- build_s5(mean_sd) %>% add_stars() %>% table_names()
write.csv(tableS5,      file.path(out_dir, "TableS5_Group_Medians.csv"), row.names = FALSE)
write.csv(tableS5_mean, file.path(out_dir, "TableS5_Group_Means_reference.csv"), row.names = FALSE)

options(width = 200)
print(tableS1); print(tableS2); print(tableS3); print(tableS5)
