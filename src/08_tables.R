# =============================================================================
# 08_tables.R -- manuscript and SI tables
# -----------------------------------------------------------------------------
# Table S1 : mean (SD) of each response by sediment history x treatment
#            (replaces the old Table_S1.R; CO2 rates removed)
# Table S2 : model fit for every vial at Km = 0.2 mg/L (least squares + MCMC)
# Table S3 : sediment extracts at collection vs after pre-incubation
# Table S4 : carbon added with the treatment solution vs water-extractable OC
# Table S5 : chemistry of the treatment solutions
# Data S1  : model fits for every vial at every Km (Supplementary Data CSV;
#            least squares at all Km, MCMC at the Sihi-range Km)
# (numbered as in the manuscript SI; the code below builds them in the order
#  S1, S2, Data S1, S4, S5, S3)
#
# Inputs : outputs/stats/microcosm_responses.csv (03_stats.R)
#          outputs/modeling/YEP_DataPrep_Log.csv, YEP_BestFit_AllKm.csv (01)
#          outputs/modeling/YEP_MCMC_Summary.csv (02)
#          data_files$water_mass, data_files$npoc_tn, data_files$ions
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

# "mean (SD)" text for tables
mean_sd <- function(x, digits) sprintf(paste0("%.", digits, "f (%.", digits, "f)"),
                                       mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

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

# ---- Table S1: group means -------------------------------------------------------
responses <- read.csv(file.path("outputs", "stats", "microcosm_responses.csv"))
responses$Treatment <- factor(responses$Treatment, levels = c("Control", "Unburned + DOC", "High Burn + DOC"))

tableS1 <- responses %>%
  group_by(Moisture, Treatment) %>%
  summarize(n = n(),
            `O2 consumption (mg O2 kg-1 h-1)` = mean_sd(O2_rate_mg_kg_h, 1),
            `Vmax (mg O2 kg-1 h-1)`           = mean_sd(Vmax_mg_kg_h, 1),
            `kL (L kg-1 h-1)`                 = mean_sd(kL_L_kg_h, 2),
            `DOC (mg C kg-1)`                 = mean_sd(DOC_mg_kg, 1),
            `NO3-N (mg L-1)`                  = mean_sd(NO3_mg_N_L, 2),
            .groups = "drop") %>%
  arrange(Moisture, Treatment)

# Pooled dry and wet rows (all treatments)
tableS1_pooled <- responses %>%
  group_by(Moisture) %>%
  summarize(Treatment = "All treatments", n = n(),
            `O2 consumption (mg O2 kg-1 h-1)` = mean_sd(O2_rate_mg_kg_h, 1),
            `Vmax (mg O2 kg-1 h-1)`           = mean_sd(Vmax_mg_kg_h, 1),
            `kL (L kg-1 h-1)`                 = mean_sd(kL_L_kg_h, 2),
            `DOC (mg C kg-1)`                 = mean_sd(DOC_mg_kg, 1),
            `NO3-N (mg L-1)`                  = mean_sd(NO3_mg_N_L, 2))
tableS1 <- bind_rows(tableS1 %>% mutate(Treatment = as.character(Treatment)), tableS1_pooled)

# Mark (*) the wet "All treatments" value where wet differs from previously
# dry, p < 0.05 (03_stats.R: Mann-Whitney, treatments pooled, n = 15 vs 15;
# Vmax and kL at Km = 0.2).
stats_dir   <- file.path("outputs", "stats")
dry_wet_sig <- bind_rows(
  read.csv(file.path(stats_dir, "MannWhitney_DryWet_TotalO2.csv")) %>% transmute(Response, p),
  read.csv(file.path(stats_dir, "MannWhitney_DryWet_byKm.csv")) %>% filter(Km_mg_L == Km_base) %>%
    transmute(Response = Parameter, p)) %>%
  transmute(Treatment = "All treatments", Response, significant = p < 0.05)

tableS1_cols <- c(O2_rate_mg_kg_h = "O2 consumption (mg O2 kg-1 h-1)",
                 Vmax_mg_kg_h    = "Vmax (mg O2 kg-1 h-1)",
                 kL_L_kg_h       = "kL (L kg-1 h-1)")
for (resp in names(tableS1_cols)) {
  sig_trt <- dry_wet_sig$Treatment[dry_wet_sig$Response == resp & dry_wet_sig$significant]
  rows <- tableS1$Moisture == "Wet" & tableS1$Treatment %in% sig_trt
  tableS1[[tableS1_cols[[resp]]]][rows] <- paste0(tableS1[[tableS1_cols[[resp]]]][rows], "*")
}

tableS1 <- tableS1 %>% table_names()
write.csv(tableS1, file.path(out_dir, "TableS1_Group_Means.csv"), row.names = FALSE)

# ---- Table S2: per-vial model fits at Km = 0.2 ----------------------------------
prep <- read.csv(file.path("outputs", "modeling", "YEP_DataPrep_Log.csv"))
best <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_AllKm.csv")) %>% filter(Km_mg_L == Km_base)
mcmc <- read.csv(file.path("outputs", "modeling", "YEP_MCMC_Summary.csv")) %>% filter(Km_mg_L == Km_base)

tableS2 <- prep %>%
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
write.csv(tableS2, file.path(out_dir, "TableS2_Model_Fits_Km0.2.csv"), row.names = FALSE)

# ---- Data S1: per-vial model fits at every Km (Supplementary Data) -----------
# Same fits as Table S2, for every Km. Plain numbers (no formatting), one row
# per vial x Km. Least-squares estimates exist for every Km (01_model_DO.R);
# MCMC was run only for Km in the Sihi et al. (2020) range (0.13-1.22 mg/L),
# so the MCMC columns are empty at the other Km values.
best_all <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_AllKm.csv"))
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
stopifnot(nrow(dataS1) == nrow(best_all),
          sum(!is.na(dataS1$`kL, posterior median (L kg-1 h-1)`)) == nrow(mcmc_all))
write.csv(dataS1, file.path(out_dir, "DataS1_Model_Fits_AllKm.csv"), row.names = FALSE, na = "")

# ---- Table S4: carbon added vs water-extractable OC -------------------------
# Added solution volume per vial = total water - pore water already in the
# sediment (field-moist mass - dry mass). Added C = solution DOC x volume / dry mass.
npoc <- read_data(data_files$npoc_tn)
mass <- read.csv(data_files$water_mass) %>%
  filter(grepl("^YEP[12][ABC]_INC-[HSU]", Sample_Name)) %>%
  mutate(Moisture  = ifelse(substr(Sample_Name, 4, 4) == "1", "Dry", "Wet"),
         Treatment = dplyr::recode(sub(".*_INC-([HSU]).*", "\\1", Sample_Name),
                                   S = "Control", U = "Unburned + DOC", H = "High Burn + DOC"),
         Pore_water_mL      = as.numeric(Field_Moist_Sediment_Mass_g) - as.numeric(Dry_Sediment_Mass_g),
         Added_solution_mL  = as.numeric(Water_Mass_g) - Pore_water_mL)

# Treatment-solution DOC (OCN-SK*, -UK*, -HK* = solutions as added; mg C/L)
solution_doc <- npoc %>%
  filter(grepl("^YEP_OCN-[SUH]K", Sample_Name)) %>%
  mutate(Treatment = dplyr::recode(sub("^YEP_OCN-([SUH]).*", "\\1", Sample_Name),
                                   S = "Control", U = "Unburned + DOC", H = "High Burn + DOC"),
         DOC = flag_value(Extractable_NPOC_mg_per_L)) %>%
  group_by(Treatment) %>% summarize(Solution_DOC_mg_L = mean(DOC))

# Water-extractable OC at collection (YEP1_/YEP2_SOC) and after pre-incubation (YEP1A-C/YEP2A-C SOC-1..3)
weoc <- npoc %>%
  filter(grepl("^YEP[12][ABC]?_SOC-[123]$", Sample_Name)) %>%
  mutate(Moisture = ifelse(substr(Sample_Name, 4, 4) == "1", "Dry", "Wet"),
         Stage    = ifelse(grepl("^YEP[12]_", Sample_Name), "Field", "Pre-incubation"),
         WEOC     = flag_value(Extractable_NPOC_mg_per_kg)) %>%
  group_by(Moisture, Stage) %>% summarize(WEOC_mg_kg = mean(WEOC), .groups = "drop") %>%
  pivot_wider(names_from = Stage, values_from = WEOC_mg_kg, names_prefix = "WEOC_mg_kg_")

tableS4 <- mass %>%
  left_join(solution_doc, by = "Treatment") %>%
  mutate(Added_C_mg_kg = Solution_DOC_mg_L * (Added_solution_mL / 1000) / (as.numeric(Dry_Sediment_Mass_g) / 1000)) %>%
  group_by(Moisture, Treatment) %>%
  summarize(Solution_DOC_mg_L = mean(Solution_DOC_mg_L), Added_solution_mL = mean(Added_solution_mL),
            Added_C_mg_kg = mean(Added_C_mg_kg), .groups = "drop") %>%
  left_join(weoc, by = "Moisture") %>%
  mutate(Added_C_pct_of_field_WEOC          = 100 * Added_C_mg_kg / `WEOC_mg_kg_Field`,
         Added_C_pct_of_preincubation_WEOC  = 100 * Added_C_mg_kg / `WEOC_mg_kg_Pre-incubation`) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Unburned + DOC", "High Burn + DOC"))) %>%
  arrange(Moisture, Treatment) %>%
  mutate(across(where(is.numeric), ~ round(., 1))) %>%
  table_names()
write.csv(tableS4, file.path(out_dir, "TableS4_Carbon_Added_vs_WEOC.csv"), row.names = FALSE)

# ---- Table S5: treatment solutions ------------------------------------------
# DOC/TN: solutions as added (OCN-SK*, -UK*, -HK*). Ions: water-only incubations
# (ION-WS*, -WU*, -WH*), i.e. the same solutions without sediment.
ions <- read_data(data_files$ions)
ion_cols <- c("00940_Cl_mg_per_L", "00945_SO4_mg_per_L_as_SO4", "00618_NO3_mg_per_L_as_N",
              "00915_Ca_mg_per_L", "00925_Mg_mg_per_L", "00935_K_mg_per_L", "00930_Na_mg_per_L")

solution_ions <- ions %>%
  filter(grepl("^YEP_ION-W[SUH]", Sample_Name)) %>%
  mutate(Treatment = dplyr::recode(sub("^YEP_ION-W([SUH]).*", "\\1", Sample_Name),
                                   S = "Control", U = "Unburned + DOC", H = "High Burn + DOC")) %>%
  mutate(across(all_of(ion_cols), flag_value)) %>%
  group_by(Treatment) %>% summarize(across(all_of(ion_cols), mean))

solution_tn <- npoc %>%
  filter(grepl("^YEP_OCN-[SUH]K", Sample_Name)) %>%
  mutate(Treatment = dplyr::recode(sub("^YEP_OCN-([SUH]).*", "\\1", Sample_Name),
                                   S = "Control", U = "Unburned + DOC", H = "High Burn + DOC"),
         TN = flag_value(Extractable_TN_mg_per_L)) %>%
  group_by(Treatment) %>% summarize(TN_mg_L = mean(TN, na.rm = TRUE))

tableS5 <- solution_doc %>% left_join(solution_tn, by = "Treatment") %>% left_join(solution_ions, by = "Treatment") %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Unburned + DOC", "High Burn + DOC"))) %>%
  arrange(Treatment) %>% mutate(across(where(is.numeric), ~ signif(., 3))) %>%
  table_names()
write.csv(tableS5, file.path(out_dir, "TableS5_Treatment_Solutions.csv"), row.names = FALSE)

# ---- Table S3: sediment extracts at collection vs after pre-incubation -------
# WEOC: mean of 3 extracts. Ions: one composite extract per date (n = 1).
extract_ions <- ions %>%
  filter(grepl("^YEP[12][ABC]?_SIN-1$", Sample_Name)) %>%
  transmute(Sample = sub("_SIN-1$", "", Sample_Name),
            NO3_N_mg_L = flag_value(`00618_NO3_mg_per_L_as_N`),
            SO4_mg_L   = flag_value(`00945_SO4_mg_per_L_as_SO4`))
extract_weoc <- npoc %>%
  filter(grepl("^YEP[12][ABC]?_SOC-[123]$", Sample_Name)) %>%
  mutate(Sample = sub("_SOC-[123]$", "", Sample_Name), WEOC = flag_value(Extractable_NPOC_mg_per_kg)) %>%
  group_by(Sample) %>% summarize(WEOC_mg_kg = mean(WEOC), WEOC_sd = sd(WEOC))

tableS3 <- extract_weoc %>% left_join(extract_ions, by = "Sample") %>%
  mutate(Sediment = ifelse(substr(Sample, 4, 4) == "1", "Dry", "Wet"),
         Stage = ifelse(nchar(Sample) == 4, "Field (collection)",
                        paste0("Pre-incubated, incubation day ", match(substr(Sample, 5, 5), LETTERS))),
         across(where(is.numeric), ~ signif(., 3))) %>%
  select(Sediment, Stage, Sample, WEOC_mg_kg, WEOC_sd, NO3_N_mg_L, SO4_mg_L) %>%
  arrange(Sediment, Sample)
write.csv(tableS3, file.path(out_dir, "TableS3_PreIncubation_Change.csv"), row.names = FALSE)

print(tableS1, width = 200)
print(tableS4)
print(tableS5)
print(tableS3)
