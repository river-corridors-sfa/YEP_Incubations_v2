# data/v4: corrected files for the YEP data package v4

This folder holds **only** the files that change between the published v3 package
(`data/v3/`, doi:10.15485/3001969) and v4, under the names the processing scripts
give them. Unchanged inputs are read from `data/v3/` (see `src/00_data_paths.R`).
When v4 is published, download it, point `src/00_data_paths.R` at it, and archive
`data/v3/` and this folder.

| # | What changes | Why | Script | File the script produces | Replaces in package (v3 → v4 name) | Status |
|---|---|---|---|---|---|---|
| 1 | `Water_Mass_g` for the 30 incubation (INC) rows; dry masses unchanged | Vials were topped up to the brim, not filled with exactly 30 mL. New water = brim-full vial volume − grain volume (dry mass ÷ 2.65) | Mass script (`vial_brim_volume_mL = 41.56`) | `YEP_Merged_Wet_and_Dry_Mass_<YYYY-MM-DD>.csv` | `v2_YEP_Sediment_Water_Mass_Volume.csv` → `v3_` | Done: `data/working/YEP_Merged_Wet_and_Dry_Mass_2026-10-30.csv` (drop STR rows) |
| 2 | Raw DO for YEP2A_INC-H1 and -U1 starts at the true incubation start (~11:08 and ~11:00); +42 and +30 rows | Start-time correction: published series began 3–4 min late | None (already made) | `Merged_DO_Firesting_Data_2026-10-07.csv` (in `data/working/`) | `YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temp.csv` → `v2_` | Ready |
| 3 | Per-L rates and regression stats for H1 (−21.38 → −26.47 mg L⁻¹ h⁻¹) and U1 (−29.61 → −40.72) | Start-time fix (#2) | `YEP_Calculate_Respiration.R` | `YEP_Respiration.csv`* | `v2_YEP_Sediment_Incubations_Respiration_Rates_pH_Temp.csv` → `v3_` (per-L columns) | Done: `data/working/YEP_Respiration_10-30-26.csv` |
| 4 | Per-kg rates for all 30 vials | New water mass (#1) | `YEP_Calculate_Normalized_Respiration.R` | `YEP_Respiration.csv`* | Same file as #3 (per-kg column) | Done: same file; add pH, Temperature_degC, IGSN, Material from v3 and package formatting (">0.99", "<0.001") |
| 5 | DO-vs-time plots for H1 and U1 | Start-time fix (#2) | `YEP_Calculate_Respiration.R` | `<Sample_Name>_segmented.pdf` (one per vial) | `YEP_Sediment_Incubations_DO_vs_Incubation_Time_Plots.pdf` → `v2_` (combined) | After #3 |
| 6 | Per-kg NPOC and TN for the 30 post-incubation samples (SOC-H/S/U); mg/L values and pre-incubation SOC-1/2/3 unchanged | New water mass (#1) | NPOC script | `YEP_Sediment_NPOC_Mass_<YYYY-MM-DD>.csv`** | `YEP_Sediment_NPOC_TN.csv` → `v2_` | Done: `data/working/YEP_Sediment_NPOC_Mass_2026-10-30.csv` (flags preserved) |
| 7 | Readme: change history for #1–6 and #8–9; the "6 hours" stop applies only to water-only controls; Note 5 (H1/U2 on Day 2) applies to dry sediment only; Day 1 dry H1/U2 were dropped and repeated on Day 2; YEP2A_SIR-S2 and YEP2B_SIR-S3 had low signal at 0.075 s and were reacquired at 0.2 s (only the 0.2 s spectra were processed; both .xml sets provided) | Accuracy | — | — | `v3_readme_YEP.pdf` → `v4_` | After #1–6 |
| 8 | Rewrite method descriptions `DWA_T_DP_003` (water mass) and `GAS_T_DP_000` (remove rate calculation) | Text describes the old 30 mL calculation and the removed rates | — | — | `v3_YEP_Methods_Codes.csv` → `v4_` | Done (Vanessa) |
| 9 | Remove CO₂, CH₄ and N₂O rate columns, and their rows in the data dictionary (`v3_YEP_dd.csv` → `v4_`) | Rates are end-point ÷ incubation time, so they mostly track how long each vial ran | — | — | `v3_YEP_Sediment_CO2_CH4_N2O.csv` → `v4_` | Brie |

\* #3 and #4 both write `YEP_Respiration.csv` to the same folder; rename them (e.g., `YEP_Respiration_per_L.csv`, `YEP_Respiration_Normalized.csv`) before running.
\*\* The NPOC script currently names its output `HP_…`; set `study.code` and the prefix to `YEP`.

Not changing: gravimetric moisture, ions, FTICR files (the 0.2 s reruns are kept; see readme row), sample name metadata, field metadata, IGSN mapping, field protocol, site photos.

## Vial brim volume (measured 2026-10-30)

| Replicate | Empty vial + cap (g) | Full (g) | Capacity (mL)* |
|---|---|---|---|
| 1 | 28.54 | 70.33 | 41.87 |
| 2 | 28.83 | 69.65 | 40.90 |
| 3 | 28.37 | 70.18 | 41.89 |
| Mean ± SD | | | 41.56 ± 0.57 |

\* (full − empty) ÷ 0.998 g/mL. Empty and full weights include cap and septum. "Full" = water left after
filling to the brim, inserting the probe through the septum with the lid on (overflow released), and removing
the probe, i.e. the liquid volume present during an incubation. The probe displacement (0.65 mL, graduated
cylinder) is therefore already excluded and is not subtracted again.
Brim volume used for water mass = **41.56 mL**.

Resulting water in incubation vials (mean): dry 37.8 g (was 30.4; ×1.24), wet 36.3 g (was 36.6; ×0.99).
Solution actually added: dry ~37.4 mL, wet ~29.7 mL.
