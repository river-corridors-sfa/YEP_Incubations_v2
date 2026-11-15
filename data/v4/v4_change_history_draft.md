# YEP data package v4: change history draft and file status (for Brie)

Corrected source files are in the YEP_Incubations repo (github.com/river-corridors-sfa/YEP_Incubations), folder `data/working/`. They carry the names the processing scripts produced, not data package names. Full tracker: `data/v4/CHANGES.md`.

**Version 4, [Month 2026]**

MODIFIED FILES

**YEP_Sediment_Water_Mass_Volume.csv (v3)**
• Corrected Water_Mass_g for incubation (INC) samples. Water mass was previously calculated assuming 30 mL of treatment solution was added to each vial. Vials were filled until there was no headspace, so water mass is now calculated as the measured brim-full vial volume (41.56 mL, with septum cap and probe) minus the volume of sediment grains (dry sediment mass ÷ 2.65 g/cm³). Dry sediment masses did not change.
- Status: Done (Vanessa)
- Source file: `data/working/YEP_Merged_Wet_and_Dry_Mass_2026-10-30.csv`
- Formatting needed: drop the STR rows (not in the package); add IGSN, Material and header rows.

**YEP_Sediment_Respiration_Raw_Dissolved_Oxygen_Temp.csv (v2)**
• Corrected start times for YEP2A_INC-H1 and YEP2A_INC-U1. The time series now begin at the start of the incubation; previously the first 3–4 minutes were missing.
- Status: Done (Vanessa)
- Source file: `data/working/Merged_DO_Firesting_Data_2026-10-07.csv` (corrected start times are in `data/working/Merged_Firesting_Mapping_2026-10-07.csv`)
- Formatting needed: only YEP2A_INC-H1 (+42 rows) and YEP2A_INC-U1 (+30 rows) differ from v3; add IGSN, Material, Methods_Deviation and header rows.

**YEP_Sediment_Incubations_Respiration_Rates_pH_Temp.csv (v3)**
• Corrected respiration rates and regression statistics for YEP2A_INC-H1 and YEP2A_INC-U1 using the corrected time series.
• Corrected normalized respiration rates (per kg dry sediment) for all incubation samples using the corrected water masses.
- Status: Done (Vanessa)
- Source file: `data/working/YEP_Respiration_10-30-26.csv`
- Formatting needed: carry over the pH and Temperature_degC columns from v3 (unchanged); add IGSN, Material and header rows; use the package formats (R² of 1 → ">0.99", p of 0 → "<0.001").

**YEP_Sediment_Incubations_DO_vs_Incubation_Time_Plots.pdf (v2)**
• Updated plots for YEP2A_INC-H1 and YEP2A_INC-U1.
- Status: Done (Vanessa)
- Source: individual plots in Vanessa's SharePoint plots folder; combine into one PDF.

**YEP_Sediment_NPOC_TN.csv (v2)**
• Corrected Extractable_NPOC_mg_per_kg and Extractable_TN_mg_per_kg for post-incubation samples using the corrected water masses. Concentrations in mg/L and pre-incubation extractions did not change.
- Status: Done (Vanessa)
- Source file: `data/working/YEP_Sediment_NPOC_Mass_2026-10-30.csv`
- Formatting needed: drop the Randomized_ID column; add IGSN, Material and header rows. Flags and Methods_Deviation already match v3.

**YEP_Sediment_CO2_CH4_N2O.csv (v4)**
• Removed the CO₂, CH₄, and N₂O production rate columns. Rates were calculated as end-point concentration divided by incubation duration and were sensitive to incubation length. Partial pressures and dissolved concentrations did not change.
- Status: To do (Brie)
- Source file: v3 `v3_YEP_Sediment_CO2_CH4_N2O.csv`; delete the three Rate_* columns.

**YEP_Methods_Codes.csv (v4)**
• Updated DWA_T_DP_003 to describe the revised water mass calculation.
• Updated GAS_T_DP_000 to remove the production rate calculation, describe the solubility constants consistently with the readme, and correct the concentration units to moles per liter.
- Status: Done (Vanessa)
- To check (Brie): whether a revised method keeps its ID or gets a new one (e.g., DWA_T_DP_004).

**YEP_flmd.csv (v4)**
• Updated file names and version numbers for modified files.
- Status: To do (Brie)

**YEP_dd.csv (v4)**
• Removed rows for the gas production rate columns.
- Status: To do (Brie)

**readme_YEP.pdf (v4)**
• Updated versioning information and the version number in the data package title.
• Clarified that the 6-hour incubation limit applied only to treatment-solution-only incubations.
• Clarified that replicates H1 and U2 were incubated on Day 2 for dry sediments only. *(added)*
• Added that two dry-sediment Day 1 incubations (YEP1A_INC-H1 and YEP1A_INC-U2) were not retained and were repeated on Day 2. *(added)*
• Added that YEP2A_SIR-S2 and YEP2B_SIR-S3 had low signal at 0.075 s ion accumulation and were reacquired at 0.2 s; only the 0.2 s spectra were processed through CoreMS and included in the processed data files; both sets of .xml files are provided. *(added)*
- Status: To do (Brie)
