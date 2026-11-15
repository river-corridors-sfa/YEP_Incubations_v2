# =============================================================================
# 00_data_paths.R -- single source of truth for every input file
# -----------------------------------------------------------------------------
# Every analysis script starts with:  source("src/00_data_paths.R")
# and reads data only through data_files$<name>. Paths are relative to the
# repo root (open YEP_Incubations.Rproj so the working directory is correct).
#
# Folders:
#   data/v3/       published ESS-DIVE package v3 (doi:10.15485/3001969). Never edit.
#   data/v4/       ONLY corrected files, with the names the processing scripts
#                  give them (not data-package names). See data/v4/CHANGES.md.
#   data/working/  inputs that are not in the data package
#
# When a corrected file lands in data/v4/, change ONLY its line below.
# When v4 is published: download it, point every line at it, and archive
# data/v3/ and data/v4/.
# =============================================================================

v3_root     <- file.path("data", "v3")
v3_sample   <- file.path(v3_root, "v3_YEP_Sample_Data")
v3_fticr    <- file.path(v3_sample, "FTICR")
v4_root     <- file.path("data", "v4")
working     <- file.path("data", "working")

data_files <- list(

  # ---- Files corrected for v4 (interim names; swap to published v4 names later) ----
  # Water mass, corrected (brim volume 41.56 mL); includes STR rows not in the package
  water_mass        = file.path(working, "YEP_Merged_Wet_and_Dry_Mass_2026-10-30.csv"),

  # Respiration rates, corrected (H1/U1 start-time fix; per-kg with corrected water)
  rates             = file.path(working, "YEP_Respiration_10-30-26.csv"),

  # Published v3 rates: use ONLY for the pH and Temperature_degC columns
  rates_published   = file.path(v3_sample, "v2_YEP_Sediment_Incubations_Respiration_Rates_pH_Temp.csv"),

  # NPOC/TN, corrected per-kg for post-incubation samples (flags preserved)
  npoc_tn           = file.path(working, "YEP_Sediment_NPOC_Mass_2026-10-30.csv"),

  # Raw DO with the H1/U1 start-time fix (not yet in the data package)
  raw_do            = file.path(working, "Merged_DO_Firesting_Data_2026-10-07.csv"),

  # ---- Files NOT changing in v4 ----
  ions              = file.path(v3_sample, "YEP_Sediment_Ions.csv"),
  gas               = file.path(v3_sample, "v3_YEP_Sediment_CO2_CH4_N2O.csv"),  # rate columns to be removed in v4; do not use them
  moisture          = file.path(v3_sample, "YEP_Sediment_Gravimetric_Moisture.csv"),
  methods_codes     = file.path(v3_sample, "v3_YEP_Methods_Codes.csv"),
  sample_metadata   = file.path(v3_root,   "v2_YEP_Sample_Name_Metadata.csv"),
  field_metadata    = file.path(v3_root,   "YEP_Field_Metadata.csv"),
  igsn_mapping      = file.path(v3_root,   "YEP_Metadata_IGSN_Mapping.csv"),
  fticr_data        = file.path(v3_fticr,  "YEP_Sediment_CoreMS_Processed_ICR_Data.csv"),
  fticr_mol         = file.path(v3_fticr,  "YEP_Sediment_CoreMS_Processed_ICR_Mol.csv"),
  fticr_calibration = file.path(v3_fticr,  "YEP_Sediment_CoreMS_Processed_ICR_Calibration.csv"),

  # ---- Working files (not in the data package) ----
  firesting_mapping = file.path(working, "Merged_Firesting_Mapping_2026-10-07.csv")   # H1/U1 start times corrected
)

# Stop early if any input is missing or misnamed -- never fall back silently.
missing <- names(data_files)[!file.exists(unlist(data_files))]
if (length(missing) > 0) {
  stop("Missing input file(s): ",
       paste0(missing, " -> ", unlist(data_files[missing]), collapse = "; "))
}

# Helper for ESS-DIVE formatted CSVs (metadata header rows, -9999 = missing).
read_essdive <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr   <- grep("^Field_Name,", lines)[1]
  meta  <- "^(Unit|Unit_Basis|MethodID_|Analysis_|Data_Status|#End_Data)"
  body  <- lines[(hdr + 1):length(lines)]
  body  <- sub("^#Start_Data,", "N/A,", body)
  body  <- body[!grepl(meta, body) & nzchar(body)]
  df <- utils::read.csv(text = c(lines[hdr], body), check.names = FALSE,
                        na.strings = c("-9999", "NA"), stringsAsFactors = FALSE)
  df[names(df) != "Field_Name"]
}

# Read either a plain CSV (script outputs in data/working, data/v4) or an
# ESS-DIVE formatted CSV (published package files with metadata header rows).
read_data <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  if (grepl("^#Columns", first)) read_essdive(path)
  else utils::read.csv(path, check.names = FALSE, na.strings = c("-9999", "NA"),
                       stringsAsFactors = FALSE)
}
