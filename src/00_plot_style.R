# =============================================================================
# 00_plot_style.R -- labels, colors and shapes shared by every figure
# -----------------------------------------------------------------------------
# Figure scripts start with:  source("src/00_plot_style.R")
#
# This script was created to rename the treatment categories in the figures.
# The data files and model outputs keep the original treatment codes
# ("Unburned + DOC", "High Burn + DOC"). Figures and statistics tables show the
# names used in the manuscript text. To rename a treatment everywhere, change
# ONLY the right-hand side of treatment_labels below.
# =============================================================================

# code in the data  =  name shown in figures and tables
treatment_labels <- c("Control"         = "Control",
                      "Unburned + DOC"  = "Unburned DOM",
                      "High Burn + DOC" = "Burned DOM")

# Colors follow the displayed names, in the same order as treatment_labels
treatment_colors <- setNames(c("#0072B2", "darkgreen", "#8B4513"), treatment_labels)

moisture_shapes <- c("Dry" = 16, "Wet" = 17)
moisture_lines  <- c("Dry" = "dashed", "Wet" = "solid")

# Convert treatment codes to the displayed names (ordered factor: Control first)
relabel_treatment <- function(x) {
  factor(treatment_labels[as.character(x)], levels = unname(treatment_labels))
}
