# =============================================================================
# 04_fig1_km_sensitivity.R -- Fig. 1 and the full-range SI figure
# -----------------------------------------------------------------------------
# Vmax and kL (per kg dry sediment) against the fixed Km used in the fit.
#   Fig. 1        : Km within the Sihi et al. (2020) posterior range (0.13-1.22 mg/L)
#   Fig. S (full) : Km 0.05-10 mg/L, including Km = 10 mg/L used in the
#                   submitted manuscript (for the response to reviewers)
#
# Four panels: rows = parameter (Vmax, kL), columns = sediment history
# (previously dry, wet). Every panel has its own y axis.
#   Stars   (Fig. S only, top of each panel): previously dry vs wet at that Km,
#           DOM treatments pooled, Mann-Whitney test (03_stats.R). The same
#           star is repeated in the dry and the wet panel of a row.
#           *** p < 0.001, ** p < 0.01, * p < 0.05, ns p >= 0.05.
#           Fig. 1 has no stars: dry and wet differ at every Km it shows
#           (p < 0.001), which is stated in the caption instead. The script
#           stops if that is ever not true.
#   Letters (above the points): DOM treatments within that panel. Shown only
#           where the Kruskal-Wallis test among the three treatments is
#           significant (p < 0.05) at that Km. Treatments that share a letter
#           do not differ (Dunn's test, no p adjustment). Letters are in the
#           same left-to-right order as the points (control, unburned, burned).
#           No letters = Kruskal-Wallis not significant.
# Large symbols: median; bars: interquartile range; small symbols: single vials.
#
# Inputs : outputs/modeling/YEP_BestFit_AllKm.csv (01_model_DO.R)
#          outputs/stats/MannWhitney_DryWet_byKm.csv (03_stats.R)
#          outputs/stats/Treatment_Letters_byKm.csv  (03_stats.R)
# Outputs: outputs/figures/Fig1_Km_Sensitivity_SihiRange.pdf/.png
#          outputs/figures/FigS_Km_Sensitivity_FullRange.pdf/.png
# =============================================================================

rm(list = ls(all = TRUE))
library(dplyr)
library(tidyr)
library(ggplot2)
source("src/00_plot_style.R")   # treatment names, colors, shapes

fig_dir <- file.path("outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sihi_range   <- c(0.13, 1.22)  # Sihi et al. (2020) posterior range, mg/L

# Panel titles (plotmath): sediment history on the first line, parameter and
# units on the second
history_titles   <- c(Dry = "bold('Previously dry sediments')", Wet = "bold('Wet sediments')")
parameter_titles <- c(Vmax_mg_kg_h = "V[max]~(mg~O[2]~kg^{-1}~h^{-1})",
                      kL_L_kg_h    = "k[L]~(L~kg^{-1}~h^{-1})")
panel_title <- function(parameter, moisture) {
  paste0("atop(", history_titles[moisture], ",", parameter_titles[parameter], ")")
}
# Order of the panels: Vmax dry, Vmax wet, kL dry, kL wet
panel_order <- c(panel_title("Vmax_mg_kg_h", "Dry"), panel_title("Vmax_mg_kg_h", "Wet"),
                 panel_title("kL_L_kg_h", "Dry"),    panel_title("kL_L_kg_h", "Wet"))

# Text sizes (ggplot text size is in mm; 5 is about 14 pt)
base_size   <- 16
star_size   <- 5
letter_size <- 5

# Horizontal spacing between the three treatments at the same Km (log10 units of Km)
offset_treatment <- 0.030   # points
offset_letters   <- 0.036   # letters above the points

# ---- 1. Data ------------------------------------------------------------------
best_fits <- read.csv(file.path("outputs", "modeling", "YEP_BestFit_AllKm.csv"))
mw        <- read.csv(file.path("outputs", "stats", "MannWhitney_DryWet_byKm.csv"))
letters   <- read.csv(file.path("outputs", "stats", "Treatment_Letters_byKm.csv"))

# One row per vial x Km x parameter, with its horizontal position
vials <- best_fits %>%
  select(Sample_Name, Moisture, Treatment, Km_mg_L, Vmax_mg_kg_h, kL_L_kg_h) %>%
  mutate(Treatment = as.character(relabel_treatment(Treatment))) %>%
  pivot_longer(c(Vmax_mg_kg_h, kL_L_kg_h), names_to = "Parameter", values_to = "Value") %>%
  mutate(Panel = panel_title(Parameter, Moisture),
         x = Km_mg_L * 10^((match(Treatment, treatment_labels) - 2) * offset_treatment))

# Median and interquartile range of the 5 vials in each group
groups <- vials %>%
  group_by(Panel, Parameter, Moisture, Treatment, Km_mg_L, x) %>%
  summarize(median = median(Value), q25 = quantile(Value, 0.25), q75 = quantile(Value, 0.75),
            .groups = "drop")

# Fixed factor orders for every data frame
set_levels <- function(d) {
  d <- d %>% mutate(Panel = factor(Panel, levels = panel_order))
  if ("Treatment" %in% names(d)) d$Treatment <- factor(d$Treatment, levels = unname(treatment_labels))
  if ("Moisture" %in% names(d))  d$Moisture  <- factor(d$Moisture, levels = names(moisture_shapes))
  d
}

# ---- 2. Plot function (used for both figures) -------------------------------
# compact = TRUE (many Km values close together): stars and letters are written
# vertically, and the three letters of a Km become one label ("a a b").
# show_stars = FALSE leaves out the dry-vs-wet stars (allowed only when every
# Km in the figure has p < 0.001, so the caption can say so).
km_plot <- function(km_min, km_max, shade_sihi, compact = FALSE, show_stars = TRUE) {
  in_range <- function(d) d %>% filter(Km_mg_L >= km_min, Km_mg_L <= km_max)
  v <- in_range(vials)
  g <- in_range(groups)

  # Top and height of each panel, used to place letters and stars above the data
  panels <- v %>% group_by(Panel, Parameter, Moisture) %>%
    summarize(panel_top = max(Value), panel_range = diff(range(Value)), .groups = "drop") %>%
    mutate(panel_range = ifelse(panel_range == 0, 1, panel_range))
  line_height <- if (compact) 0.30 else 0.09   # one line of text, as a fraction of the panel height

  # Letters: DOM treatments, in one row above all vials at that Km
  km_top <- v %>% group_by(Panel, Km_mg_L) %>% summarize(km_top = max(Value), .groups = "drop")
  l <- in_range(letters %>% filter(!is.na(Km_mg_L))) %>%
    mutate(Panel = panel_title(Response, Moisture),
           x = Km_mg_L * 10^((match(Treatment, treatment_labels) - 2) * offset_letters)) %>%
    left_join(km_top, by = c("Panel", "Km_mg_L")) %>%
    left_join(panels %>% select(Panel, panel_top, panel_range), by = "Panel")
  l <- l %>% mutate(y = km_top + 0.07 * panel_range)
  # Two Km values that sit very close together (1.00 and 1.22): if their letters
  # would land at about the same height, lift the second set by one line so the
  # letters do not run into each other
  if (!compact) {
    heights <- l %>% distinct(Panel, Km_mg_L, y, panel_range) %>% arrange(Panel, Km_mg_L)
    for (i in seq_len(nrow(heights))[-1]) {
      same_panel <- heights$Panel[i] == heights$Panel[i - 1]
      close_km   <- log10(heights$Km_mg_L[i] / heights$Km_mg_L[i - 1]) < 0.11
      one_line   <- line_height * heights$panel_range[i]
      if (same_panel && close_km && abs(heights$y[i] - heights$y[i - 1]) < one_line) {
        heights$y[i] <- heights$y[i - 1] + one_line
      }
    }
    l <- l %>% select(-y) %>% left_join(heights %>% select(Panel, Km_mg_L, y), by = c("Panel", "Km_mg_L"))
  }
  if (compact) {
    l <- l %>% arrange(match(Treatment, treatment_labels)) %>%
      group_by(Panel, Km_mg_L, y, panel_top, panel_range) %>%
      summarize(Letter = paste(Letter, collapse = " "), x = unique(Km_mg_L), .groups = "drop")
  }

  if (!show_stars) {
    mw_shown <- in_range(mw)
    if (any(is.na(mw_shown$p) | mw_shown$p >= 0.001)) {
      stop("Dry vs wet is not p < 0.001 at every Km in this figure: use show_stars = TRUE")
    }
  }

  # Stars: dry vs wet, in a row above everything else in the panel
  highest <- l %>% group_by(Panel) %>% summarize(letters_top = max(y), .groups = "drop")
  s <- in_range(mw) %>%
    mutate(label = ifelse(is.na(stars), "", stars)) %>%
    tidyr::crossing(Moisture = names(moisture_shapes)) %>%
    mutate(Panel = panel_title(Parameter, Moisture)) %>%
    left_join(panels %>% select(Panel, panel_top, panel_range), by = "Panel") %>%
    left_join(highest, by = "Panel") %>%
    mutate(y = pmax(panel_top, letters_top + panel_range * line_height, na.rm = TRUE) + 0.08 * panel_range)
  if (!show_stars) s <- s[0, ]

  text_angle <- if (compact) 90 else 0
  text_hjust <- if (compact) 0 else 0.5      # vertical text starts at y and reads upward
  top_space  <- if (compact) 0.12 else 0.06  # empty room above the stars

  v <- set_levels(v); g <- set_levels(g); l <- set_levels(l); s <- set_levels(s)

  p <- ggplot(g, aes(x, median, color = Treatment, shape = Moisture,
                     group = interaction(Treatment, Moisture)))
  if (shade_sihi) {
    p <- p + annotate("rect", xmin = sihi_range[1], xmax = sihi_range[2],
                      ymin = -Inf, ymax = Inf, fill = "gray85", alpha = 0.5)
  }
  p +
    geom_point(data = v, aes(y = Value), size = 1.5, alpha = 0.35) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, alpha = 0.7, linewidth = 0.7) +
    geom_line(aes(linetype = Moisture), linewidth = 0.7) +
    geom_point(size = 3.4) +
    geom_text(data = s, aes(x = Km_mg_L, y = y, label = label), inherit.aes = FALSE,
              size = star_size, angle = text_angle, hjust = text_hjust) +
    geom_text(data = l, aes(x = x, y = y, label = Letter), inherit.aes = FALSE,
              size = letter_size, fontface = "bold", angle = text_angle, hjust = text_hjust) +
    facet_wrap(~ Panel, ncol = 2, scales = "free_y", labeller = label_parsed) +
    scale_x_log10(breaks = sort(unique(g$Km_mg_L))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, top_space))) +
    scale_color_manual(values = treatment_colors) +
    scale_shape_manual(values = moisture_shapes, guide = "none") +
    scale_linetype_manual(values = moisture_lines, guide = "none") +
    labs(x = expression(K[m]~(mg~L^{-1}*","~log~scale)), y = NULL, color = "DOM treatment") +
    theme_bw(base_size = base_size) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = base_size), legend.title = element_text(size = base_size),
          strip.text = element_text(size = base_size - 1, margin = margin(5, 0, 5, 0)),
          axis.text = element_text(size = base_size - 2, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = if (compact) base_size - 5 else base_size - 2),
          axis.title = element_text(size = base_size + 1),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4, linewidth = 1)))
}

# ---- 3. Fig. 1: Sihi et al. (2020) posterior range --------------------------
fig1 <- km_plot(sihi_range[1], sihi_range[2], shade_sihi = FALSE, show_stars = FALSE)
ggsave(file.path(fig_dir, "Fig1_Km_Sensitivity_SihiRange.pdf"), fig1, width = 12, height = 10)
ggsave(file.path(fig_dir, "Fig1_Km_Sensitivity_SihiRange.png"), fig1, width = 12, height = 10, dpi = 300)

# ---- 4. Fig. S: full Km range, including the submitted value -----------------
figS <- km_plot(0, Inf, shade_sihi = TRUE, compact = TRUE)
ggsave(file.path(fig_dir, "FigS_Km_Sensitivity_FullRange.pdf"), figS, width = 14, height = 10.5)
ggsave(file.path(fig_dir, "FigS_Km_Sensitivity_FullRange.png"), figS, width = 14, height = 10.5, dpi = 300)
