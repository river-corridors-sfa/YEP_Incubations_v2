# =============================================================================
# 05_fig_total_O2.R -- measured (total) O2 consumption by DOM treatment
# -----------------------------------------------------------------------------
# Measured O2 consumption rate per kg dry sediment for each DOM treatment,
# within previously dry and within wet sediments. This rate comes from the
# DO time series directly, so it does not depend on Km.
#   Text in each panel: Kruskal-Wallis test among the three treatments.
#   Letters: shown only where Kruskal-Wallis is significant (p < 0.05).
#            Treatments that share a letter do not differ (Dunn's test, no p
#            adjustment).
# Large symbols: median; bars: interquartile range; small symbols: single vials.
#
# Inputs : outputs/stats/microcosm_responses.csv           (03_stats.R)
#          outputs/stats/KruskalWallis_Treatment_byKm.csv  (03_stats.R)
#          outputs/stats/Treatment_Letters_byKm.csv        (03_stats.R)
# Outputs: outputs/figures/Fig_Total_O2_by_Treatment.pdf/.png
# =============================================================================

rm(list = ls(all = TRUE))
library(dplyr)
library(ggplot2)
source("src/00_plot_style.R")   # treatment names, colors, shapes

fig_dir <- file.path("outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

history_titles <- c(Dry = "Previously dry sediments", Wet = "Wet sediments")
base_size <- 16

# ---- 1. Data ------------------------------------------------------------------
vials <- read.csv(file.path("outputs", "stats", "microcosm_responses.csv")) %>%
  transmute(ID, Moisture, Treatment = relabel_treatment(Treatment), Value = O2_rate_mg_kg_h)

kw <- read.csv(file.path("outputs", "stats", "KruskalWallis_Treatment_byKm.csv")) %>%
  filter(Response == "O2_rate_mg_kg_h")
letters <- read.csv(file.path("outputs", "stats", "Treatment_Letters_byKm.csv")) %>%
  filter(Response == "O2_rate_mg_kg_h")

# Median and interquartile range of the 5 vials in each group
groups <- vials %>%
  group_by(Moisture, Treatment) %>%
  summarize(median = median(Value), q25 = quantile(Value, 0.25), q75 = quantile(Value, 0.75),
            top = max(Value), .groups = "drop")

# Height of each panel, used to place the text above the data
panels <- vials %>% group_by(Moisture) %>%
  summarize(panel_top = max(Value), panel_range = diff(range(Value)), .groups = "drop")

# Kruskal-Wallis result, top left of each panel
kw_labels <- kw %>% left_join(panels, by = "Moisture") %>%
  mutate(label = sprintf("Kruskal-Wallis p = %.3f", p), y = panel_top + 0.22 * panel_range)

# Letters above each treatment (only where Kruskal-Wallis is significant)
letter_labels <- letters %>%
  mutate(Treatment = factor(Treatment, levels = unname(treatment_labels))) %>%
  left_join(groups %>% select(Moisture, Treatment, top), by = c("Moisture", "Treatment")) %>%
  left_join(panels, by = "Moisture") %>%
  mutate(y = top + 0.07 * panel_range)

# ---- 2. Plot ------------------------------------------------------------------
set.seed(1)   # same horizontal jitter of the small points every run
fig <- ggplot(groups, aes(Treatment, median, color = Treatment, shape = Moisture)) +
  geom_jitter(data = vials, aes(y = Value), width = 0.12, height = 0, size = 2.2, alpha = 0.45) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, linewidth = 0.8) +
  geom_point(size = 4.5) +
  geom_text(data = letter_labels, aes(y = y, label = Letter), color = "black",
            size = 6, fontface = "bold", show.legend = FALSE) +
  geom_text(data = kw_labels, aes(x = 0.6, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, size = 5) +
  facet_wrap(~ Moisture, scales = "free_y", labeller = labeller(Moisture = history_titles)) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  scale_shape_manual(values = moisture_shapes, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(x = NULL, y = expression(Measured~O[2]~consumption~(mg~O[2]~kg^{-1}~h^{-1}))) +
  theme_bw(base_size = base_size) +
  theme(strip.text = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size - 2, color = "black"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "Fig_Total_O2_by_Treatment.pdf"), fig, width = 11, height = 5.5)
ggsave(file.path(fig_dir, "Fig_Total_O2_by_Treatment.png"), fig, width = 11, height = 5.5, dpi = 300)
