# ==============================================================================
# 04_pool_sankey_figure.R
# Vasopressor Escalation in Refractory Distributive Shock - CLIF Consortium
# Coordinating site: alluvial figure of complete paths (poster figure)
#
# Input: PATHS_POOLED from 00_pool_load.R (one row per unique path, summed
#   across sites). Each band follows the same encounters from T0 to 48 h.
#   States at 12/24/36 h = most severe state reached in the preceding 12 h.
# ==============================================================================

# Requires: 00_pool_load.R

library(ggplot2)
library(ggalluvial)

if (nrow(PATHS_POOLED) == 0) {
  stop("PATHS_POOLED is empty. Check sites/{site}/upload_to_box/sankey_paths_{site}.csv.", call. = FALSE)
}

dir_fig = here("output", "figures")
if (!dir.exists(dir_fig)) dir.create(dir_fig, recursive = TRUE)

message("\n== Alluvial figure ==")

# ==============================================================================
# STEP 1: Factor levels (top to bottom) and labels
# ==============================================================================

state_levels = c(
  "S3",
  "S2b",
  "S2a",
  "S0",
  "S1",
  "W",
  "D"
)

outcome_levels = c(
  "alive_48h_survived",
  "alive_48h_died_later",
  "dead_48h"
)

outcome_labels = c(
  "Survived",
  "Died after 48 h",
  "Died by 48 h"
)

state_cols = c(
  "state_0h",
  "state_12h",
  "state_24h",
  "state_36h"
)

# ==============================================================================
# STEP 2: Plot data
# ==============================================================================

plot_dt = copy(PATHS_POOLED)

for (v in state_cols) {
  plot_dt[, (v) := factor(get(v), levels = state_levels)]
}

plot_dt[, outcome_48h := factor(
  outcome_48h,
  levels = outcome_levels,
  labels = outcome_labels
)]

unexpected = plot_dt[
  is.na(state_0h) | is.na(state_12h) | is.na(state_24h) | is.na(state_36h) | is.na(outcome_48h)
]

if (nrow(unexpected) > 0) {
  print(unexpected)
  stop("Paths contain unexpected state or outcome values (shown above).", call. = FALSE)
}

n_total = sum(plot_dt$n)

outcome_summary = plot_dt[, .(n = sum(n)), by = outcome_48h][order(outcome_48h)]
outcome_summary[, pct := round(n / n_total * 100, 1)]

message("  Encounters: ", format_n(n_total), " | unique paths: ", format_n(nrow(plot_dt)))
print(outcome_summary)

# ==============================================================================
# STEP 3: Figure
# ==============================================================================
# Bands colored by 48-h outcome: one-hue ordinal ramp (validated, light mode)
# from survived (light) to died by 48 h (dark).
# geom_flow + aes.bind = "flows": between each pair of time points, encounters
# with the same from-state, to-state, and outcome form ONE band, and bands of
# the same outcome sit together inside each state. (geom_alluvium drew every
# unique path as its own ribbon, which produced the striped T0 column.)

outcome_fill = c(
  "Survived"        = "#86b6ef",
  "Died after 48 h" = "#256abf",
  "Died by 48 h"    = "#0d366b"
)

text_primary   = "#0b0b0b"
text_secondary = "#52514e"
stratum_fill   = "#f0efec"

fig_caption = paste(
  "States: S0, vasopressin + one adrenergic agent; S1, one agent;",
  "S2a, two or more adrenergic agents; S2b, angiotensin II + norepinephrine and/or vasopressin;",
  "S3, four or more agents, or angiotensin II + epinephrine, phenylephrine, or dopamine;",
  "W, no vasopressor (includes discharged alive); D, died or discharged to hospice.",
  "Outcomes: survived = discharged alive; died after 48 h = in-hospital death or hospice after 48 h.",
  "States at 12, 24, and 36 h are the most severe state reached in the preceding 12 h.",
  sep = "\n"
)

fig = ggplot(
  plot_dt,
  aes(
    axis1 = state_0h,
    axis2 = state_12h,
    axis3 = state_24h,
    axis4 = state_36h,
    axis5 = outcome_48h,
    y     = n
  )
) +
  # only arguments confirmed to render on ggalluvial 0.12.6: `color = NA` (and
  # possibly the curve arguments) made ggplot drop every flow as missing
  geom_flow(
    aes(fill = outcome_48h),
    width    = 0.3,
    alpha    = 0.9,
    aes.bind = "flows"
  ) +
  geom_stratum(
    width     = 0.3,
    fill      = stratum_fill,
    color     = text_secondary,
    linewidth = 0.3
  ) +
  geom_text(
    stat  = "stratum",
    aes(label = ifelse(
      after_stat(prop) >= 0.03,
      paste0(after_stat(stratum), "\n", round(after_stat(prop) * 100), "%"),
      NA
    )),
    size       = 3,
    lineheight = 0.9,
    color = text_primary,
    na.rm = TRUE
  ) +
  # continuous x (axes sit at 1:5): a discrete scale with limits drops the
  # geom_flow layer, which left only the strata in the last version
  scale_x_continuous(
    breaks = 1:5,
    labels = c("T0", "12 h", "24 h", "36 h", "48 h"),
    expand = expansion(add = 0.25)
  ) +
  scale_fill_manual(
    values = outcome_fill,
    name   = NULL
  ) +
  labs(
    title    = "Vasopressor state over the first 48 hours of refractory distributive shock",
    subtitle = paste0("Encounters, n = ", format_n(n_total), "; all encounters start in S0 (norepinephrine + vasopressin); bands colored by 48-hour outcome"),
    x        = NULL,
    y        = NULL,
    caption  = fig_caption
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(color = text_primary, size = 11),
    plot.title      = element_text(color = text_primary, face = "bold"),
    plot.subtitle   = element_text(color = text_secondary),
    plot.caption    = element_text(color = text_secondary, hjust = 0, size = 7.5, lineheight = 1.1),
    legend.position = "top",
    legend.text     = element_text(color = text_primary),
    plot.background = element_rect(fill = "#fcfcfb", color = NA)
  )

# ==============================================================================
# STEP 4: Save
# ==============================================================================

ggsave(
  filename = file.path(dir_fig, paste0("alluvial_paths_", today, ".png")),
  plot     = fig,
  width    = 11,
  height   = 7,
  dpi      = 300
)

ggsave(
  filename = file.path(dir_fig, paste0("alluvial_paths_", today, ".pdf")),
  plot     = fig,
  width    = 11,
  height   = 7
)

fwrite(outcome_summary, file.path(dir_fig, paste0("alluvial_outcome_summary_", today, ".csv")))

message("  Saved: output/figures/alluvial_paths_", today, ".png and .pdf")
message("\n== Alluvial figure complete ==")
