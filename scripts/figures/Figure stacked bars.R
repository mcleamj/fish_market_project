library(tidyverse)
library(patchwork)
library(stringr)

# -----------------------------------------------------------------------------
# 1. Setup Factor Levels and Palette
# -----------------------------------------------------------------------------
diet_levels <- c("Grazer", "Microphage", "Planktivore", "Omnivore", "Invertivore", "Piscivore")

diet_cols <- c(
  "Grazer"      = "#0034F5", # Blue
  "Microphage"  = "#34866D", # Dark Green
  "Planktivore" = "#78B214", # Light Green
  "Omnivore"    = "#E9CB2A", # Yellow
  "Invertivore" = "#FF9115", # Orange
  "Piscivore"   = "#FD3000"  # Red
)

# -----------------------------------------------------------------------------
# 2. Clean Data (Clean dots, Capitalize, Set Factor Levels)
# -----------------------------------------------------------------------------
reef_diet <- reef_diet %>%
  mutate(
    location = gsub(".", " ", location, fixed = TRUE),
    diet     = str_to_title(diet),
    diet     = factor(diet, levels = diet_levels)
  )

landings_diet <- landings_diet %>%
  mutate(
    location = gsub(".", " ", location, fixed = TRUE),
    diet     = str_to_title(diet),
    diet     = factor(diet, levels = diet_levels)
  )

# -----------------------------------------------------------------------------
# 3. Create Plots (Default Axis Lines Restored)
# -----------------------------------------------------------------------------
# Panel A: Reef Observations
p1 <- ggplot(reef_diet, aes(fill = diet, x = location, y = value)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = diet_cols, name = NULL) +
  labs(x = NULL, y = "Proportion of Biomass", title = "Reef Observations") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# Panel B: Fisheries Landings (Y-axis stripped)
p2 <- ggplot(landings_diet, aes(fill = diet, x = location, y = value)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = diet_cols, name = NULL) +
  labs(x = NULL, y = NULL, title = "Fisheries Landings") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.title   = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# -----------------------------------------------------------------------------
# 4. Combine into Two-Panel Figure (Bold A and B Tags)
# -----------------------------------------------------------------------------
combined_plot <- (p1 + p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 16))
  )

# Display plot
combined_plot

ggsave(
  filename = "figures/stacked_bars.png",
  plot     = combined_plot,
  width    = 12,
  height   = 6,
  units    = "in",
  dpi      = 300
)
