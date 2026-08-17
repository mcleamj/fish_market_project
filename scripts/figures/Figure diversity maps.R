library(ggplot2)
library(sf)
library(dplyr)
library(cowplot)

load("~/fish_market_project/data/loaded_workspace.RData")

# ------------------------------------------------------------------------------
# 1. SETUP BOUNDING BOX & DATA
# ------------------------------------------------------------------------------
xlim_coords <- c(162.882, 163.052)
ylim_coords <- c(5.242, 5.39)

kos_polys_sf     <- st_as_sf(kos_polys)
kos_shoreline_sf <- st_as_sf(kos_shoreline)
kos_reefs_sf     <- st_as_sf(kos_reefs)

map_data_centroid <- kos_polys_sf %>%
  mutate(
    reef_val     = PC1_centroid_reefs[match(geographic, names(PC1_centroid_reefs))],
    landings_val = PC1_centroid_landings[match(geographic, names(PC1_centroid_landings))]
  )
map_data_entropy <- kos_polys_sf %>%
  mutate(
    reef_val     = reef_entropy[match(geographic, names(reef_entropy))],
    landings_val = landings_entropy[match(geographic, names(landings_entropy))]
  )

lims_centroid <- range(c(PC1_centroid_reefs, PC1_centroid_landings), na.rm = TRUE)
lims_entropy  <- range(c(reef_entropy, landings_entropy), na.rm = TRUE)

pc1_diff <- data.frame(
  geographic = names(PC1_centroid_landings),
  diff       = PC1_centroid_reefs - PC1_centroid_landings
) %>% mutate(geographic = factor(geographic, levels = unique(geographic)))

entropy_diff <- data.frame(
  geographic = names(landings_entropy),
  diff       = reef_entropy - landings_entropy
) %>% mutate(geographic = factor(geographic, levels = unique(geographic)))


# ------------------------------------------------------------------------------
# 2. INDIVIDUAL MAP PLOTS
# ------------------------------------------------------------------------------

# PANEL A
p_A <- ggplot(map_data_centroid) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = adjustcolor("lightblue", alpha.f = 0.5)) +
  geom_sf(aes(fill = reef_val), color = "black", linewidth = 0.3) +
  geom_sf(data = kos_reefs_sf, fill = "grey90", color = NA) +
  geom_sf(data = kos_shoreline_sf, fill = "grey", color = "black", linewidth = 0.3) +
  coord_sf(xlim = xlim_coords, ylim = ylim_coords, expand = FALSE, datum = st_crs(4326)) +
  scale_fill_distiller(
    palette = "Blues", direction = 1, limits = lims_centroid, name = "PC1 Centroid",
    guide = guide_colorbar(
      barheight = unit(6.5, "cm"), # Increased height
      barwidth = unit(0.5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.5,
      title.position = "right",
      title.theme = element_text(angle = 90, hjust = 0.5, size = 10, face = "bold")
    )
  ) +
  labs(title = "Reef Observations", y = "Latitude") +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(b = 2)),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.x = element_blank(),
    plot.margin  = margin(1, 1, 1, 1)
  )

# PANEL B
p_B <- ggplot(map_data_centroid) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = adjustcolor("lightblue", alpha.f = 0.5)) +
  geom_sf(aes(fill = landings_val), color = "black", linewidth = 0.3) +
  geom_sf(data = kos_reefs_sf, fill = "grey90", color = NA) +
  geom_sf(data = kos_shoreline_sf, fill = "grey", color = "black", linewidth = 0.3) +
  coord_sf(xlim = xlim_coords, ylim = ylim_coords, expand = FALSE, datum = st_crs(4326)) +
  scale_fill_distiller(
    palette = "Blues", direction = 1, limits = lims_centroid, name = "PC1 Centroid",
    guide = guide_colorbar(
      barheight = unit(6.5, "cm"), # Increased height
      barwidth = unit(0.5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.5,
      title.position = "right",
      title.theme = element_text(angle = 90, hjust = 0.5, size = 10, face = "bold")
    )
  ) +
  labs(title = "Fisheries Landings", y = "Latitude") + # Restored Y label to hold space
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(b = 2)),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "transparent"), # Transparent text locks map size to match A
    plot.margin  = margin(1, 1, 1, 1)
  )

# PANEL D
p_D <- ggplot(map_data_entropy) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = adjustcolor("lightblue", alpha.f = 0.5)) +
  geom_sf(aes(fill = reef_val), color = "black", linewidth = 0.3) +
  geom_sf(data = kos_reefs_sf, fill = "grey90", color = NA) +
  geom_sf(data = kos_shoreline_sf, fill = "grey", color = "black", linewidth = 0.3) +
  coord_sf(xlim = xlim_coords, ylim = ylim_coords, expand = FALSE, datum = st_crs(4326)) +
  scale_fill_distiller(
    palette = "Oranges", direction = 1, limits = lims_entropy, name = "Functional Entropy",
    guide = guide_colorbar(
      barheight = unit(6.5, "cm"), # Increased height
      barwidth = unit(0.5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.5,
      title.position = "right",
      title.theme = element_text(angle = 90, hjust = 0.5, size = 10, face = "bold")
    )
  ) +
  labs(title = "Reef Observations", x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(b = 2)),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin  = margin(1, 1, 1, 1)
  )

# PANEL E
p_E <- ggplot(map_data_entropy) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = adjustcolor("lightblue", alpha.f = 0.5)) +
  geom_sf(aes(fill = landings_val), color = "black", linewidth = 0.3) +
  geom_sf(data = kos_reefs_sf, fill = "grey90", color = NA) +
  geom_sf(data = kos_shoreline_sf, fill = "grey", color = "black", linewidth = 0.3) +
  coord_sf(xlim = xlim_coords, ylim = ylim_coords, expand = FALSE, datum = st_crs(4326)) +
  scale_fill_distiller(
    palette = "Oranges", direction = 1, limits = lims_entropy, name = "Functional Entropy",
    guide = guide_colorbar(
      barheight = unit(6.5, "cm"), # Increased height
      barwidth = unit(0.5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.5,
      title.position = "right",
      title.theme = element_text(angle = 90, hjust = 0.5, size = 10, face = "bold")
    )
  ) +
  labs(title = "Fisheries Landings", x = "Longitude", y = "Latitude") + # Restored Y label to hold space
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(b = 2)),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(color = "transparent"), # Transparent text locks map size to match D
    plot.margin  = margin(1, 1, 1, 1)
  )


# ------------------------------------------------------------------------------
# 3. INDIVIDUAL BARPLOTS
# ------------------------------------------------------------------------------
p_C <- ggplot(pc1_diff, aes(x = geographic, y = diff)) +
  geom_col(fill = "red", color = "black", linewidth = 0.2) +
  scale_y_continuous(position = "right") +
  labs(x = NULL, y = "Difference in PC1 Centroid") +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "grey90"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin        = margin(1, 1, 1, 1)
  )

p_F <- ggplot(entropy_diff, aes(x = geographic, y = diff)) +
  geom_col(fill = "red", color = "black", linewidth = 0.2) +
  scale_y_continuous(position = "right") +
  labs(x = NULL, y = "Difference in Functional Entropy") +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "grey90"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin        = margin(1, 1, 1, 1)
  )


# ------------------------------------------------------------------------------
# 4. EXTRACT LEGENDS & ASSEMBLE GRID
# ------------------------------------------------------------------------------
leg_top    <- get_legend(p_B + theme(legend.position = "right"))
leg_bottom <- get_legend(p_E + theme(legend.position = "right"))

p_A_clean <- p_A + theme(legend.position = "none") + labs(tag = "A")
p_B_clean <- p_B + theme(legend.position = "none") + labs(tag = "B")
p_C_clean <- p_C + labs(tag = "C")

p_D_clean <- p_D + theme(legend.position = "none") + labs(tag = "D")
p_E_clean <- p_E + theme(legend.position = "none") + labs(tag = "E")
p_F_clean <- p_F + labs(tag = "F")

row1 <- plot_grid(
  p_A_clean, p_B_clean, leg_top, p_C_clean,
  ncol = 4,
  rel_widths = c(1, 1, 0.25, 1),
  align = "h",
  axis = "tb"
)

row2 <- plot_grid(
  p_D_clean, p_E_clean, leg_bottom, p_F_clean,
  ncol = 4,
  rel_widths = c(1, 1, 0.25, 1),
  align = "h",
  axis = "tb"
)

final_figure <- plot_grid(
  row1, row2,
  ncol = 1,
  rel_heights = c(1, 1)
)

final_figure

# ------------------------------------------------------------------------------
# 5. SAVE WITH MATCHING MAP CANVAS ASPECT RATIO (COLLAPSES MAP PADDING)
# ------------------------------------------------------------------------------
ggsave(
  filename = "final_6panel_figure.pdf",
  plot = final_figure,
  width = 12.5,
  height = 6.2,
  units = "in",
  dpi = 300
)