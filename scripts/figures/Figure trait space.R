#######################################
## TRAIT SPACE FIGURE WITH CENTROIDS ##
#######################################

#load("~/fish_market_project/data/loaded_workspace.RData")

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(vegan)

############################
## PCA AXIS LABELS         ##
############################

PC1_variance <- summary(trait_space)$importance[2,1]*100
PC2_variance <- summary(trait_space)$importance[2,2]*100

pc1_label <- paste0("PC1 (", round(PC1_variance,1), "%)")
pc2_label <- paste0("PC2 (", round(PC2_variance,1), "%)")


############################
## DATA PREPARATION        ##
############################

############################
## DATA PREPARATION        ##
############################

trait_df <- data.frame(
  PC1 = PC1,
  PC2 = PC2,
  Diet = species_info$Diet
)

# Place your updated color mapping here:
diet_cols <- c(
  "Piscivore"   = "#FD3000", # Red
  "Microphage"  = "#34866D", # Dark Green
  "Planktivore" = "#78B214", # Light Green
  "Invertivore" = "#FF9115", # Orange
  "Grazer"      = "#0034F5", # Blue
  "Omnivore"    = "#E9CB2A"  # Yellow
)

############################
## PANEL A - TRAIT SPACE   ##
############################

hull_all <- trait_df[chull(trait_df$PC1, trait_df$PC2), ]

pA <- ggplot(trait_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull_all,
    aes(PC1, PC2),
    fill = "grey",
    alpha = 0.3,
    color = "grey",
    linewidth = 1
  ) +
  geom_point(
    shape = 21,
    color = "black",
    fill = "grey",
    size = 4
  ) +
  labs(
    x = pc1_label,
    y = pc2_label,
    tag = "A"
  ) +
  theme_test() +
  theme(
    plot.tag = element_text(
      face="bold",
      size=18
    ),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )


############################
## PANEL B - VECTORS       ##
############################

vector_df <- as.data.frame(scores(vectors_1_2, display = "vectors"))
vector_df$trait <- rownames(vector_df)
vector_df <- vector_df %>%
  mutate(
    trait = recode(
      trait,
      "relative_eye" = "rel. eye size",
      "relative_maxillary" = "rel. maxillary length",
      "pectoral_fin_size" = "pect. fin size",
      "caudal_throttling" = "caudal throttling",
      "body_lateral_shape" = "body lateral shape",
      "eye_vertical" = "eye vertical pos.",
      "pectoral_fin_vertical" = "pect. fin vertical pos.",
      "oral_gape_trait" = "oral gape",
      "body_elongation" = "body elongation",
      "max_length" = "max. length"
    )
  )

pB <- ggplot(trait_df, aes(PC1, PC2)) +
  geom_point(
    alpha = 0
  ) +
  geom_segment(
    data = vector_df,
    aes(
      x = 0,
      y = 0,
      xend = PC1 * 3,
      yend = PC2 * 3
    ),
    color = "black",
    linewidth = 0.5,
    arrow = arrow(
      length = unit(0.15, "cm"),
      type = "closed"
    )
  ) +
  geom_text(
    data = vector_df,
    aes(
      x = PC1 * 3.25,
      y = PC2 * 3.75,
      label = trait
    ),
    size = 4,
    hjust = 0.5,
    vjust = -0.5
  )  +
  labs(
    x = pc1_label,
    y = pc2_label,
    tag = "B"
  ) +
  theme_test() +
  theme(
    plot.tag = element_text(
      face="bold",
      size=18
    ),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )


############################
## PANEL C - DIET HULLS    ##
############################

diet_hulls <- trait_df %>%
  group_by(Diet) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

pC <- ggplot(trait_df, aes(PC1, PC2)) +
  geom_polygon(
    data = diet_hulls,
    aes(
      fill = Diet,
      color = Diet,
      group = Diet
    ),
    alpha = 0.1,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    aes(color = Diet),
    size = 4
  ) +
  scale_fill_manual(
    values = diet_cols
  ) +
  scale_color_manual(
    values = diet_cols
  ) +
  guides(
    color = guide_legend(nrow = 2)
  ) +
  labs(
    x = pc1_label,
    y = pc2_label,
    tag = "C"
  ) +
  theme_test() +
  theme(
    legend.position = c(0.5, 1.02),
    legend.justification = c(0.5, 0),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank(),
    plot.margin = margin(t = 55, r = 10, b = 10, l = 10, unit = "pt"),
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.75
    )
  ) +
  coord_cartesian(clip = "off")


############################
## PANEL D - CENTROIDS     ##
############################

centroid_df <- data.frame(
  PC1 = c(
    PC1_centroid_reefs,
    PC1_centroid_landings
  ),
  PC2 = c(
    PC2_centroid_reefs,
    PC2_centroid_landings
  ),
  Type = rep(c(
    "Reef Observations",
    "Market Landings"),
    each=length(PC1_centroid_reefs))
)

centroid_df$Type <- factor(
  centroid_df$Type,
  levels = c(
    "Reef Observations",
    "Market Landings"
  )
)

pD <- ggplot(centroid_df, aes(PC1, PC2)) +
  geom_point(
    aes(fill = Type),
    shape = 21,
    color = "black",
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Reef Observations" = "#2B5C8F",
      "Market Landings" = "#D95F02"
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 19,
        size = 5,
        color = c("#2B5C8F", "#D95F02")
      )
    )
  ) +
  coord_cartesian(
    xlim = range(
      1.5*c(
        PC1_centroid_reefs,
        PC1_centroid_landings
      )
    ),
    ylim = range(
      2*c(
        PC2_centroid_reefs,
        PC2_centroid_landings
      )
    )
  ) +
  labs(
    x = pc1_label,
    y = pc2_label,
    tag = "D",
    fill = NULL
  ) +
  theme_test() +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    ),
    axis.text = element_text(
      size = 12
    ),
    axis.title = element_text(
      size = 14
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    ),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.75, "cm"),
    legend.text = element_text(
      size = 10
    )
  )


############################
## COMBINE FIGURE          ##
############################

trait_space_figure <- (pA | pB) /
  (pC | pD)

trait_space_figure

ggsave(
  "figures/trait_space_figure_2.pdf",
  trait_space_figure,
  width = 6.5*2,
  height = 4*2,
  units = "in",
  dpi = 600,
  bg = "white"
)