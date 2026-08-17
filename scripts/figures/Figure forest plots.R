

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

reef_PC1_drivers_model <- readRDS("reef_PC1_drivers_model.rds")
reef_entropy_drivers_model <- readRDS("reef_entropy_drivers_model.rds")

# # ---------------------------
# # 1. Extract Model A (PC1 model)
# # ---------------------------
# post_PC1 <- as.data.frame(as.matrix(reef_PC1_drivers_model))
# 
# post_PC1 <- post_PC1 %>%
#   dplyr::select(
#     b_z_score_2sdtot_grav_pop,
#     b_z_score_2sdfishing_events,
#     b_z_score_2sdreef_area_5km,
#     b_z_score_2sdPC1_all,
#     b_z_score_2sdPC2_all
#   )
# 
# colnames(post_PC1) <- c(
#   "tot_grav_pop",
#   "fishing_events",
#   "reef_area_5km",
#   "PC1_all",
#   "PC2_all"
# )
# 
# post_PC1$model <- "Functional Identify (PC1)"
# 
# # ---------------------------
# # 2. Extract Model B (Entropy model)
# # ---------------------------
# post_ENT <- as.data.frame(as.matrix(reef_entropy_drivers_model))
# 
# post_ENT <- post_ENT %>%
#   dplyr::select(
#     b_z_score_2sdtot_grav_pop,
#     b_z_score_2sdfishing_events,
#     b_z_score_2sdreef_area_5km,
#     b_z_score_2sdPC1_all,
#     b_z_score_2sdPC2_all
#   )
# 
# colnames(post_ENT) <- c(
#   "tot_grav_pop",
#   "fishing_events",
#   "reef_area_5km",
#   "PC1_all",
#   "PC2_all"
# )
# 
# post_ENT$model <- "Functional Entropy"
# 
# # ---------------------------
# # 3. Combine models
# # ---------------------------
# df <- bind_rows(post_PC1, post_ENT)
# 
# # ---------------------------
# # 4. Posterior summaries
# # ---------------------------
# summ <- df %>%
#   pivot_longer(-model, names_to = "term", values_to = "value") %>%
#   group_by(model, term) %>%
#   summarise(
#     median = median(value),
#     lo50 = quantile(value, 0.25),
#     hi50 = quantile(value, 0.75),
#     lo90 = quantile(value, 0.05),
#     hi90 = quantile(value, 0.95),
#     .groups = "drop"
#   )
# 
# # ---------------------------
# # 5. Labels + grouping
# # ---------------------------
# summ <- summ %>%
#   mutate(
#     group = case_when(
#       term %in% c("fishing_events", "tot_grav_pop") ~ "Anthropogenic",
#       TRUE ~ "Environmental"
#     ),
#     term = recode(term,
#                   fishing_events = "Fishing Effort",
#                   tot_grav_pop   = "Market Gravity",
#                   reef_area_5km  = "Reef Area (5 km)",
#                   PC1_all        = "Benthic Composition (PC1)",
#                   PC2_all        = "Benthic Composition (PC2)"
#     ),
#     direction = ifelse(median > 0, "Positive", "Negative")
#   )
# 
# # ---------------------------
# # 6. Order predictors by mean absolute effect across BOTH models
# # ---------------------------
# summ$term <- factor(
#   summ$term,
#   levels = c(
#     "Benthic Composition (PC2)",
#     "Benthic Composition (PC1)",
#     "Reef Area (5 km)",
#     "Fishing Effort",
#     "Market Gravity"
#   )
# )
# 
# summ$model <- factor(
#   summ$model,
#   levels = c("Functional Identify (PC1)", "Functional Entropy")
# )
# 
# # ---------------------------
# # 7. Forest plot (A/B panels)
# # ---------------------------
# ggplot(summ, aes(x = median, y = term)) +
#   
#   # 90% interval
#   geom_errorbarh(aes(xmin = lo90, xmax = hi90),
#                  height = 0,
#                  linewidth = 1,
#                  alpha = 0.6) +
#   
#   # 50% interval (thicker core uncertainty)
#   geom_errorbarh(aes(xmin = lo50, xmax = hi50),
#                  height = 0,
#                  linewidth = 2.5) +
#   
#   # posterior median
#   geom_point(aes(color = direction), size = 6) +
#   
#   # zero reference line
#   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
#   
#   # color scale
#   scale_color_manual(values = c("Positive" = "steelblue",
#                                 "Negative" = "firebrick")) +
#   
#   # A/B comparison panels
#   facet_wrap(~model, ncol = 2, scales = "free_x") +
#   
#   theme_bw(base_size = 15) +
#   theme(
#     legend.position = "none",
#     axis.title = element_blank(),
#     axis.text.y = element_text(size = 12),
#     
#     # 🔥 make facet titles bigger
#     strip.text = element_text(size = 16, face = "bold"),
#     
#     # 🚫 remove ALL grid lines
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )


######----------------------------------------------------

# OPTION 2

######----------------------------------------------------

# ---------------------------
# 1. Extract Model A
# ---------------------------
post_PC1 <- as.data.frame(as.matrix(reef_PC1_drivers_model))

post_PC1 <- post_PC1 %>%
  dplyr::select(
    b_z_score_2sdtot_grav_pop,
    b_z_score_2sdfishing_events,
    b_z_score_2sdreef_area_5km,
    b_z_score_2sdPC1_all,
    b_z_score_2sdPC2_all
  )

colnames(post_PC1) <- c(
  "tot_grav_pop",
  "fishing_events",
  "reef_area_5km",
  "PC1_all",
  "PC2_all"
)

post_PC1$model <- "Functional Identity (PC1)"

# ---------------------------
# 2. Extract Model B
# ---------------------------
post_ENT <- as.data.frame(as.matrix(reef_entropy_drivers_model))

post_ENT <- post_ENT %>%
  dplyr::select(
    b_z_score_2sdtot_grav_pop,
    b_z_score_2sdfishing_events,
    b_z_score_2sdreef_area_5km,
    b_z_score_2sdPC1_all,
    b_z_score_2sdPC2_all
  )

colnames(post_ENT) <- c(
  "tot_grav_pop",
  "fishing_events",
  "reef_area_5km",
  "PC1_all",
  "PC2_all"
)

post_ENT$model <- "Functional Entropy"

# ---------------------------
# 3. Combine models
# ---------------------------
df <- bind_rows(post_PC1, post_ENT)

# ---------------------------
# 4. Posterior summaries
# ---------------------------
summ <- df %>%
  pivot_longer(-model, names_to = "term", values_to = "value") %>%
  group_by(model, term) %>%
  summarise(
    median = median(value),
    lo50 = quantile(value, 0.25),
    hi50 = quantile(value, 0.75),
    lo90 = quantile(value, 0.05),
    hi90 = quantile(value, 0.95),
    .groups = "drop"
  )

# ---------------------------
# 5. Labels + grouping
# ---------------------------
summ <- summ %>%
  mutate(
    term = recode(term,
                  fishing_events = "Fishing Effort",
                  tot_grav_pop   = "Market Gravity",
                  reef_area_5km  = "Reef Area (5 km)",
                  PC1_all        = "Benthic Composition (PC1)",
                  PC2_all        = "Benthic Composition (PC2)"
    ),
    direction = ifelse(median > 0, "Positive", "Negative")
  )

# ---------------------------
# 6. FIXED predictor ordering
# ---------------------------
summ$term <- factor(
  summ$term,
  levels = rev(c(
    "Market Gravity",
    "Fishing Effort",
    "Reef Area (5 km)",
    "Benthic Composition (PC1)",
    "Benthic Composition (PC2)"
  ))
)

# ---------------------------
# 7. Shared x-axis limits
# ---------------------------
#xlims <- range(c(summ$lo90, summ$hi90))

# ---------------------------
# 8. PANEL A
# ---------------------------
p1 <- ggplot(
  filter(summ, model == "Functional Identity (PC1)"),
  aes(x = median, y = term)
) +
  
  geom_errorbarh(aes(xmin = lo90, xmax = hi90),
                 height = 0,
                 linewidth = 1,
                 alpha = 0.6) +
  
  geom_errorbarh(aes(xmin = lo50, xmax = hi50),
                 height = 0,
                 linewidth = 2.5) +
  
  geom_point(aes(color = direction), size = 6) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             linewidth = 0.8) +
  
  scale_color_manual(values = c(
    "Positive" = "steelblue",
    "Negative" = "firebrick"
  )) +
  
  #coord_cartesian(xlim = xlims) +
  
  ggtitle("Functional Identity (PC1)") +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------
# 9. PANEL B
# ---------------------------
p2 <- ggplot(
  filter(summ, model == "Functional Entropy"),
  aes(x = median, y = term)
) +
  
  geom_errorbarh(aes(xmin = lo90, xmax = hi90),
                 height = 0,
                 linewidth = 1,
                 alpha = 0.6) +
  
  geom_errorbarh(aes(xmin = lo50, xmax = hi50),
                 height = 0,
                 linewidth = 2.5) +
  
  geom_point(aes(color = direction), size = 6) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             linewidth = 0.8) +
  
  scale_color_manual(values = c(
    "Positive" = "steelblue",
    "Negative" = "firebrick"
  )) +
  
  #coord_cartesian(xlim = xlims) +
  
  ggtitle("Functional Entropy") +
  
  scale_x_continuous(n.breaks = 4) +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------
# 10. Combine panels
# ---------------------------

# (p1 | p2) +
#   plot_annotation(tag_levels = "A") &
#   theme(
#     plot.tag.position = c(0, 1),  # x, y in npc coordinates
#     plot.tag = element_text(face = "bold", size = 18)
#   )


combined <- (p1 | p2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = "topleft",
    plot.tag = element_text(face = "bold", size = 18)
  )

# Override tag position for panel A only
combined[[1]] <- combined[[1]] + theme(plot.tag.position = c(0.375, 1))
combined[[2]] <- combined[[2]] + theme(plot.tag.position = c(0, 1))

combined

#########################################----------------------------------------------------
#########################################----------------------------------------------------
## NOW THE LANDINGS FOREST PLOTS ##
#########################################----------------------------------------------------
#########################################----------------------------------------------------


landings_PC1_drivers_model <- readRDS("landings_PC1_drivers_model.rds")
landings_entropy_drivers_model <- readRDS("landings_entropy_drivers_model.rds")

# ---------------------------
# 1. Extract Model A
# ---------------------------
post_PC1 <- as.data.frame(as.matrix(landings_PC1_drivers_model))
post_PC1$spearguns <- rep(0)

post_PC1 <- post_PC1 %>%
  dplyr::select(
    b_z_score_2sdPC1_centroid_reefs,
    b_z_score_2sdtot_grav_pop,
    b_z_score_2sdfishing_events,
    b_z_score_2sdwave_energy,
    b_gearShallowBottomFishing,
    spearguns,
    b_z_score_2sdnum_fishers_lines
  )

colnames(post_PC1) <- c(
  "reef_diversity",
  "tot_grav_pop",
  "fishing_events",
  "wave_energy",
  "hook_and_line",
  "spearguns",
  "num_fishers"
)

post_PC1$model <- "Functional Identity (PC1)"

# ---------------------------
# 2. Extract Model B
# ---------------------------
post_ENT <- as.data.frame(as.matrix(landings_entropy_drivers_model))
post_ENT$spearguns <- rep(0)

post_ENT <- post_ENT %>%
  dplyr::select(
    b_z_score_2sdreef_entropyP1,
    b_z_score_2sdtot_grav_pop,
    b_z_score_2sdfishing_events,
    b_z_score_2sdwave_energy,
    b_gearShallowBottomFishing,
    spearguns,
    b_z_score_2sdnum_fishers_lines
  )

colnames(post_ENT) <- c(
  "reef_diversity",
  "tot_grav_pop",
  "fishing_events",
  "wave_energy",
  "hook_and_line",
  "spearguns",
  "num_fishers"
)

post_ENT$model <- "Functional Entropy"

# ---------------------------
# 3. Combine models
# ---------------------------
df <- bind_rows(post_PC1, post_ENT)

# ---------------------------
# 4. Posterior summaries
# ---------------------------
summ <- df %>%
  pivot_longer(-model, names_to = "term", values_to = "value") %>%
  group_by(model, term) %>%
  summarise(
    median = median(value),
    lo50 = quantile(value, 0.25),
    hi50 = quantile(value, 0.75),
    lo90 = quantile(value, 0.05),
    hi90 = quantile(value, 0.95),
    .groups = "drop"
  )

# ---------------------------
# 5. Labels + grouping
# ---------------------------
summ <- summ %>%
  mutate(
    term = recode(term,
                  reef_diversity  = "Reef Diversity Value",
                  tot_grav_pop   = "Market Gravity",
                  fishing_events = "Fishing Effort",
                  wave_energy = "Wave Energy",
                  hook_and_line  = "Hook and Line",
                  spearguns = "Spearguns",
                  num_fishers = "Number of Fishers"
    ),
    direction = ifelse(median > 0, "Positive", "Negative")
  )

# ---------------------------
# 6. FIXED predictor ordering
# ---------------------------
summ$term <- factor(
  summ$term,
  levels = rev(c(
    "Reef Diversity Value",
    "Market Gravity",
    "Fishing Effort",
    "Wave Energy",
    "Hook and Line",
    "Spearguns",
    "Number of Fishers"
  ))
)


# ---------------------------
# 8. PANEL A
# ---------------------------
p1 <- ggplot(
  filter(summ, model == "Functional Identity (PC1)"),
  aes(x = median, y = term)
) +
  
  geom_errorbarh(
    aes(xmin = lo90, xmax = hi90),
    height = 0,
    linewidth = 1,
    alpha = 0.6
  ) +
  
  geom_errorbarh(
    aes(xmin = lo50, xmax = hi50),
    height = 0,
    linewidth = 2.5
  ) +
  
  geom_point(
    data = ~ filter(.x, term != "Spearguns"),
    aes(color = direction),
    size = 6
  ) +
  
  geom_point(
    data = ~ filter(.x, term == "Spearguns"),
    shape = 22,
    size = 6,
    fill = scales::alpha("grey50", 0.5),
    color = "grey50",
    stroke = 1
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "Positive" = "steelblue",
      "Negative" = "firebrick"
    )
  ) +
  
  ggtitle("Functional Identity (PC1)") +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------
# 9. PANEL B
# ---------------------------
p2 <- ggplot(
  filter(summ, model == "Functional Entropy"),
  aes(x = median, y = term)
) +
  
  geom_errorbarh(
    aes(xmin = lo90, xmax = hi90),
    height = 0,
    linewidth = 1,
    alpha = 0.6
  ) +
  
  geom_errorbarh(
    aes(xmin = lo50, xmax = hi50),
    height = 0,
    linewidth = 2.5
  ) +
  
  geom_point(
    data = ~ filter(.x, term != "Spearguns"),
    aes(color = direction),
    size = 6
  ) +
  
  geom_point(
    data = ~ filter(.x, term == "Spearguns"),
    shape = 22,
    size = 6,
    fill = scales::alpha("grey50", 0.5),
    color = "grey50",
    stroke = 1
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "Positive" = "steelblue",
      "Negative" = "firebrick"
    )
  ) +
  
  ggtitle("Functional Entropy") +
  
  scale_x_continuous(n.breaks = 4) +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14),
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------
# 10. Combine panels
# ---------------------------

combined <- (p1 | p2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = "topleft",
    plot.tag = element_text(face = "bold", size = 18)
  )

# Override tag position for panel A only
combined[[1]] <- combined[[1]] + theme(plot.tag.position = c(0.3, 1))
combined[[2]] <- combined[[2]] + theme(plot.tag.position = c(0, 1))

combined

# ggsave(
#   "Figure7.png",
#   plot = (p1 | p2) + plot_annotation(tag_levels = "A"),
#   #width = 7.0,
#   #height = 4.0,
#   #units = "in",
#   dpi = 600,
#   type = "cairo-png"
# )
# 
# ggsave(
#   "Figure7.pdf",
#   plot = (p1 | p2) + plot_annotation(tag_levels = "A"),
#   #width = 7.0,
#   #height = 4.0,
#   #units = "in",
#   device = cairo_pdf
# )
