
library(performance)
library(ggeffects)
library(visreg)

#load("~/fish_market_project/data/loaded_workspace.RData")

dev.off()

# DISTINCTIVENESS
reef_di_preds <- data.frame(di = colMeans(brms::posterior_predict(reef_di_model)),
                            source= "reef",
                            geographic = reef_di_model$data$geographic)

market_di_preds <- data.frame(di = colMeans(brms::posterior_predict(landings_di_model)),
                              source="landings",
                              geographic = landings_di_model$data$geographic)

di_diff_preds <- rbind(reef_di_preds, market_di_preds)

di_preds_model <- glm(di ~ source * geographic , data=di_diff_preds, family = Gamma(link="log"))
summary(di_preds_model)
visreg(di_preds_model, scale="response", "source", by="geographic")

# VULNERABILITY

reef_vuln_preds <- data.frame(vuln = colMeans(brms::posterior_predict(reef_vuln_model)),
                              source= "reef",
                              geographic = reef_vuln_model$data$geographic)

market_vuln_preds <- data.frame(vuln = colMeans(brms::posterior_predict(landings_vuln_model)),
                                source="landings",
                                geographic = landings_vuln_model$data$geographic)

vuln_diff_preds <- rbind(reef_vuln_preds, market_vuln_preds)

hist(vuln_diff_preds$vuln)

vuln_preds_model <- glm(vuln ~ source * geographic , data=vuln_diff_preds, family =Gamma(link="log"))
summary(vuln_preds_model)
visreg::visreg(vuln_preds_model, "source", by="geographic")

# LET'S MAKE A RAW DATA PLOT AND A GGPREDICT PLOT FOR SUPP
library(ggeffects)
library(ggplot2)
library(patchwork)
library(ggeffects)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. GENERATE GGEFFECTS PREDICTIONS
# ==============================================================================

# Distinctiveness predictions
di_eff <- ggpredict(di_preds_model, terms = c("geographic", "source"))

# Vulnerability predictions
vuln_eff <- ggpredict(vuln_preds_model, terms = c("geographic", "source"))


# ==============================================================================
# 2. PLOTS WITH "A" AND "B" PANEL MARKERS
# ==============================================================================

source_colors <- c("reef" = "#2B5C8F", "landings" = "#D95F02")
source_labels <- c("reef" = "Reef", "landings" = "Landings")

## PLOT 1: DISTINCTIVENESS (PANEL A)
p_di <- ggplot(di_eff, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group)) +
  geom_pointrange(position = position_dodge(width = 0.5), 
                  size = 0.8, linewidth = 0.9) +
  scale_color_manual(values = source_colors, labels = source_labels, name = "Data Source") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    tag = "A",
    title = "Proportion Distinct",
    x = NULL,
    y = "Predicted Proportion (Distinctiveness)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

## PLOT 2: VULNERABILITY (PANEL B)
p_vuln <- ggplot(vuln_eff, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group)) +
  geom_pointrange(position = position_dodge(width = 0.5), 
                  size = 0.8, linewidth = 0.9) +
  scale_color_manual(values = source_colors, labels = source_labels, name = "Data Source") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    tag = "B",
    title = "Proportion Vulnerable",
    x = NULL,
    y = "Predicted Proportion (Vulnerability)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )


# ==============================================================================
# 3. COMBINE & SAVE TO PNG
# ==============================================================================

combined_plot <- p_di + p_vuln + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")

combined_plot

ggsave("supp_ggeffects_predictions.png", 
       plot = combined_plot, 
       width = 12, 
       height = 5.5, 
       units = "in",
       dpi = 300)