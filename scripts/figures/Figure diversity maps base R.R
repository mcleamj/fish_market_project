
library(pals)
library(plot3D)
library(dplyr)

#load("~/fish_market_project/data/loaded_workspace.RData")

# 1. OPEN PDF DEVICE
pdf(file = "diversity_map_figure.pdf", 
    width = 11, 
    height = 6.5)

# ------------------------------------------------------------------------------
# SETUP LAYOUT
# ------------------------------------------------------------------------------
layout_matrix <- matrix(c(1, 2, 3,
                          4, 5, 6), 
                        nrow = 2, ncol = 3, byrow = TRUE)

layout(layout_matrix, widths = c(1, 1, 0.88), heights = c(1, 1))


# ==============================================================================
# ROW 1: PC1 CENTROID
# ==============================================================================

## PANEL A: REEF OBSERVATIONS (PC1)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(PC1_centroid_reefs, 
                      col = brewer.blues(length(PC1_centroid_reefs)),
                      clim = range(c(PC1_centroid_reefs, PC1_centroid_landings)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = PC1_centroid_reefs,
          colkey = list(width = 1.8, col.clab = "transparent", col.axis = "transparent", col.ticks = "transparent"), 
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = brewer.blues(length(PC1_centroid_reefs)),
          clim = range(c(PC1_centroid_reefs, PC1_centroid_landings)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
plot(map_poly, add = TRUE, col = map_poly$color)
title("Reef Observations", font.main = 1, cex.main = 1.85)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
mtext("A", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL B: FISHERIES LANDINGS (PC1)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(PC1_centroid_landings, 
                      col = brewer.blues(length(PC1_centroid_landings)),
                      clim = range(c(PC1_centroid_reefs, PC1_centroid_landings)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = PC1_centroid_landings,
          colkey = list(width = 1.8), # Increased width to 1.8
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = brewer.blues(length(PC1_centroid_landings)),
          clim = range(c(PC1_centroid_reefs, PC1_centroid_landings)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
mtext(side = 4, line = 2.6, "PC1 Centroid", font = 1, cex = 1.25)
plot(map_poly, add = TRUE, col = map_poly$color)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
title("Fisheries Landings", font.main = 1, cex.main = 1.85)
mtext("B", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL C: BARPLOT (PC1 DIFFERENCE)
pc1_diff <- PC1_centroid_reefs - PC1_centroid_landings

par(mar = c(4.8, 2.0, 3.5, 5.0)) 

bp_c <- barplot(pc1_diff, 
                col = "red", 
                border = NA, # Removed black border around bars
                axes = FALSE, 
                names.arg = rep("", length(pc1_diff)), 
                ylim = range(c(pc1_diff, 0)) * 1.1,
                main = "")

u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "grey90", border = "black", lwd = 1.5)
barplot(pc1_diff, 
        col = "red", 
        border = NA, # Removed black border around bars
        add = TRUE, 
        axes = FALSE, 
        names.arg = rep("", length(pc1_diff)))

# Add X-axis tick marks
axis(1, at = bp_c, labels = FALSE, tcl = -0.3, lwd = 1.5)

# Y-Axis on Right
axis(4, las = 1, cex.axis = 1.25, font = 1, lwd = 1.5)
mtext("Difference in PC1 Centroid", side = 4, line = 3.6, cex = 1.25, font = 1)

# Rotated X-Axis Labels
text(x = bp_c, y = u[3] - (diff(u[3:4]) * 0.04), 
     labels = names(pc1_diff), 
     srt = 45, adj = 1, xpd = TRUE, cex = 1.20, font = 1)

mtext("C", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


# ==============================================================================
# ROW 2: FUNCTIONAL ENTROPY
# ==============================================================================

## PANEL D: REEF OBSERVATIONS (ENTROPY)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(reef_entropy, 
                      col = brewer.purples(length(reef_entropy)),
                      clim = range(c(reef_entropy, landings_entropy)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = reef_entropy,
          colkey = list(width = 1.8, col.clab = "transparent", col.axis = "transparent", col.ticks = "transparent"), 
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = brewer.purples(length(reef_entropy)),
          clim = range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
plot(map_poly, add = TRUE, col = map_poly$color)
title("Reef Observations", font.main = 1, cex.main = 1.85)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
mtext("D", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL E: FISHERIES LANDINGS (ENTROPY)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(landings_entropy, 
                      col = brewer.purples(length(landings_entropy)),
                      clim = range(c(reef_entropy, landings_entropy)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = landings_entropy,
          colkey = list(width = 1.8), # Increased width to 1.8
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = brewer.purples(length(landings_entropy)),
          clim = range(c(reef_entropy, landings_entropy)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
mtext(side = 4, line = 2.6, "Functional Entropy", font = 1, cex = 1.25)
plot(map_poly, add = TRUE, col = map_poly$color)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
title("Fisheries Landings", font.main = 1, cex.main = 1.85)
mtext("E", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL F: BARPLOT (ENTROPY DIFFERENCE)
entropy_diff <- reef_entropy - landings_entropy

par(mar = c(4.8, 2.0, 3.5, 5.0))

bp_f <- barplot(entropy_diff, 
                col = "red", 
                border = NA, # Removed black border around bars
                axes = FALSE, 
                names.arg = rep("", length(entropy_diff)), 
                ylim = range(c(entropy_diff, 0)) * 1.1,
                main = "")

u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "grey90", border = "black", lwd = 1.5)
barplot(entropy_diff, 
        col = "red", 
        border = NA, # Removed black border around bars
        add = TRUE, 
        axes = FALSE, 
        names.arg = rep("", length(entropy_diff)))

# Add X-axis tick marks
axis(1, at = bp_f, labels = FALSE, tcl = -0.3, lwd = 1.5)

# Y-Axis on Right
axis(4, las = 1, cex.axis = 1.25, font = 1, lwd = 1.5)
mtext("Difference in Functional Entropy", side = 4, line = 3.6, cex = 1.25, font = 1)

# Rotated X-Axis Labels
text(x = bp_f, y = u[3] - (diff(u[3:4]) * 0.04), 
     labels = names(entropy_diff), 
     srt = 45, adj = 1, xpd = TRUE, cex = 1.20, font = 1)

mtext("F", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)

# 2. CLOSE & SAVE DEVICE
dev.off()