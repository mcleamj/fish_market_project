

library(pals)
library(plot3D)
library(dplyr)


#load("~/fish_market_project/data/loaded_workspace.RData")

# Define custom pink color ramp for Proportion Distinct Species
#pink_ramp <- colorRampPalette(c("#FFF0F5", "#C71585"))

#dark_red_ramp <- colorRampPalette(c("#FFECEC", "#A11D21"))

yellow_ramp <- colorRampPalette(c("#FFFFE5", "#DAA520"))

# 1. OPEN PDF DEVICE
pdf(file = "di_vuln_maps_figure.pdf", 
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
# ROW 1: PROPORTION DISTINCT SPECIES
# ==============================================================================

## PANEL A: REEF OBSERVATIONS (PROPORTION DISTINCT)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(reef_di, 
                      col = yellow_ramp(length(reef_di)),
                      clim = range(c(reef_di, landings_di)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = reef_di,
          colkey = list(width = 1.8, col.clab = "transparent", col.axis = "transparent", col.ticks = "transparent"), 
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = yellow_ramp(length(reef_di)),
          clim = range(c(reef_di, landings_di)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
plot(map_poly, add = TRUE, col = map_poly$color)
title("Reef Observations", font.main = 1, cex.main = 1.85)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
mtext("A", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL B: FISHERIES LANDINGS (PROPORTION DISTINCT)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(landings_di, 
                      col = yellow_ramp(length(landings_di)),
                      clim = range(c(reef_di, landings_di)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = landings_di,
          colkey = list(width = 1.8),
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = yellow_ramp(length(landings_di)),
          clim = range(c(reef_di, landings_di)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
mtext(side = 4, line = 2.6, "Proportion Distinct Species", font = 1, cex = 1.25)
plot(map_poly, add = TRUE, col = map_poly$color)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
title("Fisheries Landings", font.main = 1, cex.main = 1.85)
mtext("B", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL C: BARPLOT (PROPORTION DISTINCT DIFFERENCE)
di_diff <- reef_di - landings_di

par(mar = c(4.8, 2.0, 3.5, 5.0)) 

bp_c <- barplot(di_diff, 
                col = "red", 
                border = NA, 
                axes = FALSE, 
                names.arg = rep("", length(di_diff)), 
                ylim = range(c(di_diff, 0)) * 1.1,
                main = "")

u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "grey90", border = "black", lwd = 1.5)
barplot(di_diff, 
        col = "red", 
        border = NA, 
        add = TRUE, 
        axes = FALSE, 
        names.arg = rep("", length(di_diff)))

# Add X-axis tick marks
axis(1, at = bp_c, labels = FALSE, tcl = -0.3, lwd = 1.5)

# Y-Axis on Right
axis(4, las = 1, cex.axis = 1.25, font = 1, lwd = 1.5)
mtext("Difference in Proportion Distinct", side = 4, line = 3.6, cex = 1.25, font = 1)

# Rotated X-Axis Labels
text(x = bp_c, y = u[3] - (diff(u[3:4]) * 0.04), 
     labels = names(di_diff), 
     srt = 45, adj = 1, xpd = TRUE, cex = 1.20, font = 1)

mtext("C", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


# ==============================================================================
# ROW 2: PROPORTION VULNERABLE SPECIES
# ==============================================================================

## PANEL D: REEF OBSERVATIONS (PROPORTION VULNERABLE)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(reef_vuln, 
                      col = brewer.oranges(length(reef_vuln)),
                      clim = range(c(reef_vuln, landings_vuln)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = reef_vuln,
          colkey = list(width = 1.8, col.clab = "transparent", col.axis = "transparent", col.ticks = "transparent"), 
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = brewer.oranges(length(reef_vuln)),
          clim = range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
plot(map_poly, add = TRUE, col = map_poly$color)
title("Reef Observations", font.main = 1, cex.main = 1.85)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
mtext("D", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL E: FISHERIES LANDINGS (PROPORTION VULNERABLE)
par(mar = c(4.8, 5.0, 3.5, 2.0), cex.axis = 1.25, cex.lab = 1.25)

map_color <- data.frame(
  geographic = geo_coords$geographic,
  color = variablecol(landings_vuln, 
                      col = brewer.oranges(length(landings_vuln)),
                      clim = range(c(reef_vuln, landings_vuln)))
)
map_poly <- kos_polys %>% left_join(map_color, by = "geographic")

scatter2D(geo_coords$Lon, geo_coords$Lat, 
          pch = 19, colvar = landings_vuln,
          colkey = list(width = 1.8),
          asp = 1,
          cex = 0, xlab = "Longitude", ylab = "Latitude",
          xlim = c(162.89, 163.045), ylim = c(5.25, 5.383),
          col = brewer.oranges(length(landings_vuln)),
          clim = range(c(reef_vuln, landings_vuln)))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
     col = adjustcolor("lightblue", alpha = 0.5), lwd = 1.5)
mtext(side = 4, line = 2.6, "Proportion Vulnerable Species", font = 1, cex = 1.25)
plot(map_poly, add = TRUE, col = map_poly$color)
plot(kos_shoreline, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey", add = TRUE)
plot(kos_reefs, xlim = range(reef_data$Lon), ylim = range(reef_data$Lat), col = "grey90", add = TRUE)
title("Fisheries Landings", font.main = 1, cex.main = 1.85)
mtext("E", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)


## PANEL F: BARPLOT (PROPORTION VULNERABLE DIFFERENCE)
vuln_diff <- reef_vuln - landings_vuln

par(mar = c(4.8, 2.0, 3.5, 5.0))

bp_f <- barplot(vuln_diff, 
                col = "red", 
                border = NA, 
                axes = FALSE, 
                names.arg = rep("", length(vuln_diff)), 
                ylim = range(c(vuln_diff, 0)) * 1.1,
                main = "")

u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "grey90", border = "black", lwd = 1.5)
barplot(vuln_diff, 
        col = "red", 
        border = NA, 
        add = TRUE, 
        axes = FALSE, 
        names.arg = rep("", length(vuln_diff)))

# Add X-axis tick marks
axis(1, at = bp_f, labels = FALSE, tcl = -0.3, lwd = 1.5)

# Y-Axis on Right
axis(4, las = 1, cex.axis = 1.25, font = 1, lwd = 1.5)
mtext("Difference in Proportion Vulnerable", side = 4, line = 3.6, cex = 1.25, font = 1)

# Rotated X-Axis Labels
text(x = bp_f, y = u[3] - (diff(u[3:4]) * 0.04), 
     labels = names(vuln_diff), 
     srt = 45, adj = 1, xpd = TRUE, cex = 1.20, font = 1)

mtext("F", side = 3, line = 1, font = 2, adj = -0.1, cex = 1.2)

# 2. CLOSE & SAVE DEVICE
dev.off()

