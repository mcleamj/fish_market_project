library(sf)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(dplyr)


# ------------------------------------------------------------------------------
# 0. ZONE ABBREVIATIONS MAPPER
# ------------------------------------------------------------------------------
zone_abbrev <- c(
  "North"            = "N",
  "Northwest Corner"  = "NWC",
  "Northwest"         = "NW",
  "West"              = "W",
  "Walung Point"      = "WP", 
  "Northeast Corner"  = "NEC",
  "Northeast"         = "NE",
  "East"              = "E",
  "Southeast"         = "SE",
  "Southeast Corner"  = "SEC"
)

# ------------------------------------------------------------------------------
# 1. LOAD DATA & SET TARGET CRS (UTM 58N - EPSG:32658)
# ------------------------------------------------------------------------------
utm_crs <- 32658 

FSM           <- readRDS("data/gadm36_FSM_0_sp.rds") |> st_as_sf() |> st_transform(utm_crs)
kos_shoreline <- st_read("data/1-01_Shoreline_and_base_layer/kos_shoreline.shp") |> st_transform(utm_crs)
kos_reefs     <- st_read("data/2-01_Reef_base_layer/kos_coral_reefs.shp")     |> st_transform(utm_crs)
kos_polys     <- st_read("data/geo_polygons_layer/kos_geo_polygons.shp")      |> st_transform(utm_crs)

kos_polys <- kos_polys |>
  mutate(abbrev = recode(location, !!!zone_abbrev))

kos_merged    <- st_union(st_geometry(kos_shoreline), st_geometry(kos_reefs))

marina_coords <- read.csv("data/marina coords.csv")
marina_sf     <- st_as_sf(marina_coords, coords = c("long", "lat"), crs = 4326) |> 
  st_transform(utm_crs)
marina_sf$port_label <- tools::toTitleCase(marina_sf$port)

# Inward label positions (Tightened offsets for Okat and Lelu)
marina_sf <- marina_sf |>
  mutate(
    nudge_x = case_when(
      port == "okat"   ~  1000,  # Reduced from 1800 to 1000 meters
      port == "lelu"   ~ -1200,  # Reduced from -1800 to -1000 meters
      port == "walung" ~  1800,
      port == "utwe"   ~   500
    ),
    nudge_y = case_when(
      port == "okat"   ~  -600,  # Reduced from -1200 to -600 meters
      port == "lelu"   ~     0,
      port == "walung" ~  -300,
      port == "utwe"   ~  1000
    )
  )

# EEZ layer
fsm_eez <- st_read(
  "https://geo.vliz.be/geoserver/MarineRegions/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=MarineRegions:eez&outputFormat=application/json&cql_filter=territory1=%27Micronesia%27"
)
# ------------------------------------------------------------------------------
# 2. MAIN MAP (Kosrae)
# ------------------------------------------------------------------------------
bbox_main <- st_bbox(kos_polys)

main_map <- ggplot() +
  # Sector Polygons (Outer layer - Light Gray fill)
  geom_sf(data = kos_polys, aes(fill = "Geographic Sections"), color = "black", linewidth = 0.3) +
  geom_sf_text(data = kos_polys, aes(label = abbrev), size = 4, fontface = "plain") +
  
  # Coral Reefs (Middle layer)
  geom_sf(data = kos_reefs, aes(fill = "Coral Reefs"), color = "black", linewidth = 0.2) +
  
  # Main Island / Shoreline (Inner layer)
  geom_sf(data = kos_shoreline, aes(fill = "Kosrae"), color = "black", linewidth = 0.3) +
  
  # Dummy layer with CRS assigned for "Lagoon & Channels" legend entry
  geom_sf(data = st_sfc(st_polygon(), crs = utm_crs), aes(fill = "Lagoon & Channels"), color = "black") +
  
  # Fill Scale: Set custom order via 'breaks'
  scale_fill_manual(
    name = NULL,
    values = c(
      "Kosrae"              = "#737373",
      "Coral Reefs"         = "#BDBDBD",
      "Lagoon & Channels"   = "#FFFFFF",
      "Geographic Sections" = "#F0F0F0"
    ),
    breaks = c("Kosrae", "Coral Reefs", "Lagoon & Channels", "Geographic Sections")
  ) +
  
  # Marina Points
  geom_sf(data = marina_sf, fill = "white", color = "black", shape = 21, size = 6, stroke = 1.2) +
  
  # Marina Text Labels
  geom_sf_text(
    data = marina_sf, 
    aes(label = port_label), 
    size = 5.2, 
    fontface = "plain",
    nudge_x = marina_sf$nudge_x,
    nudge_y = marina_sf$nudge_y
  ) +
  
  # Expanded Spatial Bounds: Added extra padding (-4000 meters) on the left (xmin)
  coord_sf(
    xlim = c(bbox_main["xmin"] - 4000, bbox_main["xmax"] + 1000),
    ylim = c(bbox_main["ymin"] - 1000, bbox_main["ymax"] + 1000),
    expand = FALSE
  ) +
  
  # North Arrow
  annotation_north_arrow(
    location = "tr", 
    which_north = "true",
    style = north_arrow_orienteering(text_col = "black", fill = c("black", "white"))
  ) +
  
  # Scale Bar
  annotation_scale(
    location = "bl", 
    width_hint = 0.3, 
    style = "bar", 
    height = unit(0.25, "cm"),
    unit_category = "metric"
  ) +
  
  # Theme Settings
  # Theme Settings
  theme_void() +
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position   = c(0.12, 0.17),
    legend.background = element_blank(),
    legend.key        = element_rect(color = "black", linewidth = 0.3),
    legend.text       = element_text(size = 9.5),
    # Keeps the panel border from being clipped off during ggsave
    plot.margin       = margin(3, 3, 3, 3, "pt")
  )

library(rnaturalearth)

# Load global world landmasses
world_land <- ne_countries(scale = "medium", returnclass = "sf")


# ------------------------------------------------------------------------------
# 3. INSET MAP (World landmasses + Enlarged Red Kosrae Box)
# ------------------------------------------------------------------------------
kos_bbox_4326 <- st_bbox(st_transform(kos_merged, 4326))
center_lon <- as.numeric((kos_bbox_4326["xmin"] + kos_bbox_4326["xmax"]) / 2)
center_lat <- as.numeric((kos_bbox_4326["ymin"] + kos_bbox_4326["ymax"]) / 2)

box_size <- 1.2 

enlarged_bbox <- c(
  xmin = center_lon - box_size,
  ymin = center_lat - box_size,
  xmax = center_lon + box_size,
  ymax = center_lat + box_size
)
class(enlarged_bbox) <- "bbox"

kos_enlarged_box <- st_as_sfc(st_set_crs(enlarged_bbox, 4326))

inset_map <- ggplot() +
  # Regional World Landmasses
  geom_sf(data = world_land, fill = "gray50", color = "black", linewidth = 0.2) +
  
  # EEZ Boundary Line
  geom_sf(data = fsm_eez, fill = NA, color = "black", linetype = "dashed", linewidth = 0.4) +
  
  # Enlarged Red Bounding Box around Kosrae
  geom_sf(data = kos_enlarged_box, fill = NA, color = "red", linewidth = 1.0) +
  
  # Text Header
  annotate("text", x = 121, y = 14.5, label = "Federated States of Micronesia EEZ", 
           hjust = -0.3, vjust = 0, size = 2.3, fontface = "plain") +
  
  # Regional Extent
  coord_sf(xlim = c(120, 175), ylim = c(-20, 17), expand = FALSE) +
  
  theme_bw() +
  theme(
    axis.text        = element_text(size = 5.5),
    axis.title       = element_blank(),
    panel.grid.major = element_line(color = "gray85", linewidth = 0.2),
    plot.background  = element_blank(),
    plot.margin      = margin(1, 1, 1, 1)
  )

# ------------------------------------------------------------------------------
# 4. COMPOSE FINAL MAP
# ------------------------------------------------------------------------------
final_plot <- ggdraw(main_map) +
  draw_plot(
    inset_map, 
    x = 0.005,   # Nudges slightly right inside the left border
    y = 0.960,   # Nudges slightly down inside the top border
    width = 0.35, 
    height = 0.28, 
    hjust = 0, 
    vjust = 1
  )

print(final_plot)

# ------------------------------------------------------------------------------
# 5. SAVE HIGH-RESOLUTION OUTPUT
# ------------------------------------------------------------------------------
# ggsave(
#   filename = "figures/Kosrae_Map_R.png",
#   plot     = final_plot,
#   width    = 10,
#   height   = 7.5,
#   units    = "in",
#   dpi      = 300
# )
