# Title: Base Map
# Author: Isma-eel Jattiem
# Date: May 2025

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(ggspatial)

sf_use_s2(FALSE)


# Define Spatial Bounds ---------------------------------------------------
# Will vary for each turtle based on the region they are released and travel 
min_lon <- 15
max_lon <- 32
min_lat <- -36
max_lat <- -28

xlim <- c(min_lon, max_lon)
ylim <- c(min_lat, max_lat)


# Load World Map ----------------------------------------------------------

world <- ne_countries(scale = "large", returnclass = "sf")


# Palettes ----------------------------------------------------------------


# Buffer ------------------------------------------------------------------

# load the countries:
africa <- ne_countries(returnclass = 'sf',
                       continent = "Africa",
                       scale = "large")


# Setting a 5km buffer around the continent
africa_buffer <- africa %>%
  st_buffer(0.05)


# Map ---------------------------------------------------------------------

ggplot() +
# geom_raster() + of desired variable 
  geom_sf(data = africa_buffer, fill = "lightblue", 
          col = "transparent") + # will change to #E6F3FF, simply for visible
  geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.3) +
  coord_sf(ylim = ylim, xlim = xlim, expand = FALSE) +
  scale_x_continuous(labels = function(x) format(abs(x), digits = 3)) +
  scale_y_continuous(labels = function(x) format(abs(x), digits = 3)) +
  annotation_scale(location = "br", width_hint = 0.2, 
                   text_cex = 0.8, bar_cols = c("black", "white")) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
  theme_minimal() +
  theme(
    # Panel styling
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#E6F3FF", color = NA),
    panel.border = element_rect(color = "grey20", fill = NA, size = 1),
    
    # Title styling (positioned top-right like NPP script)
    plot.title = element_text(size = 14, face = "bold", hjust = 0.98, vjust = -0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.98, vjust = -1.2, color = "grey20"),
    
    # Axis styling
    axis.title = element_text(size = 11, color = "grey20"),
    axis.text = element_text(size = 9, color = "grey40"),
    axis.ticks = element_line(color = "grey40"),
    
    # Margins
    plot.margin = margin(20, 20, 20, 20)) + 
  labs(
  title = "Turtle's Name - Variable being viewed",
  subtitle = "Date",
  x = "Longitude (°E)",
  y = "Latitude (°S)")


