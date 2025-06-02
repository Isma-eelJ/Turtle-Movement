# Title: Base Map
# Author: Isma-eel Jattiem
# Date: May 2025

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(viridis)

sf_use_s2(FALSE)

# Define spatial bounds
# Will vary for each turtle based on the region they are released and travel 
min_lon <- 15
max_lon <- 32
min_lat <- -36
max_lat <- -28

xlim <- c(min_lon, max_lon)
ylim <- c(min_lat, max_lat)

# Load world map
world <- ne_countries(scale = "large", returnclass = "sf")


# Buffer ------------------------------------------------------------------

# load the countries:
africa <- ne_countries(returnclass = 'sf',
                       continent = "Africa",
                       scale = "large")


# Setting a 1km buffer around the continent
buffer <- africa %>%
  st_buffer(0.01)


# Map ---------------------------------------------------------------------

ggplot() +
# geom_raster() + of desired variable 
  geom_sf(data = buffer, fill = "red", 
          col = "transparent") + # will change to lightblue, simply for visible
  geom_sf(data = world, fill = "lightgrey", col = "black", linewidth = 0.3) +
  coord_sf(ylim = ylim, xlim = xlim, expand = FALSE) +
  # scale_fill_manual(
  #   name = "variable\n(units)",
  #   values = palette,
  #   na.value = "transparent",
  #   drop = FALSE
  # ) + # Will vary for each variable
  scale_alpha(range = c(0.01, 2.3), guide = "none") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.height = unit(1.4, "cm"), 
    legend.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect("lightblue"),
    axis.title.x = element_text(vjust = -1, size = 11),
    axis.title.y = element_text(vjust = 2, size = 12),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)) +
  labs(
    x = "Longitude", 
    y = "Latitude",
    title = "Turtle's Name - Variable being viewed",
    subtitle = "Date")


