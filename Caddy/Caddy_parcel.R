# Title: Caddy Parcel Trajectory Daily Maps
# Author: Isma-eel Jattiem
# Date: June 2025

# Load necessary packages
library(tidyverse)
library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(ggspatial)  
library(viridis)    

sf_use_s2(FALSE)

# Data loading and preprocessing ------------------------------------------

# Load the parcel trajectory data
parcel_data <- read.csv("caddy/Caddy_data/Caddy_parcel_traj.csv")

cat("Original parcel data loaded:\n")
cat("Number of trajectory points:", nrow(parcel_data), "\n")
cat("Trajectory ID:", unique(parcel_data$trajectory), "\n")

# Convert time column to Date format
parcel_data$date <- as.Date(parcel_data$time)

cat("Date range:", as.character(min(parcel_data$date)), "to", as.character(max(parcel_data$date)), "\n")
cat("Number of days:", length(unique(parcel_data$date)), "\n\n")

# Verify the structure
str(parcel_data)
head(parcel_data)
tail(parcel_data)

# Setup mapping components ------------------------------------------------

# Create output directory
output_dir <- "Caddy/Figures/Caddy_parcel_daily_maps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load geographic data
world <- ne_countries(scale = "large", returnclass = "sf")
africa <- ne_countries(returnclass = 'sf', continent = "Africa", scale = "large")
africa_buffer <- africa %>% st_buffer(0.05)

# Define spatial bounds 
data_lon_range <- range(parcel_data$lon, na.rm = TRUE)
data_lat_range <- range(parcel_data$lat, na.rm = TRUE)

min_lon <- 24
max_lon <- 28
min_lat <- -35
max_lat <- -32
xlim <- c(min_lon, max_lon)
ylim <- c(min_lat, max_lat)

cat("Map bounds:\n")
cat("Longitude:", round(xlim, 3), "\n")
cat("Latitude:", round(ylim, 3), "\n\n")

# Prepare daily data structure
unique_dates <- sort(unique(parcel_data$date))

cat("Unique tracking dates:", length(unique_dates), "\n")
cat("Date range:", as.character(min(unique_dates)), "to", as.character(max(unique_dates)), "\n\n")

# Mapping function --------------------------------------------------------

create_professional_parcel_map <- function(current_date, parcel_data, world, africa_buffer, 
                                           xlim, ylim, unique_dates) {
  
  # Get cumulative track up to current date
  cumulative_track <- parcel_data[parcel_data$date <= current_date, ]
  
  # Get current day's points
  current_points <- parcel_data[parcel_data$date == current_date, ]
  
  # Create base map with ocean background
  p <- ggplot() +
    # Ocean background
    geom_sf(data = africa_buffer, fill = "#E6F3FF", col = "transparent") +
    # Coastlines
    geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.3) +
    coord_sf(ylim = ylim, xlim = xlim, expand = FALSE) +
    scale_x_continuous(labels = function(x) format(abs(x), digits = 3)) +
    scale_y_continuous(labels = function(x) format(abs(x), digits = 3))
  
  # Add cumulative tracking path
  if (nrow(cumulative_track) > 1) {
    # Path line
    p <- p + geom_path(data = cumulative_track, 
                       aes(x = lon, y = lat), 
                       color = "#FF6600", 
                       size = 2.0, 
                       alpha = 0.8,
                       linetype = "solid")
    
    # Previous points (smaller, semi-transparent)
    previous_points <- cumulative_track[cumulative_track$date < current_date, ]
    if (nrow(previous_points) > 0) {
      p <- p + geom_point(data = previous_points, 
                          aes(x = lon, y = lat), 
                          color = "#CC4400", 
                          size = 2.5, 
                          alpha = 0.6,
                          shape = 16)
    }
  }
  
  # Add current day's points (if any)
  if (nrow(current_points) > 0) {
    # Current position - larger and more prominent
    p <- p + geom_point(data = current_points, 
                        aes(x = lon, y = lat), 
                        color = "#FF6600", 
                        size = 6, 
                        alpha = 1.0,
                        shape = 16) +
      geom_point(data = current_points, 
                 aes(x = lon, y = lat), 
                 color = "white", 
                 size = 2.5, 
                 alpha = 1.0,
                 shape = 16)
    
    # Add subtle glow effect for current position
    p <- p + geom_point(data = current_points, 
                        aes(x = lon, y = lat), 
                        color = "#FF6600", 
                        size = 8, 
                        alpha = 0.3,
                        shape = 16)
  }
  
  # Add scale bar and north arrow
  p <- p + 
    annotation_scale(location = "br", width_hint = 0.2, 
                     text_cex = 0.8, bar_cols = c("black", "white")) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                           style = north_arrow_fancy_orienteering,
                           height = unit(1.2, "cm"), width = unit(1.2, "cm"))
  
  # Professional theme 
  p <- p + theme_minimal() +
    theme(
      # Panel styling
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#E6F3FF", color = NA),
      panel.border = element_rect(color = "grey20", fill = NA, size = 1),
      
      # Title styling (positioned top-right)
      plot.title = element_text(size = 14, face = "bold", hjust = 0.98, vjust = -0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.98, vjust = -1.2, color = "grey20"),
      
      # Axis styling
      axis.title = element_text(size = 11, color = "grey20"),
      axis.text = element_text(size = 9, color = "grey40"),
      axis.ticks = element_line(color = "grey40"),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Add titles (positioned top-right)
  p <- p + labs(
    title = "Caddy Parcel Trajectory",
    subtitle = format(current_date, "%d %B %Y"),
    x = "Longitude (°E)",
    y = "Latitude (°S)"
  )
  
  # Add tracking statistics as annotation (below north arrow)
  if (nrow(cumulative_track) > 1) {
    # Calculate cumulative distance properly
    cumulative_track_ordered <- cumulative_track[order(cumulative_track$date), ]
    
    # Calculate distances between consecutive points
    distances <- numeric(nrow(cumulative_track_ordered))
    if (nrow(cumulative_track_ordered) > 1) {
      for (j in 2:nrow(cumulative_track_ordered)) {
        # Calculate distance using Haversine formula
        lat1 <- cumulative_track_ordered$lat[j-1] * pi/180
        lon1 <- cumulative_track_ordered$lon[j-1] * pi/180
        lat2 <- cumulative_track_ordered$lat[j] * pi/180
        lon2 <- cumulative_track_ordered$lon[j] * pi/180
        
        dlat <- lat2 - lat1
        dlon <- lon2 - lon1
        a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
        c <- 2 * atan2(sqrt(a), sqrt(1-a))
        distances[j] <- 6371 * c  # Earth's radius in km
      }
    }
    
    total_distance <- sum(distances, na.rm = TRUE)
    current_day <- which(unique_dates == current_date)
    
    stats_text <- paste0(
      "Total distance: ", round(total_distance, 1), " km\n",
      "Day: ", current_day, " of ", length(unique_dates)
    )
    
    p <- p + annotate("text", 
                      x = xlim[1] + 0.01 * diff(xlim), 
                      y = ylim[2] - 0.1 * diff(ylim),
                      label = stats_text, 
                      hjust = 0, vjust = 1,
                      size = 4.2,  
                      color = "grey20",
                      fontface = "bold")  
  }
  
  return(p)
}

# Create all daily maps ---------------------------------------------------

cat("Creating professional daily parcel trajectory maps...\n")
cat("Output directory:", output_dir, "\n\n")

# Progress tracking
total_days <- length(unique_dates)
start_time <- Sys.time()

for (i in 1:total_days) {
  current_date <- unique_dates[i]
  cat("Processing day", i, "of", total_days, ":", format(current_date, "%Y-%m-%d"), "\n")
  
  # Create professional map
  daily_map <- create_professional_parcel_map(current_date, parcel_data, world, africa_buffer, 
                                              xlim, ylim, unique_dates)
  
  # Save with high quality
  filename <- file.path(output_dir, paste0("Caddy_parcel_map_", format(current_date, "%Y%m%d"), ".png"))
  ggsave(filename, daily_map, width = 14, height = 10, dpi = 300, 
         bg = "white", type = "cairo-png")
  
  # Progress updates
  if (i %% 5 == 0 || i == total_days) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    est_total <- elapsed * total_days / i
    remaining <- est_total - elapsed
    cat("Progress:", round(100 * i / total_days, 1), "% |", 
        "Elapsed:", round(elapsed, 1), "min |",
        "Remaining:", round(remaining, 1), "min\n")
  }
}

# Create summary visualization --------------------------------------------

# Add a sequence column for gradient coloring
parcel_data_with_seq <- parcel_data %>%
  arrange(date) %>%
  mutate(sequence = row_number(),
         gradient_color = sequence / max(sequence))

# Add trajectory statistics to summary map
total_points <- nrow(parcel_data)
total_distance_summary <- 0
if (total_points > 1) {
  for (j in 2:total_points) {
    lat1 <- parcel_data$lat[j-1] * pi/180
    lon1 <- parcel_data$lon[j-1] * pi/180
    lat2 <- parcel_data$lat[j] * pi/180
    lon2 <- parcel_data$lon[j] * pi/180
    
    dlat <- lat2 - lat1
    dlon <- lon2 - lon1
    a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    total_distance_summary <- total_distance_summary + (6371 * c)
  }
}

# Creating summary map
summary_map <- ggplot() +
  geom_sf(data = africa_buffer, fill = "#E6F3FF", col = "transparent") +
  geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.3) +
  
  # Add trajectory path with gradient coloring
  geom_path(data = parcel_data_with_seq, aes(x = lon, y = lat, color = gradient_color), 
            size = 1.2, alpha = 0.8) +
  
  # Add intermediate points with gradient (excluding first and last)
  geom_point(data = parcel_data_with_seq[2:(nrow(parcel_data_with_seq)-1), ], 
             aes(x = lon, y = lat, color = gradient_color), 
             size = 2.5, alpha = 0.8) +
  
  # Color scale for gradient from green to red
  scale_color_gradient(low = "green", high = "red", guide = "none") +
  
  # START POINT with glow effect
  # Outer glow
  geom_point(data = parcel_data_with_seq[1, ], 
             aes(x = lon, y = lat), 
             color = "green", size = 8, alpha = 0.3, shape = 16) +
  # Main point
  geom_point(data = parcel_data_with_seq[1, ], 
             aes(x = lon, y = lat), 
             color = "green", size = 6, alpha = 1.0, shape = 16) +
  # Inner white dot
  geom_point(data = parcel_data_with_seq[1, ], 
             aes(x = lon, y = lat), 
             color = "white", size = 2.5, alpha = 1.0, shape = 16) +
  
  # END POINT with glow effect 
  # Outer glow
  geom_point(data = parcel_data_with_seq[nrow(parcel_data_with_seq), ], 
             aes(x = lon, y = lat), 
             color = "red", size = 8, alpha = 0.3, shape = 16) +
  # Main point
  geom_point(data = parcel_data_with_seq[nrow(parcel_data_with_seq), ], 
             aes(x = lon, y = lat), 
             color = "red", size = 6, alpha = 1.0, shape = 16) +
  # Inner white dot
  geom_point(data = parcel_data_with_seq[nrow(parcel_data_with_seq), ], 
             aes(x = lon, y = lat), 
             color = "white", size = 2.5, alpha = 1.0, shape = 16) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  scale_x_continuous(labels = function(x) format(abs(x), digits = 3)) +
  scale_y_continuous(labels = function(x) format(abs(x), digits = 3)) +
  annotation_scale(location = "br", width_hint = 0.2, text_cex = 0.8, 
                   bar_cols = c("black", "white")) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                         style = north_arrow_fancy_orienteering,
                         height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#E6F3FF", color = NA),
    panel.border = element_rect(color = "grey20", fill = NA, size = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey20"),
    axis.title = element_text(size = 11, color = "grey20"),
    axis.text = element_text(size = 9, color = "grey40"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(title = "Caddy Parcel Trajectory",
       subtitle = paste("Tracking period:", format(min(unique_dates), "%d %B"), 
                        "to", format(max(unique_dates), "%d %B %Y")),
       x = "Longitude (°E)", y = "Latitude (°S)")

# Add annotations - start/end info in top right corner with visual dots
start_end_text <- "Start\nEnd"

summary_map <- summary_map + 
  # Add text labels
  annotate("text", 
           x = xlim[2] - 0.04 * diff(xlim), 
           y = ylim[2] - 0.01 * diff(ylim),
           label = start_end_text, 
           hjust = 1, vjust = 1,
           size = 4.2, 
           color = "grey20",
           fontface = "bold") +
  
  # Add START dot (green with glow effect) - positioned to the right of "START"
  # Outer glow for START
  annotate("point", 
           x = xlim[2] - 0.02 * diff(xlim), 
           y = ylim[2] - 0.015 * diff(ylim),
           color = "green", size = 4, alpha = 0.3, shape = 16) +
  # Main dot for START
  annotate("point", 
           x = xlim[2] - 0.02 * diff(xlim), 
           y = ylim[2] - 0.015 * diff(ylim),
           color = "green", size = 3, alpha = 1.0, shape = 16) +
  # Inner white dot for START
  annotate("point", 
           x = xlim[2] - 0.02 * diff(xlim), 
           y = ylim[2] - 0.015 * diff(ylim),
           color = "white", size = 1.3, alpha = 1.0, shape = 16) +
  
  # Outer glow for END
  annotate("point", 
           x = xlim[2] - 0.02 * diff(xlim), 
           y = ylim[2] - 0.045 * diff(ylim),
           color = "red", size = 4, alpha = 0.3, shape = 16) +
  # Main dot for END
  annotate("point", 
           x = xlim[2] - 0.02 * diff(xlim), 
           y = ylim[2] - 0.045 * diff(ylim),
           color = "red", size = 3, alpha = 1.0, shape = 16) +
  # Inner white dot for END
  annotate("point", 
           x = xlim[2] - 0.02 * diff(xlim), 
           y = ylim[2] - 0.045 * diff(ylim),
           color = "white", size = 1.3, alpha = 1.0, shape = 16) +
  # Add statistics in bottom left corner
  annotate("text", 
           x = xlim[1] + 0.01 * diff(xlim), 
           y = ylim[1] + 0.01 * diff(ylim),
           label = paste0("Total distance: ", round(total_distance_summary, 1), " km\n",
                          "Duration: ", length(unique_dates), " days"), 
           hjust = 0, vjust = 0,
           size = 4.2, 
           color = "grey20",
           fontface = "bold")

# Display the summary map
summary_map

output_dir_2 <- "Caddy/Figures"
if (!dir.exists(output_dir_2)) {
  dir.create(output_dir_2, recursive = TRUE)
}

# Save summary map
ggsave(file.path(output_dir_2, "Caddy_complete_parcel_trajectory_summary.png"), summary_map, 
       width = 14, height = 10, dpi = 300, bg = "white", type = "cairo-png")

# Final summary and cleanup -----------------------------------------------

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("\n", rep("=", 60), "\n")
cat("PROCESSING COMPLETE!\n")
cat(rep("=", 60), "\n")
cat("Total daily maps created:", total_days, "\n")
cat("Processing time:", round(total_time, 2), "minutes\n")
cat("Output directory:", output_dir, "\n")
cat("Date range:", format(min(unique_dates), "%d %B %Y"), "to", format(max(unique_dates), "%d %B %Y"), "\n")
cat("Map resolution: 14×10 inches at 300 DPI\n")
cat("Total trajectory distance:", round(total_distance_summary, 1), "km\n")
cat("Trajectory ID:", unique(parcel_data$trajectory), "\n")
cat("Summary map saved as: complete_parcel_trajectory_summary.png\n")
cat(rep("=", 60), "\n")

