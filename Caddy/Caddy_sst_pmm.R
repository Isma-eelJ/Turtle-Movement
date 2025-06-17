# Title: Sea Surface Temperature Daily Maps
# Author: Isma-eel Jattiem
# Date: June 2025

# Load required libraries
library(ncdf4)
library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(lubridate)
library(ggspatial)  

sf_use_s2(FALSE)

# Load tracking data ------------------------------------------------------

tracking_file <- "Caddy/Caddy_data/Caddy_pmm_path.csv"
tracking_data <- read.csv(tracking_file)
tracking_data$date <- as.Date(tracking_data$date)

cat("Tracking data loaded:\n")
cat("Date range:", as.character(min(tracking_data$date)), "to", as.character(max(tracking_data$date)), "\n")
cat("Number of tracking points:", nrow(tracking_data), "\n\n")

# Load NetCDF data --------------------------------------------------------

nc_file <- "Caddy/Caddy_data/harmonized_sst.nc"
nc_data <- nc_open(nc_file)

# Extract dimensions
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
dim_time <- ncvar_get(nc_data, "time")

# Enhanced time conversion
t_units_attr <- ncatt_get(nc_data, "time", "units")
t_units <- t_units_attr$value

print(paste("Time units:", t_units))

if (grepl("since", t_units)) {
  parts <- strsplit(t_units, " since ")[[1]]
  time_unit <- trimws(parts[1])
  ref_date_str <- trimws(parts[2])
  
  if (grepl(":", ref_date_str)) {
    ref_date <- as.POSIXct(ref_date_str, tz = "UTC")
  } else {
    ref_date <- as.POSIXct(paste(ref_date_str, "00:00:00"), tz = "UTC")
  }
  
  print(paste("Reference date:", ref_date))
  print(paste("Time unit:", time_unit))
  
  if (grepl("day", time_unit, ignore.case = TRUE)) {
    date <- ref_date + days(dim_time)
  } else if (grepl("hour", time_unit, ignore.case = TRUE)) {
    date <- ref_date + hours(dim_time)
  } else if (grepl("second", time_unit, ignore.case = TRUE)) {
    date <- ref_date + seconds(dim_time)
  } else {
    warning("Unknown time unit, assuming days")
    date <- ref_date + days(dim_time)
  }
} else {
  warning("Unusual time format, using fallback method")
  t_ustr <- strsplit(t_units, " ")
  t_dstr <- strsplit(unlist(t_ustr)[3], "-")
  date <- ymd(t_dstr) + dseconds(dim_time)
}

print(paste("Converted dates range:", min(date), "to", max(date)))

# Data info
cat("Data dimensions:\n")
cat("Longitude:", length(lon), "points\n")
cat("Latitude:", length(lat), "points\n")
cat("Time:", length(date), "time steps\n")
cat("Date range:", as.character(min(date)), "to", as.character(max(date)), "\n\n")

# Load and analyze SST data for consistent scaling ------------------------

sst_data <- ncvar_get(nc_data, "analysed_sst")

# Convert from Kelvin to Celsius 
sst_data_celsius <- sst_data - 273.15

# Calculate global data range for consistent color scaling
cat("Analyzing SST data range...\n")
sst_range <- range(sst_data_celsius, na.rm = TRUE)
cat("SST data range:", round(sst_range[1], 2), "to", round(sst_range[2], 2), "°C\n")

# Define quantile-based breaks for better data distribution
sst_quantiles <- quantile(sst_data_celsius, probs = seq(0, 1, 0.1), na.rm = TRUE)
cat("SST quantiles (10% intervals):\n")
print(round(sst_quantiles, 2))

# Create optimized color scale and breaks
breaks <- c(sst_range[1], sst_quantiles[2:10], sst_range[2])
breaks <- unique(round(breaks, 1))  # Remove duplicates and round

# Create labels for the breaks
labels <- character(length(breaks) - 1)
for (i in 1:(length(breaks) - 1)) {
  if (i == 1) {
    labels[i] <- paste0("<", breaks[i + 1], "°C")
  } else if (i == length(breaks) - 1) {
    labels[i] <- paste0(">", breaks[i], "°C")
  } else {
    labels[i] <- paste0(breaks[i], "–", breaks[i + 1], "°C")
  }
}

# Enhanced color palette for temperature (cool to warm)
palette <- colorRampPalette(c(
  "#000080",  # Dark blue (coldest)
  "#0000FF",  # Blue
  "#4169E1",  # Royal blue
  "#00BFFF",  # Deep sky blue
  "#87CEEB",  # Sky blue
  "#FFE4B5",  # Moccasin (neutral)
  "#FFA500",  # Orange
  "#FF4500",  # Orange red
  "#FF0000",  # Red
  "#8B0000"   # Dark red (warmest)
))(length(labels))

cat("Color scale created with", length(labels), "categories\n\n")

# Setup mapping components ------------------------------------------------

output_dir <- "Caddy/Figures/Caddy_sst_pmm_daily_maps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Map extent - use exact data extent
lon_range <- range(lon, na.rm = TRUE)
lat_range <- range(lat, na.rm = TRUE)
xlim <- c(lon_range[1]+ 0.05, lon_range[2]- 0.05)
ylim <- c(lat_range[1]+ 0.05, lat_range[2]- 0.05)

# Geographic data
world <- ne_countries(returnclass = 'sf', scale = "large")
africa <- ne_countries(returnclass = 'sf', continent = "Africa", scale = "large")
africa_buffer <- africa %>% st_buffer(0.05)

# Enhanced mapping function -----------------------------------------------

create_professional_map <- function(day_index, sst_data, lon, lat, date, world, africa_buffer, 
                                    xlim, ylim, tracking_data, breaks, labels, palette) {
  
  # Prepare SST data
  sst_df <- expand.grid(lon = lon, lat = lat)
  sst_df$analysed_sst <- as.vector(sst_data[, , day_index])
  
  # Convert from Kelvin to Celsius
  sst_df$analysed_sst <- sst_df$analysed_sst - 273.15
  sst_df <- sst_df[!is.na(sst_df$analysed_sst), ]
  
  # Apply consistent binning - ensure all categories are represented
  sst_df$sst_binned <- cut(sst_df$analysed_sst, breaks = breaks, labels = labels, 
                           include.lowest = TRUE, right = FALSE)
  
  # Convert to factor with all levels to ensure consistent legend
  sst_df$sst_binned <- factor(sst_df$sst_binned, levels = labels)
  
  # Tracking data up to current date
  current_date <- as.Date(date[day_index])
  cumulative_track <- tracking_data[tracking_data$date <= current_date, ]
  
  # Create base map
  p <- ggplot() +
    geom_raster(data = sst_df, aes(x = lon, y = lat, fill = sst_binned), 
                interpolate = TRUE, alpha = 0.9) +
    geom_sf(data = africa_buffer, fill = "#E6F3FF", col = "transparent") +
    geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.2) +
    coord_sf(ylim = ylim, xlim = xlim, expand = FALSE)
  
  # Add consistent color scale - fixed legend with all categories
  p <- p + scale_fill_manual(
    name = "Temperature\n(°C)",
    values = setNames(palette, labels),  
    na.value = "transparent",
    drop = FALSE,  # Keep all levels even if not present in data
    limits = labels,  # Force all categories to appear
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0,
      keywidth = unit(1.2, "cm"),
      keyheight = unit(0.8, "cm"),
      reverse = FALSE,
      ncol = 1,
      override.aes = list(alpha = 1)  # Ensure legend colors are fully opaque
    )
  ) +
    scale_x_continuous(labels = function(x) format(abs(x), digits = 3)) +
    scale_y_continuous(labels = function(x) format(abs(x), digits = 3))
  
  # Add tracking path
  if (nrow(cumulative_track) > 0) {
    if (nrow(cumulative_track) > 1) {
      p <- p + geom_path(data = cumulative_track, 
                         aes(x = lon, y = lat), 
                         color = "#FF0000", 
                         size = 1.5, 
                         alpha = 0.8,
                         linetype = "solid")
    }
    
    # Previous points
    if (nrow(cumulative_track) > 1) {
      previous_points <- cumulative_track[-nrow(cumulative_track), ]
      p <- p + geom_point(data = previous_points, 
                          aes(x = lon, y = lat), 
                          color = "#8B0000", 
                          size = 2.5, 
                          alpha = 0.7,
                          shape = 16)
    }
    
    # Current position
    current_point <- cumulative_track[cumulative_track$date == current_date, ]
    if (nrow(current_point) > 0) {
      p <- p + geom_point(data = current_point, 
                          aes(x = lon, y = lat), 
                          color = "#FF0000", 
                          size = 5, 
                          alpha = 1.0,
                          shape = 16) +
        geom_point(data = current_point, 
                   aes(x = lon, y = lat), 
                   color = "white", 
                   size = 2, 
                   alpha = 1.0,
                   shape = 16)
    }
  }
  
  # Add scale bar and north arrow - repositioned
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
      # Legend styling
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey20", size = 0.5),
      legend.margin = margin(10, 10, 10, 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      
      # Panel styling
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#E6F3FF", color = NA),
      panel.border = element_rect(color = "grey20", fill = NA, size = 1),
      
      # Title styling (positioned top-right, subtitle moved higher)
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
    title = "Caddy Sea Surface Temperature",
    subtitle = format(current_date, "%d %B %Y"),
    x = "Longitude (°E)",
    y = "Latitude (°s)"
  )
  
  return(p)
}

# Process all daily maps --------------------------------------------------

cat("Creating professional daily SST maps...\n")
cat("Temperature range being used:", round(sst_range[1], 2), "to", round(sst_range[2], 2), "°C\n")
cat("Number of color categories:", length(labels), "\n\n")

# Progress tracking
total_days <- length(date)
start_time <- Sys.time()

for (i in 1:total_days) {
  cat("Processing day", i, "of", total_days, ":", format(date[i], "%Y-%m-%d"), "\n")
  
  # Create professional map
  daily_map <- create_professional_map(i, sst_data, lon, lat, date, world, africa_buffer, 
                                       xlim, ylim, tracking_data, breaks, labels, palette)
  
  # Save with high quality
  filename <- file.path(output_dir, paste0("Caddy_sst_map_", format(date[i], "%Y%m%d"), ".png"))
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

# Cleanup and summary -----------------------------------------------------

nc_close(nc_data)

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("\n" , rep("=", 60), "\n")
cat("PROCESSING COMPLETE!\n")
cat(rep("=", 60), "\n")
cat("Total maps created:", total_days, "\n")
cat("Processing time:", round(total_time, 2), "minutes\n")
cat("Output directory:", output_dir, "\n")
cat("SST data range:", round(sst_range[1], 2), "to", round(sst_range[2], 2), "°C\n")
cat("Color categories:", length(labels), "\n")
cat("Map resolution: 14×10 inches at 300 DPI\n")
cat(rep("=", 60), "\n")

# Create a sample legend plot for reference

output_dir_2 <- "Caddy/Figures/Legends"
if (!dir.exists(output_dir_2)) {
  dir.create(output_dir_2, recursive = TRUE)
}

legend_plot <- ggplot(data.frame(x = 1, y = 1:length(labels), 
                                 fill = factor(labels, levels = labels))) +
  geom_tile(aes(x = x, y = y, fill = fill), width = 0.8, height = 0.8) +
  scale_fill_manual(values = palette, name = "Temperature (°C)") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(ncol = 1, keywidth = unit(1.5, "cm"), 
                             keyheight = unit(1, "cm")))

ggsave(file.path(output_dir_2, "temperature_scale_reference.png"), legend_plot, 
       width = 6, height = 8, dpi = 300, bg = "white")

cat("Temperature scale reference saved as 'temperature_scale_reference.png'\n")

