# Title: Net Primary Production Daily Maps
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

nc_file <- "Caddy/Caddy_data/harmonized_npp.nc"
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

# Load and analyze NPP data for consistent scaling ------------------------

npp_data <- ncvar_get(nc_data, "npp")

# Calculate global data range for consistent color scaling
cat("Analyzing NPP data range...\n")
npp_range <- range(npp_data, na.rm = TRUE)
cat("NPP data range:", round(npp_range[1], 2), "to", round(npp_range[2], 2), "mg.m^-2.day^-1\n")

# Define quantile-based breaks for better data distribution
npp_quantiles <- quantile(npp_data, probs = seq(0, 1, 0.1), na.rm = TRUE)
cat("NPP quantiles (10% intervals):\n")
print(round(npp_quantiles, 2))

# Create optimized color scale and breaks
breaks <- c(npp_range[1], npp_quantiles[2:10], npp_range[2])
breaks <- unique(round(breaks, 1))  # Remove duplicates and round

# Create labels for the breaks
labels <- character(length(breaks) - 1)
for (i in 1:(length(breaks) - 1)) {
  if (i == 1) {
    labels[i] <- paste0("<", breaks[i + 1])
  } else if (i == length(breaks) - 1) {
    labels[i] <- paste0(">", breaks[i])
  } else {
    labels[i] <- paste0(breaks[i], "–", breaks[i + 1])
  }
}

# Enhanced color palette
palette <- colorRampPalette(c(
  "#000080",  # Dark blue
  "#0000FF",  # Blue
  "#4169E1",  # Royal blue
  "#00BFFF",  # Deep sky blue
  "#00FFFF",  # Cyan
  "#7FFFD4",  # Aquamarine
  "#90EE90",  # Light green
  "#32CD32",  # Lime green
  "#228B22",  # Forest green
  "#006400"   # Dark green
))(length(labels))

cat("Color scale created with", length(labels), "categories\n\n")

# Setup mapping components ------------------------------------------------

output_dir <- "Caddy/Figures/Caddy_npp_pmm_daily_maps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Map extent - use exact data extent
lon_range <- range(lon, na.rm = TRUE)
lat_range <- range(lat, na.rm = TRUE)
xlim <- c(lon_range[1], lon_range[2])
ylim <- c(lat_range[1], lat_range[2])

# Geographic data
world <- ne_countries(returnclass = 'sf', scale = "large")
africa <- ne_countries(returnclass = 'sf', continent = "Africa", scale = "large")
africa_buffer <- africa %>% st_buffer(0.05)

# Enhanced mapping function -----------------------------------------------

create_professional_map <- function(day_index, npp_data, lon, lat, date, world, africa_buffer, 
                                    xlim, ylim, tracking_data, breaks, labels, palette) {
  
  # Prepare NPP data
  npp_df <- expand.grid(lon = lon, lat = lat)
  npp_df$npp <- as.vector(npp_data[, , day_index])
  npp_df <- npp_df[!is.na(npp_df$npp), ]
  
  # Apply consistent binning - ensure all categories are represented
  npp_df$npp_binned <- cut(npp_df$npp, breaks = breaks, labels = labels, 
                           include.lowest = TRUE, right = FALSE)
  
  # Force all levels to be present, even if empty for the current day's data
  npp_df$npp_binned <- factor(npp_df$npp_binned, levels = labels)
  
  # Create a small data frame with one entry for each possible label
  dummy_legend_data <- data.frame(
    lon = xlim[1], # Just a placeholder coordinate, outside visible range or irrelevant
    lat = ylim[1], # Just a placeholder coordinate
    npp_binned = factor(labels, levels = labels) # All labels as factor levels
  )
  
  # Tracking data up to current date
  current_date <- as.Date(date[day_index])
  cumulative_track <- tracking_data[tracking_data$date <= current_date, ]
  
  # Create base map
  p <- ggplot() +
    # Add the dummy data layer first, making it completely transparent
    geom_raster(data = dummy_legend_data, aes(x = lon, y = lat, fill = npp_binned), 
                alpha = 0) + # Set alpha to 0 to make it invisible
    geom_raster(data = npp_df, aes(x = lon, y = lat, fill = npp_binned), 
                interpolate = TRUE, alpha = 0.9) +
    geom_sf(data = africa_buffer, fill = "#E6F3FF", col = "transparent") +
    geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.2) +
    coord_sf(ylim = ylim, xlim = xlim, expand = FALSE)
  
  # Add consistent color scale with all categories
  p <- p + scale_fill_manual(
    name = "NPP\n(mg.m⁻².day⁻¹)",
    values = setNames(palette, labels),
    na.value = "transparent",
    drop = FALSE,
    limits = labels,
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0,
      keywidth = unit(1.2, "cm"),
      keyheight = unit(0.8, "cm"),
      reverse = FALSE,
      ncol = 1,
      override.aes = list(alpha = 1)
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
    title = "Caddy Net Primary Production",
    subtitle = format(current_date, "%d %B %Y"),
    x = "Longitude (°E)",
    y = "Latitude (°S)"
  )
  
  return(p)
}

# Process all daily maps --------------------------------------------------

cat("Creating professional daily maps...\n")
cat("Data range being used:", round(npp_range[1], 2), "to", round(npp_range[2], 2), "mg.m^-2.day^-1\n")
cat("Number of color categories:", length(labels), "\n\n")

# Progress tracking
total_days <- length(date)
start_time <- Sys.time()

for (i in 1:total_days) {
  cat("Processing day", i, "of", total_days, ":", format(date[i], "%Y-%m-%d"), "\n")
  
  # Create professional map
  daily_map <- create_professional_map(i, npp_data, lon, lat, date, world, africa_buffer, 
                                       xlim, ylim, tracking_data, breaks, labels, palette)
  
  # Save with high quality
  filename <- file.path(output_dir, paste0("Caddy_npp_map_", format(date[i], "%Y%m%d"), ".png"))
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
cat("NPP data range:", round(npp_range[1], 2), "to", round(npp_range[2], 2), "mg.m^-2.day^-1\n")
cat("Color categories:", length(labels), "\n")
cat("Map resolution: 14×10 inches at 300 DPI\n")
cat(rep("=", 60), "\n")

output_dir_2 <- "Caddy/Figures/Legends"
if (!dir.exists(output_dir_2)) {
  dir.create(output_dir_2, recursive = TRUE)
}

# Create a sample legend plot for reference
legend_plot <- ggplot(data.frame(x = 1, y = 1:length(labels), 
                                 fill = factor(labels, levels = labels))) +
  geom_tile(aes(x = x, y = y, fill = fill), width = 0.8, height = 0.8) +
  scale_fill_manual(values = palette, name = "NPP (mg.m⁻².day⁻¹)") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(ncol = 1, keywidth = unit(1.5, "cm"), 
                             keyheight = unit(1, "cm")))

ggsave(file.path(output_dir_2, "npp_color_scale_reference.png"), legend_plot, 
       width = 6, height = 8, dpi = 300, bg = "white")

cat("Color scale reference saved as 'color_scale_reference.png'\n")