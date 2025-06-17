# Title: Current Data Processing and Daily Mapping
# Author: Isma-eel Jattiem
# Date: June 2025

# Libraries
library(ncdf4)
library(fst)
library(metR)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(viridis)
library(oce)
library(ocedata)
library(ggspatial)

sf_use_s2(FALSE)

# Load tracking data ------------------------------------------------------

tracking_file <- "Caddy/Caddy_data/Caddy_pmm_path.csv"  

# Read tracking data
tracking_data <- read.csv(tracking_file)
tracking_data$date <- as.Date(tracking_data$date)

# Print tracking data info
cat("Tracking data loaded:\n")
cat("Date range:", as.character(min(tracking_data$date)), "to", as.character(max(tracking_data$date)), "\n")
cat("Number of tracking points:", nrow(tracking_data), "\n\n")

# Load NetCDF current data -----------------------------------------------

# Define file path to your NetCDF file
nc_file <- "Caddy/Caddy_data/Regridded/harmonized_current.nc"  

# Open the NetCDF file
nc_data <- nc_open(nc_file)

# Print NetCDF structure to understand variable names
print(nc_data)

# Extract dimension information
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
dim_time <- ncvar_get(nc_data, "time")

# Time conversion 
t_units_attr <- ncatt_get(nc_data, "time", "units")
t_units <- t_units_attr$value

print(paste("Time units:", t_units))

# Parse time units
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
  warning("Unusual time format, using fallback parsing")
  t_ustr <- strsplit(t_units, " ")
  t_dstr <- strsplit(unlist(t_ustr)[3], "-")
  date <- ymd(t_dstr) + dseconds(dim_time)
}

print(paste("Converted dates range:", min(date), "to", max(date)))

# Get current data
# Check available variables first
print("Available variables:")
print(names(nc_data$var))

# Check if depth dimension exists
depth_exists <- "depth" %in% names(nc_data$dim) || "lev" %in% names(nc_data$dim) || "level" %in% names(nc_data$dim)

if (depth_exists) {
  # Try common depth dimension names
  if ("depth" %in% names(nc_data$dim)) {
    depth <- ncvar_get(nc_data, "depth")
    depth_name <- "depth"
  } else if ("lev" %in% names(nc_data$dim)) {
    depth <- ncvar_get(nc_data, "lev")
    depth_name <- "lev"
  } else if ("level" %in% names(nc_data$dim)) {
    depth <- ncvar_get(nc_data, "level")
    depth_name <- "level"
  }
  
  cat("Depth levels found:", length(depth), "levels\n")
  cat("Depth values:", head(depth), "...\n")
  
  # Extract surface level (typically first level, index 1)
  surface_index <- 1
  cat("Using surface level at index", surface_index, "(depth =", depth[surface_index], ")\n")
  
  # Extract U and V components for surface level only
  uo_data <- ncvar_get(nc_data, "uo", start = c(1, 1, surface_index, 1), 
                       count = c(-1, -1, 1, -1))  # Get all lon, lat, surface depth, all time
  vo_data <- ncvar_get(nc_data, "vo", start = c(1, 1, surface_index, 1), 
                       count = c(-1, -1, 1, -1))  # Get all lon, lat, surface depth, all time
  
  # Remove the depth dimension (which is now size 1)
  uo_data <- drop(uo_data)
  vo_data <- drop(vo_data)
  
} else {
  # No depth dimension - extract normally
  uo_data <- ncvar_get(nc_data, "uo")  # East-west current component
  vo_data <- ncvar_get(nc_data, "vo")  # North-south current component
}

# Check dimensions of the current data after processing
cat("U current data dimensions:", dim(uo_data), "\n")
cat("V current data dimensions:", dim(vo_data), "\n")

# Determine if data is 2D (lon, lat) or 3D (lon, lat, time)
is_3d <- length(dim(uo_data)) == 3

# Print data dimensions
cat("Data dimensions:\n")
cat("Longitude:", length(lon), "points\n")
cat("Latitude:", length(lat), "points\n")
cat("Time:", length(date), "time steps\n")
cat("Date range:", as.character(min(date)), "to", as.character(max(date)), "\n\n")

# Base map setup ----------------------------------------------------------

# Create output directory
output_dir <- "Caddy/Figures/Caddy_current_pmm_daily_maps"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# Get data extent for map zooming (like in extraction script)
lon_range <- range(lon, na.rm = TRUE)
lat_range <- range(lat, na.rm = TRUE)
xlim <- c(lon_range[1], lon_range[2])  # Add small buffer
ylim <- c(lat_range[1], lat_range[2])  # Add small buffer

# Load world map
world <- ne_countries(scale = "large", returnclass = "sf")

# Buffer
africa <- ne_countries(returnclass = 'sf',
                       continent = "Africa",
                       scale = "large")

africa_buffer <- africa %>%
  st_buffer(0.05)


# Function to process current data for a single day -----------------------

process_current_day <- function(day_index, uo_data, vo_data, lon, lat, is_3d) {
  
  # Extract current data for the specific day
  if (is_3d) {
    # 3D data: extract specific time slice
    uo_day <- uo_data[, , day_index]
    vo_day <- vo_data[, , day_index]
  } else {
    # 2D data: use the same data for all days
    uo_day <- uo_data
    vo_day <- vo_data
  }
  
  # Create data frame
  current_df <- expand.grid(lon = lon, lat = lat)
  current_df$uo <- as.vector(uo_day)
  current_df$vo <- as.vector(vo_day)
  
  # Remove NA values and filter to reasonable bounds
  current_df <- current_df %>% 
    filter(!is.na(uo), !is.na(vo)) %>%
    # Filter to data extent (using xlim/ylim bounds)
    filter(lon >= xlim[1], lon <= xlim[2],
           lat >= ylim[1], lat <= ylim[2])
  
  if(nrow(current_df) == 0) {
    return(NULL)
  }
  
  # Convert to spatial points for gridding
  current_sf <- current_df %>% 
    st_as_sf(coords = c("lon", "lat")) %>%
    st_set_crs(4326)
  
  # Create grid (coarser grid for broader patterns)
  current_grid <- current_sf %>% 
    st_make_grid(n = c(40, 30)) %>%
    st_sf()
  
  # Grid the data
  current_gridded <- current_grid %>% 
    mutate(
      id = 1:n(), 
      contained = lapply(st_contains(st_sf(geometry), current_sf), identity),
      obs = sapply(contained, length),
      u = sapply(contained, function(x) {median(current_sf[x,]$uo, na.rm = TRUE)}),
      v = sapply(contained, function(x) {median(current_sf[x,]$vo, na.rm = TRUE)})
    ) %>%
    select(obs, u, v) %>% 
    filter(!is.na(u), !is.na(v), obs > 0)
  
  if(nrow(current_gridded) == 0) {
    return(NULL)
  }
  
  # Extract centroid coordinates
  coordinates <- current_gridded %>% 
    st_centroid() %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    rename(x = X, y = Y)
  
  # Remove geometry
  st_geometry(current_gridded) <- NULL
  
  # Combine coordinates with current data
  current_combined <- coordinates %>% 
    bind_cols(current_gridded)
  
  # Interpolate U component
  u_interp <- tryCatch({
    interpBarnes(
      x = current_combined$x, 
      y = current_combined$y, 
      z = current_combined$u,
      xg = seq(xlim[1], xlim[2], length.out = 50),
      yg = seq(ylim[1], ylim[2], length.out = 40)
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if(is.null(u_interp)) {
    return(NULL)
  }
  
  # Interpolate V component
  v_interp <- tryCatch({
    interpBarnes(
      x = current_combined$x, 
      y = current_combined$y, 
      z = current_combined$v,
      xg = seq(xlim[1], xlim[2], length.out = 50),
      yg = seq(ylim[1], ylim[2], length.out = 40)
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if(is.null(v_interp)) {
    return(NULL)
  }
  
  # Get dimensions
  dimension <- data.frame(lon = u_interp$xg, u_interp$zg) %>% dim()
  
  # Create U component table
  u_tb <- data.frame(lon = u_interp$xg, u_interp$zg) %>% 
    gather(key = "lata", value = "u", 2:dimension[2]) %>% 
    mutate(lat = rep(u_interp$yg, each = dimension[1])) %>% 
    select(lon, lat, u) %>% 
    as_tibble()
  
  # Create V component table
  v_tb <- data.frame(lon = v_interp$xg, v_interp$zg) %>% 
    gather(key = "lata", value = "v", 2:dimension[2]) %>% 
    mutate(lat = rep(v_interp$yg, each = dimension[1])) %>% 
    select(lon, lat, v) %>% 
    as_tibble()
  
  # Combine U and V components and calculate velocity
  uv_combined <- u_tb %>% 
    bind_cols(v_tb %>% select(v)) %>% 
    mutate(vel = sqrt(u^2 + v^2)) %>%
    filter(!is.na(u), !is.na(v), !is.infinite(u), !is.infinite(v))
  
  return(uv_combined)
}


# Calculate global min/max current speeds for fixed legend scale ------------------------------------
cat("Calculating global min/max current speeds for fixed legend...\n")
all_vel_values <- c()
for (i in 1:length(date)) {
  daily_data <- process_current_day(i, uo_data, vo_data, lon, lat, is_3d)
  if (!is.null(daily_data) && nrow(daily_data) > 0) {
    all_vel_values <- c(all_vel_values, daily_data$vel)
  }
}
global_min_vel <- min(all_vel_values, na.rm = TRUE)
global_max_vel <- max(all_vel_values, na.rm = TRUE)
cat(paste("Global min velocity:", round(global_min_vel, 3), "\n"))
cat(paste("Global max velocity:", round(global_max_vel, 3), "\n\n"))


# Function to create daily current map ------------------------------------

create_daily_current_map <- function(day_index, uo_data, vo_data, lon, lat, date, world, africa_buffer, xlim, ylim, tracking_data, is_3d, global_min_vel, global_max_vel) {
  
  # Process current data for this day
  uv_combined <- process_current_day(day_index, uo_data, vo_data, lon, lat, is_3d)
  
  if(is.null(uv_combined) || nrow(uv_combined) == 0) {
    # Create empty plot if no data
    p <- ggplot() +
      geom_sf(data = africa_buffer, fill = "#E6F3FF", col = "transparent") +
      geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.2) +
      coord_sf(ylim = ylim, xlim = xlim, expand = FALSE) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "lightblue")
      ) +
      labs(
        title = "Ocean Current Flow Patterns - No Data Available",
        subtitle = paste(format(date, "%Y-%m-%d")),
        x = "Longitude",
        y = "Latitude"
      )
    return(p)
  }
  
  # Get cumulative tracking data up to current date
  current_date <- as.Date(date)
  cumulative_track <- tracking_data[tracking_data$date <= current_date, ]
  
  # Create current streamline plot
  p <- ggplot() +
    geom_streamline(
      data = uv_combined, 
      aes(x = lon, y = lat, dx = u, dy = v, 
          color = sqrt(after_stat(dx^2) + after_stat(dy^2)), 
          alpha = after_stat(step)),
      L = 2,
      res = 2,
      n = 15,
      lineend = "butt",
      linewidth = 0.5) +
    geom_sf(data = africa_buffer, fill = "#E6F3FF", col = "transparent") +
    geom_sf(data = world, fill = "grey85", col = "grey30", linewidth = 0.2) +
    coord_sf(ylim = ylim, xlim = xlim, expand = FALSE) +
    scale_color_viridis_c(
      name = "Current\nSpeed\n(m.s⁻¹)",
      option = "plasma",
      trans = "sqrt",
      limits = c(global_min_vel, global_max_vel)) + # Fixed limits
    scale_alpha(range = c(0.01, 2.3), guide = "none") +
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
      plot.background = element_rect(colour = "white"),
      
      # Axis styling
      axis.title = element_text(size = 11, color = "grey20"),
      axis.text = element_text(size = 9, color = "grey40"),
      axis.ticks = element_line(color = "grey40"),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Add titles (positioned top-right)
  p <- p + labs(
    title = "Caddy Ocean Current",
    subtitle = format(current_date, "%d %B %Y"),
    x = "Longitude (°E)",
    y = "Latitude (°S)"
  )
  
  return(p)
}

# Main processing loop ----------------------------------------------------

cat("Processing daily current maps...\n")

# Progress tracking
total_days <- length(date)
start_time <- Sys.time()

for (i in 1:total_days) {
  cat("Processing day", i, "of", total_days, ":", format(date[i], "%Y-%m-%d"), "\n")
  
  # Create map for current day
  daily_map <- create_daily_current_map(i, uo_data, vo_data, lon, lat, date[i], world, africa_buffer, xlim, ylim, tracking_data, is_3d, global_min_vel, global_max_vel)
  
  # Save the map
  filename <- file.path(output_dir, paste0("Caddy_current_map_", format(date[i], "%Y%m%d"), ".png"))
  ggsave(filename, daily_map, width = 14, height = 10, dpi = 300, bg = "white", type = "cairo-png")
  
  # Display progress every 5 days
  if (i %% 5 == 0 || i == total_days) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    est_total <- elapsed * total_days / i
    remaining <- est_total - elapsed
    cat("Progress:", round(100 * i / total_days, 1), "% |", 
        "Elapsed:", round(elapsed, 1), "min |",
        "Remaining:", round(remaining, 1), "min\n")
  }
}

# Close the NetCDF file
nc_close(nc_data)

cat("\nAll current maps have been saved to the '", output_dir, "' directory.\n")

