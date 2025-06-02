# Title: Processing Current Data Of Caddy Smaller data, quicker trouble shooting
# Author: Isma-eel Jattiem
# Date: May 2025

library(ncdf4)
library(lubridate)
library(dplyr)
library(progressr)
library(foreach)
library(doParallel)
library(fst)

# Enable progress reporting
handlers(global = TRUE)
handlers("progress")

# Load in Data
currents <- nc_open("Caddy/Caddy_data/Caddy_current_data.nc")

# Extract Dimensions
dim_lon <- ncvar_get(currents, "longitude")
dim_lat <- ncvar_get(currents, "latitude")
dim_depth <- ncvar_get(currents, "depth")
dim_time <- ncvar_get(currents, "time")

# Calculate grid sizes for better memory estimation
n_lon <- length(dim_lon)
n_lat <- length(dim_lat)
n_depth <- length(dim_depth)
n_time <- length(dim_time)
total_points <- n_lon * n_lat * n_depth * n_time
cat("Grid size:", n_lon, "x", n_lat, "x", n_depth, "x", n_time, 
    "=", total_points, "points\n")

# Extract Variables - this is more memory efficient than expanding the whole grid
uo <- ncvar_get(currents, "uo", collapse_degen=FALSE)
vo <- ncvar_get(currents, "vo", collapse_degen=FALSE)

# Time Conversion
t_units <- ncatt_get(currents, "time", "units")
t_ustr <- t_units$value

# Print the original time units string to better understand its format
cat("Time units string:", t_ustr, "\n")

# More robust time parsing - this handles common NetCDF time formats
if (grepl("since", t_ustr, ignore.case = TRUE)) {
  # Extract the reference date part after "since"
  date_part <- sub(".*since\\s+([^\\s]+).*", "\\1", t_ustr)
  
  # Try to parse the reference date
  ref_date <- try(as.Date(date_part), silent = TRUE)
  
  if (inherits(ref_date, "try-error") || is.na(ref_date)) {
    # If standard parsing fails, try various common formats
    formats <- c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y")
    for (fmt in formats) {
      ref_date <- try(as.Date(date_part, format = fmt), silent = TRUE)
      if (!inherits(ref_date, "try-error") && !is.na(ref_date)) {
        break
      }
    }
  }
  
  # If we still couldn't parse the date, print a warning
  if (inherits(ref_date, "try-error") || is.na(ref_date)) {
    warning("Could not parse reference date: ", date_part)
    # Default to a placeholder date
    ref_date <- as.Date("1970-01-01")
  }
  
  # Determine the time unit (days, hours, seconds, etc.)
  time_unit <- "seconds" # default
  if (grepl("day", t_ustr, ignore.case = TRUE)) {
    time_unit <- "days"
  } else if (grepl("hour", t_ustr, ignore.case = TRUE)) {
    time_unit <- "hours"
  } else if (grepl("minute", t_ustr, ignore.case = TRUE)) {
    time_unit <- "minutes"
  } else if (grepl("second", t_ustr, ignore.case = TRUE)) {
    time_unit <- "seconds"
  }
  
  # Convert time values based on the unit
  cat("Using reference date:", as.character(ref_date), "with unit:", time_unit, "\n")
  date <- ref_date
  if (time_unit == "days") {
    date <- ref_date + dim_time
  } else if (time_unit == "hours") {
    date <- ref_date + dhours(dim_time)
  } else if (time_unit == "minutes") {
    date <- ref_date + dminutes(dim_time)
  } else {
    date <- ref_date + dseconds(dim_time)
  }
} else {
  warning("Time units string does not contain 'since'. Using time values as is.")
  date <- dim_time
}

# Verify the date conversion worked
cat("First few dates:", head(date), "\n")
# Process data in chunks rather than all at once

# Set up parallel processing
n_cores <- min(6, parallel::detectCores() - 1)
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores for parallel processing\n")

# Function to process data in chunks
process_chunks <- function() {
  # Process by depth levels to save memory
  with_progress({
    p <- progressor(steps = n_depth)
    
    # Explicitly capture all needed variables for parallel processing
    dim_lon_p <- dim_lon
    dim_lat_p <- dim_lat
    dim_depth_p <- dim_depth
    date_p <- date
    uo_p <- uo
    vo_p <- vo
    n_lon_p <- n_lon
    n_lat_p <- n_lat
    n_time_p <- n_time
    n_depth_p <- n_depth
    
    result_list <- foreach(d_idx = 1:n_depth_p, .combine = rbind, 
                           .packages = c("dplyr", "lubridate")) %dopar% {
                             # Update progress
                             p(message = sprintf("Processing depth level %d of %d", d_idx, n_depth_p))
                             
                             # Pre-allocate a data frame for this depth level
                             depth_df <- data.frame(
                               lon = rep(dim_lon_p, each = n_lat_p * n_time_p),
                               lat = rep(rep(dim_lat_p, each = n_time_p), n_lon_p),
                               depth = dim_depth_p[d_idx],
                               time = rep(date_p, n_lon_p * n_lat_p),
                               uo = as.vector(uo_p[, , d_idx, ]),
                               vo = as.vector(vo_p[, , d_idx, ])
                             )
                             
                             return(depth_df)
                           }
    
    return(result_list)
  })
}

# Run the chunked processing
current_df <- process_chunks()

# Save Files
saveRDS(current_df, "Caddy_processed_currents.rds")
# Or as CSV (though RDS preserves data types better)
# write.csv(current_df, "processed_currents.csv", row.names = FALSE)

# Convert once in a separate R script
#current_df <- readRDS("processed_currents_new.rds")
write_fst(current_df, "Caddy_processed_currents.fst")


# Analyse Data
current_analysis <- current_df %>%
  filter(depth == min(depth)) %>% 
  group_by(time) %>% 
  mutate(
    current_speed = sqrt(uo^2 + vo^2),
    current_direction = (90 - (atan2(vo, uo) * 180/pi)) %% 360)

