# Title: Sea Surface Height (SSH) Preprocessing 
# Author: Isma-eel Jattiem
# Date: May 2025

library(ncdf4)
library(RNetCDF)
library(lubridate)
library(doParallel)
library(progressr)
library(tidyverse)
library(ggplot2)
library(fst)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

sf_use_s2(FALSE)
handlers(global = TRUE)

# Load in Data
ssh <- nc_open("Caddy/Caddy_data/Caddy_ssh_data.nc")

# Extract Dimensions
dim_lon <- ncvar_get(ssh, "longitude")
dim_lat <- ncvar_get(ssh, "latitude")
dim_time <- ncvar_get(ssh, "time")

# Extract Variables 
ssh_var <- ncvar_get(ssh, "zos", collapse_degen=FALSE)

# Improved Time Conversion
t_units_attr <- ncatt_get(ssh, "time", "units")
t_units <- t_units_attr$value

print(paste("Time units:", t_units))

# More robust parsing of time units
# Common formats: "days since YYYY-MM-DD" or "seconds since YYYY-MM-DD HH:MM:SS"
if (grepl("since", t_units)) {
  # Split on "since" to get the reference date
  parts <- strsplit(t_units, " since ")[[1]]
  time_unit <- trimws(parts[1])  # e.g., "days" or "seconds"
  ref_date_str <- trimws(parts[2])  # e.g., "1950-01-01" or "1950-01-01 00:00:00"
  
  # Parse reference date (handle both date-only and datetime formats)
  if (grepl(":", ref_date_str)) {
    # Has time component
    ref_date <- as.POSIXct(ref_date_str, tz = "UTC")
  } else {
    # Date only
    ref_date <- as.POSIXct(paste(ref_date_str, "00:00:00"), tz = "UTC")
  }
  
  print(paste("Reference date:", ref_date))
  print(paste("Time unit:", time_unit))
  print(paste("Sample time values:", paste(head(dim_time), collapse = ", ")))
  
  # Convert based on time unit
  if (grepl("day", time_unit, ignore.case = TRUE)) {
    date <- ref_date + days(dim_time)
  } else if (grepl("hour", time_unit, ignore.case = TRUE)) {
    date <- ref_date + hours(dim_time)
  } else if (grepl("second", time_unit, ignore.case = TRUE)) {
    date <- ref_date + seconds(dim_time)
  } else {
    # Default to days if unclear
    warning("Unknown time unit, assuming days")
    date <- ref_date + days(dim_time)
  }
} else {
  # Fallback to your original method if "since" not found
  warning("Unusual time format, using original parsing method")
  t_ustr <- strsplit(t_units, " ")
  t_dstr <- strsplit(unlist(t_ustr)[3], "-")
  date <- ymd(t_dstr) + dseconds(dim_time)
}

print(paste("Converted dates range:", min(date), "to", max(date)))

# Create Dimensions grid
coords <- as.matrix(expand.grid(dim_lon, dim_lat, date))

# Combine Variables and Dimensions into a single Data Frame
# Fixed: using ssh_var instead of npp_var
ssh_df <- data.frame(cbind(coords, as.vector(ssh_var)))

# Renaming columns to standardized format
names(ssh_df) <- c("lon", "lat", "time", "ssh")

# Convert time back to POSIXct (it gets converted to numeric in the cbind)
ssh_df$time <- as.POSIXct(ssh_df$time, origin = "1970-01-01", tz = "UTC")

summary(ssh_df)

# Save Files
saveRDS(ssh_df, "Caddy_processed_ssh.rds")

# Convert to fst format
write_fst(ssh_df, "Caddy_processed_ssh.fst")