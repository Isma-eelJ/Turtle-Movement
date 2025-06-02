# Title: Sea Surfacce Temperature (SST) Preprocessing 
# Author: Isma-eel Jattiem
# Date: May 2025

library(ncdf4)
library(RNetCDF) # when expecting data
library(lubridate)
library(doParallel)
library(progressr)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

sf_use_s2(FALSE)

handlers(global = TRUE)

# Load in Data
sst <- nc_open("Caddy/Caddy_data/Caddy_sst_data.nc")

# From RNetCDF
sst_1 <- open.nc("Caddy/Caddy_data/Caddy_sst_data.nc")

# inspect data
file.inq.nc(sst_1)

# View
print.nc(sst_1)

# No depth
# Moving back to wotking with ssh

# Extract Dimensions
dim_lon <- ncvar_get(sst, "longitude")
dim_lat <- ncvar_get(sst, "latitude")
dim_time <- ncvar_get(sst, "time")

# Extract Variables 
sst_var <- ncvar_get(sst, "analysed_sst", collapse_degen=FALSE)

# Time Conversion
t_units <- ncatt_get(sst, "time", "units")
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + dseconds(dim_time)

# Create Dimensions grid
coords <- as.matrix(expand.grid(dim_lon, dim_lat, date))

# Combine Variables and Dimensions into a single Data Frame
sst_df <- data.frame(cbind(coords, sst_var))

# Renaming columns to stardadized format
names(sst_df) <- c("lon", "lat", "time", "temp")

summary(sst_df)

# Save Files
saveRDS(sst_df, "Caddy_processed_sst.rds")

# Convert once in a separate R script
library(fst)
sst_df <- readRDS("Caddy_processed_sst.rds")
write_fst(sst_df, "Caddy_processed_sst.fst")