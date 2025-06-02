library(ncdf4)
library(RNetCDF)
library(lubridate)
library(doParallel)
library(progressr)
library(tidyverse)
library(fst)
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
npp <- nc_open("Caddy/Caddy_data/Caddy_npp_data.nc")

# Extract Dimensions
dim_lon <- ncvar_get(npp, "longitude")
dim_lat <- ncvar_get(npp, "latitude")
#dim_depth <- ncvar_get(npp, "depth") # No depth
dim_time <- ncvar_get(npp, "time")

# Extract Variables 
npp_var <- ncvar_get(npp, "npp", collapse_degen=FALSE)

# Time Conversion
t_units <- ncatt_get(npp, "time", "units")
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + dseconds(dim_time)
date

# Create Dimensions grid
coords <- as.matrix(expand.grid(dim_lon, dim_lat, date))

# Combine Variables and Dimensions into a single Data Frame
npp_df <- data.frame(cbind(coords, npp_var))

# Renaming columns to stardadized format
names(npp_df) <- c("lon", "lat", "time", "npp")

# Save Files
saveRDS(npp_df, "Caddy_processed_npp.rds")

# Convert once in a separate R script
library(fst)
npp_df <- readRDS("Caddy_processed_npp.rds")
write_fst(npp_df, "Caddy_processed_npp.fst")
