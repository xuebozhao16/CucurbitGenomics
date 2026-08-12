# ==========================================================
# Title: Extract Environmental Variables for Cucumber Samples
# Author: Xuebo Zhao
# Description:
#   This script extracts bioclimatic, precipitation, solar radiation,
#   and elevation data from WorldClim raster layers for cucumber sample locations.
# ==========================================================

# Load required package
library(raster)

# ==========================================================
# 1. Load WorldClim Raster Layers
# ==========================================================
setwd("/Users/xuebozhao/Documents/CucumberProject/EnvironmentData/")

# Load BIO variables (1–19)
bio_files <- paste0("wc2.1_30s_bio_", 1:19, ".tif")
bio_list <- lapply(bio_files, raster)

# Load Precipitation (monthly 1–12)
prec_files <- paste0("wc2.1_30s_prec_", sprintf("%02d", 1:12), ".tif")
prec_list <- lapply(prec_files, raster)

# Load Solar Radiation (monthly 1–12)
srad_files <- paste0("wc2.1_30s_srad_", sprintf("%02d", 1:12), ".tif")
srad_list <- lapply(srad_files, raster)

# Load Elevation
elev <- raster("wc2.1_30s_elev.tif")

# ==========================================================
# 2. Load Sample Locations
# ==========================================================
setwd("/Users/xuebozhao/Documents/CucumberProject/Environment/")

# Input file should contain columns: lon, lat
location <- read.table("cucumber_sample_locations.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

lon <- location$lon
lat <- location$lat
samples <- data.frame(lon, lat)

# ==========================================================
# 3. Extract Environmental Variables for Each Sample
# ==========================================================
temp.data <- samples

# Extract bioclimatic variables
for (i in 1:19) {
  temp.data[[paste0("bio", i)]] <- extract(bio_list[[i]], samples)
}

# Extract precipitation variables
for (i in 1:12) {
  temp.data[[paste0("prec", sprintf("%02d", i))]] <- extract(prec_list[[i]], samples)
}

# Extract solar radiation variables
for (i in 1:12) {
  temp.data[[paste0("srad", sprintf("%02d", i))]] <- extract(srad_list[[i]], samples)
}

# Extract elevation
temp.data$elevation <- extract(elev, samples)

# ==========================================================
# 4. Save Output
# ==========================================================
write.table(temp.data, "cucumber_environment_variables.txt", sep="\t", row.names=FALSE, quote=FALSE)

cat("Extraction completed successfully. Results saved to cucumber_environment_variables.txt\n")
