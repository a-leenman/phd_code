# Script to convert LAS files to DEMs
# 
# Anya Leenman
# 14 Feb 2019
#
# Ammended 9 April 2019 while processing extra data from 10 min intervals (e.g. 10, 20, 40, 50 mins)
# Ammended 16 Nov while processing extra data at 1 min intervals. 
#   Major changes: now using KNN/IDW instead of TIn in grid_terrain. slightly faster and smoother.

# load required packages
library(lidR)
library(rgdal)
#library(raster) # already loaded by lidR

rm(list = ls())

# set run and rep
run <- 2
repe <- 1

# Set wd
setwd(paste0("L:/Experiments/Processing/Run", run, "/Run", run, "_rep", repe, "/_1min_intervals/pointclouds"))

# specify mask of 1 cm from walls (warping caused noise at edges)
xy_min <- 0.010 # Use these lines if specifying mask as filter in readLAS function.
xy_max <- 2.430
z_max <- -0.1 # max elevation on actual fan. Ideally, this removes dangling orange cable in SE corner.
mask <- paste("-keep_xy", xy_min, xy_min, xy_max, xy_max, "-drop_z_above", z_max, sep = " ")

# mask <- readOGR("F:/Experiments/Mask/Mask_5mm.shp") # Use this line to import mask as shapefile

# Set resolution for conversion:
resolution <- 0.001 # in m

# file list
f_list <- list.files(getwd(), pattern = ".las")

for (i in f_list){ # ADJUST INDICES! Should be all :)

  # Import LAS; note threshold set by "filter" option in readLAS function.
  las_f <- readLAS(i, select = "xyzc", filter = mask) 
  las_f$Classification <- 2L # Change classification to "2" (i.e. ground) to make DEM
  
  # get file number (i.e. time)
  t <- paste0(strsplit(i,'')[[1]][11:14], collapse = '') # get time, as character
  
  # if mask is a shapefile, and las was NOT filtered in readLAS: mask the las
  # las_f <- lasclip(las_f, mask)
  
  # convert to DEM
  dem_f <- grid_terrain(las_f, res = resolution, algorithm = knnidw(k = 2L, p = 2), keep_lowest = FALSE) # 30s
  
  #plot(dem_f, main = t)
  
  # Create DEM file name
  dem_name <- paste0("L:/Experiments/Processing/Run", run, "/Run", run, "_rep", repe, "/_1min_intervals/DEM/DEM", t, ".tif")
  
  # write DEM
  writeRaster(dem_f, filename = dem_name, format = "GTiff", overwrite=TRUE)
  print(paste0(t, " is finished!"))
}

