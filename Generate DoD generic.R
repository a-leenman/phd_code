# script to calculate a DoD.
# SMoothes DEMs w' 7 x 7 window first.
# Output DoDs and DEMs saved to "summary_data/relevantfolders". Might want to change?
#
# Designed to make 1-min DoDs; at end of loop, t2 DEM gets renamed to t1 DEM so that you're not
# reimporting the same raster twice. However, if you want to use a longer time interval I
# have a diff script for this that I can send too!
#
# Anya Leenman
# 22 Nov 2019
# 20 Jan 2020 Edits to streamline, graph agg/deg volumes
#----------------------------
# Housekeeping
rm(list = ls())

library(raster)
library(dplyr)
library(scales) # contains alpha function
library(RColorBrewer)
library(igraph)

# Specify run and rep
run <- 2
repe <- 1
drive <- "L:" # path of hard drive to write/read to

# Set working directory (WD)
setwd(paste0(drive, "/Experiments/Processing/Run", run, "/Run", run, "_rep", repe, "/_1min_intervals"))

# file list
f_list <-
  list.files(paste0("./DEM"), pattern = ".tif") %>% # DEMS, unsmoothed
  `[`(2:length(.))# adjust so that lists have same timestamp at start and end (mine usually diff lengths)
m_list <-
  list.files(paste0("./DoD_binary"), pattern = ".tif") # Masks

# error trap
if (length(f_list) != length(m_list)) {
  print("oopsies, lists need to have same length")
  break
}

first_index <- 1 # index to start at; default is 1
w_size <- 7 # must be odd number. Window size for DEM smoothing filter. Default is 7

# empty df to save aggradation and degradation volumes, if you need them
agg_deg_vols <- data.frame(t1 = first_index:(length(f_list) - 1), 
                           agg = first_index:(length(f_list) - 1), 
                           deg = first_index:(length(f_list) - 1))

# set up for first file:
f1 <- f_list[first_index]
t1 <- paste0(strsplit(f1, '')[[1]][4:7], collapse = '') # get time, as char
t1_check <- paste0(strsplit(m_list[first_index], '')[[1]][11:14], collapse = '') # get time, as char

# error check
if(t1 != t1_check){
  print("Make sure both lists start at same t-step!")
  break
}
# if you get an error here because either t1 or t1_check is not a character variable made of 4 numbers, 
# adjust the indices in the final square brackets used to define them (e.g. 4:7) so that you're actually
# reading hte t-step from the file name.

# Make t1, before starting loop.
t1_DEM <- raster(paste0("./DEM/", f1)) # make raster
t1_DEM <- focal(t1_DEM, w = matrix(1, w_size, w_size), fun = mean) # smooth it

# write to file
DEM_name <- paste0("./summary_data/DEM_7mm_smoothed/", t1, ".tif")
writeRaster(t1_DEM, filename = DEM_name, format = "GTiff", overwrite=TRUE)


for (i in (first_index + 1):length(f_list)) {
  f2 <- f_list[i] # t2 DEM
  t2 <- paste0(strsplit(f2, '')[[1]][4:7], collapse = '') # get time, as character
  
  m <- m_list[i] # t2 mask
  t2_DEM <- raster(paste0("./DEM/", f2)) # make rasters
  msk <- raster(paste0("./DoD_binary/", m))
  
  # smooth DEM w/moving window
  t2_DEM <- focal(t2_DEM, w = matrix(1, w_size, w_size), fun = mean)
  
  # Difference
  DoD <- t2_DEM - t1_DEM %>%
    `*` (msk) # multiply by mask to extract values only on fan
  threshold <- -0.002 # threshold the DoD, in m. Express as negative value
  DoD[DoD > threshold & DoD < abs(threshold)] <- NA
  # # plot to check if sensible
  # plot(DoD, col = brewer.pal(11, "RdBu"), main = paste0(t2, " - ", t1),
  #      legend.args = list(text = "Elevation change (m)", side = 2))
  
  #-----------------------------
  # remove clumps smaller than 2 cm2:
  # uses code from https://gis.stackexchange.com/questions/130993/remove-clumps-of-pixels-in-r
  rc <- clump(DoD)
  
  minsize <- 2 # smallest allowable clump size, in cm2
  
  #extract IDs of clumps TO BE REMOVED
  clump2cm <- data.frame(freq(rc))
  clump2cm <- clump2cm[clump2cm$count < (minsize * 100),] #remove clump observations with frequency smaller than 200
  clump2cm <- as.vector(clump2cm$value) # record IDs from clumps which met the criteria in previous step
  
  rc[rc %in% clump2cm] <- NA #remove cells with IDs which do not belong to the group of interest
  rc[!is.na(rc)] <- 1 # set all other clumps to value 1 for multiplication
  
  DoD <- DoD * rc # multiply clumps w/ value 1 to remove small clumps from DOD
  # # plot to check if sensible
  # plot(DoD, col = brewer.pal(11, "RdBu"), main = paste0(t2, " - ", t1),
  #      legend.args = list(text = "Elevation change (m)", side = 2))
  
  #-----------------------------------
  
  # save output rasters
  DoD_name <- paste0("./summary_data/DoD/", t2, "-", t1, ".tif")
  DEM_name <- paste0("./summary_data/DEM_7mm_smoothed/", t2, ".tif")
  
  writeRaster(DoD, filename = DoD_name, format = "GTiff", overwrite=TRUE)
  writeRaster(t2_DEM, filename = DEM_name, format = "GTiff", overwrite=TRUE)
  
  # save volumetric data: agg, then deg
  agg_deg_vols$t1[i - first_index] <- t1
  agg_deg_vols$agg[i - first_index] <- sum(DoD[DoD > 0]) / 1000000 # convert to change in m3
  agg_deg_vols$deg[i - first_index] <-sum(DoD[DoD < 0]) / 1000000
  
  
  # progress tracker
  print(paste0("t1 = ", t1, " and t2 = ", t2))
  
  # rename t2 to t1 for next iter.
  f1 <- f2
  t1 <- t2
  t1_DEM <- t2_DEM
}

# write df to file:
write.csv(agg_deg_vols, file = paste0("./t_series/Run", run, "_rep",
                            repe, "_agg_deg_vols.csv"),
row.names = F)

# reload if script closed cause you paused for lunch:
# agg_deg_vols <- read.csv(file = paste0("./t_series/Run",
#                                        run, "_rep", repe, "_agg_deg_vols.csv"))
# colours
cols <- brewer.pal(5, "RdBu")

# plot volumes of agg/deg:
plot(
  1:nrow(agg_deg_vols),
  agg_deg_vols$agg,
  ylim = c(-0.001, 0.001),
  type = "n",
  xlab = "Time since first_index (mins)",
  ylab = "Volume eroded / deposited"
)

lines(agg_deg_vols$agg, col = cols[5])
lines(agg_deg_vols$deg, col = cols[1])

# comparison - maybe you get higher vols of agg and deg at same t-steps? 
# i.e. generally more activity during floods?
# or maybe they are negatively correlated - e.g. less deg when more agg? 
# Only tried this on a subset of my data and looks like the relation is weak, 
# but your floods might create interesting pattern.
# note deg values expressed as negatives so have changed them to pos for this plot.
plot(agg_deg_vols$agg,
     abs(agg_deg_vols$deg),
     pch = 16,
     col = alpha("black", 0.5),
     xlab = "Aggradation volume (m3)",
     ylab = "Degradation volume (m3)")
