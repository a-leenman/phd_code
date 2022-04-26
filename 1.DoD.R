# Script to make DEMs of Difference
# 
# Anya Leenman
# 19 Feb 2019
#
# Ammended 9 April while processing extra data from 10 min intervals (i.e. 10, 20, 40, 50 mins)
# Ammended 5 June for ease of use.
# Ammended 17 Nov when re-running with 1min intervals.
    # changes: commented out smoothing step.

rm(list = ls())

# load required packages
library(rgdal)
library(raster)
#library(bfastSpatial) # contains areaSieve().

# This version just uses code for areaSieve() function, from Luc Dutrieux's github/R package documentation.
# https://www.rdocumentation.org/packages/bfastSpatial/versions/0.6.2/source

# define areaSieve():
areaSieve <- function(x, thresh=5000, directions=8, verbose=FALSE, keepzeros=FALSE, cores=1, ...)
  
{
  
  
  
  require(igraph)
  
  # convert thresh from area to pixel threshold
  
  # TODO: make this applicable to all projections
  
  thresh <- ceiling(thresh/(res(x)[1]*res(x)[2]))
  
  if(verbose)
    
    cat("Converted threshold to ", thresh, " pixels.\n", sep="")
  
  
  
  # generic sieve function
  
  sieve <- function(inp, ...){
    
    # derive a forest clump raster from unitRaster
    
    clumps <- clump(inp, directions=directions)
    
    
    
    # calculate pixel frequency for each clumpID
    
    clumpFreq <- as.data.frame(freq(clumps))
    
    
    
    # clumpID to be excluded from output raster
    
    excludeID <- clumpFreq$value[which(clumpFreq$count < thresh)]
    
    
    
    # function to assign NA to x wherever a clump with ID %in% excludeID is found
    
    subNA <- function(a, b){
      
      a[b %in% excludeID] <- NA
      
      return(a)
      
    }
    
    
    
    # apply sieve to unitRaster
    
    if(!keepzeros){
      
      y <- overlay(inp, clumps, fun=subNA, ...)
      
    } else {
      
      y <- overlay(inp, clumps, fun=subNA)
      
    }
    
    
    
    return(y)
    
  }
  
  
  
  if(nlayers(x) > 1){
    
    require(doMC)
    
    registerDoMC(cores=cores)
    
    y <- foreach(i = 1:nlayers(x)) %dopar% {
      
      if(keepzeros){
        
        unitRaster <- x[[i]]
        
        unitRaster[!is.na(unitRaster)] <- 1
        
        
        
        # apply sieve on unitRaster
        
        y <- sieve(unitRaster)
        
        
        
        # use sieved unitRaster to mask input raster
        
        y <- mask(x[[i]], y)
        
        
        
      } else {
        
        y <- sieve(x[[i]])
        
      }
      
      return(y)
      
    }
    
    y <- do.call("brick", y)
    
    names(y) <- names(x)
    
    if(hasArg(filename))
      
      writeRaster(x, ...)
    
  } else {
    
    # create a unit raster if keepzeros==TRUE
    
    if(keepzeros){
      
      unitRaster <- x
      
      unitRaster[!is.na(unitRaster)] <- 1
      
      
      
      # apply sieve on unitRaster
      
      y <- sieve(unitRaster)
      
      
      
      # use sieved unitRaster to mask input raster
      
      y <- mask(x, y, ...)
      
      
      
    } else {
      
      y <- sieve(x, ...)
      
    }
    
  }
  
  
  
  return(y)
  
}

# Specify run and rep
run <- 2
rep <- 1
drive <- "L:" # path of hard drive to write/read to
first_index <- 2 # index to start at (usually 5)

# Set working directory (WD)
setwd(paste0(drive, "/Experiments/Processing/Run", run, "/Run", run, "_rep", rep, "/_1min_intervals"))

# file list
f_list <- list.files(paste0("./DEM"), pattern = ".tif")

# import t0 DEM
t0_DEM <- raster(paste0("./DEM/", f_list[1]))

# Set threshold for DoD classification into binary
threshold <- 0.006 # set threshold limit (in m)
#w_size <- 11 # window size for smoothing binary layer; must be odd number. Larger window = slower processing. was 25!

for (i in f_list[first_index:length(f_list)]){ # Index [1] is 0000 data. Adjust index until you get result.
  # import tn DEM
  t <- paste0(strsplit(i,'')[[1]][4:7], collapse = '') # get time, as character
  print(paste0(t, " is begun!"))
  tn_DEM <- raster(paste0(getwd(), "/DEM/", i))
  
  # Difference
  DoD <- tn_DEM - t0_DEM
  DoD_name <- paste0("./DoD/DoD", t, ".tif")
  
  # write output as DoD.tif
  # writeRaster(DoD, filename = DoD_name, format = "GTiff", overwrite=TRUE) # as of May 2021: no longer write! Doesn't need to be saved. Intermediate step.
  
  # reclassify DoD into binary format, based on threshold set above
  rv <- c((DoD@data@min - 1), threshold, 0, threshold, DoD@data@max, 1) # Prepare the reclassification matrix, as a vector
  rm <- matrix(rv, ncol = 3, byrow = TRUE) # reshape into matrix
  DoD_binary <- reclassify(DoD, rm)
  binary_name <- paste0("./DoD_binary/DoD_binary", t, ".tif")
  
  # smoothing filter
  #DoD_binary_2 <- focal(DoD_binary, w = matrix(1, w_size, w_size), fun = modal)
  
  # remove "islands" of noise above drain:
  isize <- ifelse(as.numeric(t) < 400, 0.003 * as.numeric(t), 1)# set island size(isize) as function of time
  DoD_binary[DoD_binary == 0] <- NA # change 0 values to NA
  DoD_binary_3 <- areaSieve(DoD_binary, thresh = isize, directions = 8, verbose = FALSE, keepzeros = TRUE)
  
  # visual test:
  # plot(DoD_binary_3, main = t)
  
  # error trap
  if(DoD_binary_3@data@max != 1){
    print("Error in areaSieve: no data")
    break
  }
  
  # write output to file
  writeRaster(DoD_binary_3, filename = binary_name, format = "GTiff", overwrite=TRUE)
  print(paste0(t, " is finished!"))
}

