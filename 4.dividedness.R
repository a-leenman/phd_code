# Script to calculate flow "dividedness" (number of channels) on fan experiments
#
# Anya Leenman
# 18 March 2019
#
#######################################################################################################
# Housekeeping
### Housekeeping
rm(list = ls())

run <- 2 # run number
rep <- 1 # Experimental repeat number
drive <- "L:" # drive to write to
#first_index_f <- 1 #first index of binary map
#first_index_a <- 1 # first index of fan mask. Adjust so both lists start at same t-value


# Set working directory
setwd(paste0(drive, "/Experiments/Processing/Run", run, "/Run", run, "_rep", rep, "/_1min_intervals"))

# load packages
library(raster)
#library(bfastSpatial) # contains areaSieve

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

######################################################################################################
# organizational stuff
# load data
f_list <- list.files(paste0("./wet_binary"), pattern = ".tif")
a_list <- read.csv(paste0("./t_series/Run", run, "_rep", rep, ".csv"))
#f_list <- f_list[first_index_f:length(f_list)]
#a_list <- a_list[first_index_a:length(a_list[[1]]), ] 

# empty vec for outputs
# dividedness <- matrix(ncol = 3, nrow = length(f_list))
#######################################################################################################

for (i in 1:length(f_list)){
  # load binary channel map
  f <- raster(paste0("./wet_binary/", f_list[i]))
  t <- paste0(strsplit(f_list[i],'')[[1]][11:14], collapse = '') # get time, as character

  # smooth
  f[is.na(f)] <- 0 # NA to 0, so modal function runs correctly
  w_size <- 21 # window size, in number of cells (each cell = 1 mm accross). ODD number
  f2 <- focal(f, w = matrix(1, w_size, w_size), fun = modal) # matrix must be specified in this format
  # or function will edit raster values themselves...
  
  # remove islands
  f2 <- areaSieve(f2, thresh = 0.001, directions = 8, verbose = FALSE)
  f2[is.na(f2)] <- 0
  # # estimate fan radius:
  # radius <- sqrt((a_list$fan.area[[i]] * 4)/pi)
  # arc <- c(0.5, 0.8) # proportions downfan for which cross-sections req'd
  # channels <- arc # vector with dimensions that match input cases
  # 
  # for (r in 1:2){
  #   radius2 <- radius * arc[r] # to get line to intersect with lower fan.
  # # get circle at given radius
  #   x_line <- seq(0, radius2, (radius2/1001))
  #   y_line <- sqrt(radius2 ^ 2 - x_line ^ 2)
  #   cds <- cbind(x_line, y_line)
  #   l <- spLines(cds)
  # 
  #   plot(f2, main = t)
  #   lines(l)
  # 
  #   # calculate number of intersections with line
  #   prof <- extract(f2, l)
  #   prof <- prof[[1]]
  #   prof_lag <- as.vector(dplyr::lag(prof))
  #   prof_dif <- prof - prof_lag
  #   channels[r] <- sum(prof_dif == 1, na.rm = T) # sum all cases where dif = 1 i.e. new channel edge
  #   # (will only count 1 bank; dif == -1 for other bank.)
  # }
  # dividedness[i, 1] <- t
  # dividedness[i, 2:3] <- channels
  
  # save output raster
  binary_name <- paste0("./wet_binary_smoothed/wbs", t, ".tif")
  writeRaster(f2, filename = binary_name, format = "GTiff", overwrite=TRUE)
  print(paste0(t, " = done!"))
}

# dividedness <- data.frame(dividedness)
# names(dividedness)[2:3] <- as.character(arc)
# names(dividedness)[1] <- "time"

# write.csv(dividedness, 
#           paste0("./t_series/Run", run, "_rep", rep, "_flow_dividedness.csv"), 
#           row.names = FALSE)
