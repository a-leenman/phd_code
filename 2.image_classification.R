##################################################################################
# Script to classify orthomosaic TIF files from alluvial fan table experiments.  #
#                                                                                #
# Anya Leenman                                                                   #
# 8 Feb 2019
# Editted 5 June 2019 for streamlining
##################################################################################
##################################################################################
### Housekeeping
rm(list = ls())

run <- 2 # run number
repe <- 1 # Experimental repeat number
drive <- "L:" # harddrive to write/read to
first_index <- 2 # index of first orthomosaic. Normaly 5. 
# Adjust so f_list starts at same t-value as m_list.

# Set working directory
setwd(paste0(drive, "/Experiments/Processing/Run", run, "/Run", run, "_rep", repe, "/_1min_intervals"))

# load required packages
library(raster)
library(rgdal)
library(magrittr) # for pipes
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

##############################################################################
# set up loop
m_list <- list.files(paste0("./DoD_binary"), pattern = "tif")
f_list <- list.files(paste0("./orthomosaics"), pattern = "tif")
f_list <- f_list[first_index:length(f_list)] # tweak so that both lists start at same value of t

# empty vectors for area calcs:
t_v        <- numeric(length(f_list))
wet_v      <- numeric(length(f_list))
f_area_v   <- numeric(length(f_list))
wet_frac_v <- numeric(length(f_list))
threshold  <- numeric(length(f_list))

for (i in 1:length(f_list)){ 
  # set time
  t <- paste0(strsplit(f_list[i],'')[[1]][11:14], collapse = '') # get time, as character
  t_v[i] <- t
  
  # Masking
  #-------------------------
  # read in mask of fan area
  m <- paste0("./DoD_binary/", m_list[i])
  msk <- raster(m)
  
  # read in orthomosaic tif; mask it
  f <- paste0("./orthomosaics/", f_list[i])
  tf <- brick(f) %>% # import file w/ all 3 bands in single object
    `*` (msk) # multiply by mask to extract values only on fan.
  
  # Classify based on RGB values
  #-------------------------------
  combo <- (tf[[2]] + tf[[3]] - tf[[1]]) / (tf[[1]] + tf[[2]] + tf[[3]])
  thresh <- mean(combo@data@values, na.rm = TRUE) * 1.05 # set threshold as 5% above mean
  combo <- combo > thresh 
  combo[combo != 1] <- NA
  combo <- areaSieve(combo, thresh = 0.001, directions = 8, verbose = FALSE)  # apply area sieve to remove "noise
  
  # write output
  combo_name <- paste0("./wet_binary/wet_binary", t)
  writeRaster(combo, filename = combo_name, format = "GTiff", overwrite=TRUE)
  
  # area calculations
  #-------------------------------
  # total wet area
  wet <- cellStats(combo, "sum") * 0.000001 # constant is conversion to m^2 (each cell is 0.001 * 0.001 m in area)
  wet_v[i] <- wet # write output at correct element of vector
  
  # total fan area (fan mask)
  total <- cellStats(msk, "sum") * 0.000001
  f_area_v[i] <- total
  
  # wet fraction
  wet_fr <- wet/total
  wet_frac_v[i] <- wet_fr
  
  # threshold
  threshold[i] <- thresh
  
  # par(mar = c(2, 2, 2, 2), mfrow = c(1, 2))
  # plot(combo, main = t)
  # plotRGB(tf, mar = TRUE)
  print(paste0(t, " is done!"))
}

# write output vectors to dataframe
t_series <- data.frame("time" = t_v, "wet area" = wet_v, "fan area" = f_area_v, "wet fraction" = wet_frac_v, "threshold" = threshold)
write.csv(t_series, file = paste0("./t_series/Run", run, "_rep", repe, ".csv"))

#################################################################################
# Some other (mostly rejected) options to classify based on RGB

# Red:Green ratio
#red_green_ratio <- tf[[1]] / tf[[3]] 

# Red * Green
#red_green_product <- tf[[1]] * tf[[3]]

# log red
#red_log <- log(tf[[1]])

# red exp
#red_exp <- tf[[1]] ^ 2

# rgb product
#rgb_prod <- tf[[1]] * tf[[2]] * tf[[3]]

# rgb division
#rgb_ratio <- tf[[1]] / tf[[2]] / tf[[3]]

# log blue
#blue_log <- log(tf[[3]])

# red threshold - this seems to work well.
#red <- tf[[1]] < 110
#plot(red)
#red[red == 0] <- NA
#plot(red)

# Vegetation Index Upadhyay et al 2016 : https://ieeexplore.ieee.org/document/7724745
# (blue + red - green) / (blue + red + green)
#VI <- (tf[[3]] + tf[[1]] - tf[[2]]) / (tf[[3]] + tf[[1]] + tf[[2]])
#VI_thresh <- VI < 0.3
#plot(VI_thresh) # Good result.

# Water index from Upadhyay et al 2016
# (Red - blue + green) / (blue + red + green)
#WI <- (tf[[1]] - tf[[3]] + tf[[2]]) / (tf[[3]] + tf[[1]] + tf[[2]])
#WI_thresh <- WI < 0.38
#plot(WI_thresh) # Not so great - doesn't pick up channel on left side of fan.

# Normalized red threshold also good: 
# red / (red + green + blue)

# blue clustering
# blue / (red _ green _ blue)
#blue <- tf[[3]] / (tf[[1]] + tf[[2]] + tf[[3]])
#blue <- blue > 0.303 # set threshold
#blue[blue != 1] <- NA
#blue <- areaSieve(blue, thresh = 0.001, directions = 8, verbose = FALSE) 

#green <- tf[[2]] / (tf[[1]] + tf[[2]] + tf[[3]])
#green <- green > 0.36
#green[green != 1] <- NA
#green <- areaSieve(green, thresh = 0.001, directions = 8, verbose = FALSE) 


# based on these results, best identification at high values of blue and green (norm), low values of red (norm).
# (blue + green - red) / (blue + green + red)
#combo <- (tf[[2]] + tf[[3]] - tf[[1]]) / (tf[[1]] + tf[[2]] + tf[[3]])
#thresh <- mean(combo@data@values, na.rm = TRUE)
#combo <- combo > thresh # set threshold
#combo[combo != 1] <- NA
#combo <- areaSieve(combo, thresh = 0.001, directions = 8, verbose = FALSE) 

# Amount of dye changes through time - didn't want to use supervised machine learning as training sample 
# would have to be redrawn everytime new dye added.
