# Script to calculate overlap (or lack thereof) between successive binary maps of channel. 
# Uses smoothed binary maps as input.
#
# Anya Leenman
# 29 May 2019
#
#######################################################################################
# Housekeeping
rm(list = ls())

run <- 2 # run number
repe <- 1 # Experimental repeat number
drive <- "L:" # drive to write/read to

# Set working directory
setwd(paste0(drive, "/Experiments/Processing/Run", run, "/Run", run, "_rep", repe, "/_1min_intervals"))

# load packages
library(raster)
library(dplyr)
library(lubridate)

#######################################################################################
# load data
f_list <- list.files("./wet_binary_smoothed", pattern = ".tif")
area_list <- read.csv(paste0("./t_series/Run", run, "_rep", repe, ".csv")) %>%
  dplyr::select(-X)

for( j in 1) {
#for( j in 1:20) {
  lag_val <- j

  # Specify indices to start/end at. 
  first_index <- lag_val + 1 # Must be >= lag + 1.
  final_index <- length(f_list) # must be <= length(f_list)
  
  # empty vectors to store loop outputs
  #---------------------------------------
  t_char <- numeric(final_index) # time as character
  
  new_area_raw <- numeric(final_index) # new flow area at t2 (in square m)
  new_area_normby_t2_flow <- numeric(final_index) # new flow area at t2,
  # normalized by flow area at t2.
  new_area_normby_t2_fan <- numeric(final_index) # new flow area at t2, 
  # normalized by fan area at t2.
  
  abandoned_area_raw <- numeric(final_index) # area of flow at t1 that's abandoned at t2 (in sq m)
  abandoned_area_normby_t1_flow <- numeric(final_index) # abandoned flow area at t1, 
  # normalized by flow area at t1.
  abandoned_area_normby_t2_fan <- numeric(final_index) # abandoned flow area at t1,
  # normalized by fan area at t1.
  
  ov_raw <- numeric(final_index) # area of flow overlap between t1 and t2 (in sq m)
  ov_normby_t2_flow <- numeric(final_index) # overlap, normalized by flow area at t2
  ov_normby_t2_fan <- numeric(final_index) # overlap, normalized by fan area at t2.
  #----------------------------
  
  for (i in first_index:final_index){ # first index should be 2 (working out diff from 1-2.)
    # time wrangling
    t_char[i] <- paste0(strsplit(f_list[i],'')[[1]][4:7], collapse = '') # time, as character 
    
    # # trap to catch outlier at t = 1221
    # if (paste0(strsplit(f_list[i],'')[[1]][4:7], collapse = '') == "1221" | 
    #     paste0(strsplit(f_list[(i - lag_val)],'')[[1]][4:7], collapse = '') == "1221") {
    #   next
    # }
    
    # set up t1 raster
    f_t1 <- raster(paste0("./wet_binary_smoothed/", f_list[i - lag_val])) # wet binary at t1
    #f_t1[is.na(f_t1)] <- 0 # use for layer w/ data saved as NA = non-channel.
    
    f_t2 <- raster(paste0("./wet_binary_smoothed/", f_list[i])) # wet binary at t2
    f_t2[f_t2 == 1] <- 2 # set channel to "2" for raster math
    #f_t2[is.na(f_t2)] <- 0 # use for layer w/ data saved as NA = non-channel.
    
    f_ov <- f_t1 + f_t2 # overlap matrix. 
    # key: 1 = old channel only. 2 = new channel. 3 = overlap
    #plot(f_ov)
    
    # case 1/3/4: newly created flow area
    # raw:
    # Code looks messy but is more efficient than running a frequency count on f_ov.
    new_area <- f_ov # create raster copy
    new_area[new_area != 2] <- NA
    new_area[new_area == 2] <- 1 # select only values of interest; change to 1 for summation
    new_area_raw[i] <- cellStats(new_area, "sum") * 0.000001 # constant for area conversion
    
    # normalized by t2 flow area:
    f_t2[f_t2 == 2] <- 1 # change to 1 for summation; set values for next iteration
    new_area_normby_t2_flow[i] <- cellStats(new_area, "sum") / cellStats(f_t2, "sum")
    
    # normalized by t2 fan area
    new_area_normby_t2_fan[i] <- cellStats(new_area, "sum") * 0.000001 / 
      area_list$fan.area[area_list$time == as.vector(t_char[i])]
    
    # case 2: one channel thread captures all flow
    # i.e. former channel is abandoned
    # raw
    abandoned_area <- f_ov
    abandoned_area[abandoned_area != 1] <- NA
    abandoned_area_raw[i] <- cellStats(abandoned_area, "sum") * 0.000001
    
    # normalized by t1 flow area
    abandoned_area_normby_t1_flow[i] <- cellStats(abandoned_area, "sum") / cellStats(f_t1, "sum")
    # normalized by t2 fan area
    abandoned_area_normby_t2_fan[i] <- cellStats(abandoned_area, "sum") * 0.000001 / 
      area_list$fan.area[area_list$time == as.vector(t_char[i])]
    
    overlap <- f_ov
    overlap[overlap != 3] <- NA
    overlap[overlap == 3] <- 1
    ov_raw[i] <- cellStats(overlap, "sum") * 0.000001 
    # constant is conversion to m^2 (each cell is 0.001 * 0.001 m in area)
    ov_normby_t2_flow[i] <- cellStats(overlap, "sum") / cellStats(f_t2, "sum") # total cells in overlap divided by 
    # total cells in t1 (i.e. overlap as fraction of initial channel pattern)
    ov_normby_t2_fan[i] <- cellStats(overlap, "sum") * 0.000001 / 
      area_list$fan.area[area_list$time == as.vector(t_char[i])]
    
    # error catch
    if (ov_normby_t2_flow[i] + new_area_normby_t2_flow[i] != 1){
      print("Yikes, flows don't add up")
      break
    }
    
    print(t_char[i])
  }
  
  # write output vectors to dataframe
  t_series <- data.frame("time" = as.numeric(t_char), new_area_raw, new_area_normby_t2_flow, new_area_normby_t2_fan,
                         abandoned_area_raw, abandoned_area_normby_t1_flow, abandoned_area_normby_t2_fan,
                         ov_raw, ov_normby_t2_flow, ov_normby_t2_fan)
  t_series[t_series == 0] <- NA # change no-data value to NA
  
  # convert time character to period
  
  t_series <- mutate(t_series, t_h = ifelse((time %% 100) == 0, time / 100, floor (time / 100))) %>% # get hour
    mutate(t_m = ifelse((time %% 100) == 0, 0, time - t_h * 100)) %>% # get minute
    filter(complete.cases(.)) %>% # get rid of rows with NA.
    mutate(time = paste(t_h, t_m, sep = ":")) %>% # concatenate to hh:mm format
    mutate(time = hm(time)) %>% # convert to time format using lubridate
    mutate(t_sec = period_to_seconds(time)) %>% # Convert to number of seconds since start
    mutate(t_min = t_sec / 60) %>% # time in minutes
    mutate(t_hr = t_min / 60) %>% # time in hours
    dplyr::select(-t_h, -t_m) # remove unnecessary columns
  
  
  write.csv(t_series, file = paste0("./t_series/change_detection/Run", run, "_rep", repe, "_lag", lag_val, "_overlap_and_antioverlap.csv"), row.names = FALSE)
  
  print(paste0("lag ", j, " is complete!"))
}
# plotting
#plot(t_series$t_hr, t_series$new_area_normby_t2_fan, type = "l", ylim = c(0.03, 0.09))
#lines(t_series$t_hr, t_series$abandoned_area_normby_t2_fan, col = "Red")

