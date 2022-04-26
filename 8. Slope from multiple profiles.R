# Script to estimate fan slope from set of 20 down-fan profiles.
# November 2019, sometime.

# Housekeeping
rm(list = ls())

#library(lidR)
library(rgeos)
library(SDraw)
library(raster)
library(dplyr)

#--------------------------------------------------------------
# specify drive, run, repe and number of profiles
run <- 2
repe <- 1
drive <- "L:"
nprof <- 90

# set working directory
setwd(paste0(drive, "/Experiments/Processing/Run", run, "/Run", run, 
             "_rep", repe, "/_1min_intervals"))

# specify max elevation (removes dangling orange cable in corner)
#z_max <- -0.1

# file lists:
m_list <- list.files("./DoD_binary", pattern = ".tif")
# f_list <- list.files("./pointclouds", pattern = ".las")
# f_list <- f_list[(length(f_list) - length(m_list) + 1):length(f_list)]
d_list <- list.files(paste0("./DEM"), pattern = ".tif|.asc")
d_list <- d_list[(length(d_list) - length(m_list) + 1):length(d_list)]

# work out profile lines:
# Get quarter-circle curve that approximates fan toe
rad <- 3 # max poss radius of fan.
x_line <- seq(0, rad, 0.0001)
y_line <- sqrt(rad ^ 2 - x_line ^ 2)
cds <- cbind(x_line, y_line)
l <- spLines(cds)
#lines(l)

l_length <- (pi * 2 * rad / 4) # length of arc line

# interpolate evenly spaced pts long that line
pts <- rgeos::gInterpolate(l, seq(0.01, (l_length - 0.01), length.out = nprof), 
                    normalized = FALSE)
#plot(pts)

# vector of evenly spaced dist vals for naming each profile:
prof_names <- paste0("p", as.character(seq(1, nprof, by = 1)))

# generate profile between 0.01, 0.01 and each pt on line; save to list
l_list <- list()
for (j in 1:nprof){
  line.end <- pts@coords[j, ]
  cds <- rbind(c(0.01, 0.01), line.end) # coords of line start and end
  rad_prof <- spLines(cds)
  #lines(rad_prof)
  l_pts <- rgeos::gInterpolate(rad_prof, d = seq(0, rad-0.01, 0.001), normalized = FALSE) # extract pts at 1mm intervals
  #points(l_pts, col = "red")
  l_list[[j]] <- l_pts
}

# empty global-vectors to store outputs:
slope_mean    <- numeric(length(m_list))
length_mean   <- numeric(length(m_list))
curv_old_mean <- numeric(length(m_list))
curv_new_mean <- numeric(length(m_list))
curv_new_norm_mean <- numeric(length(m_list))

slope_st_dev    <- numeric(length(m_list))
length_st_dev   <- numeric(length(m_list))
curv_old_st_dev <- numeric(length(m_list))
curv_new_st_dev <- numeric(length(m_list))
curv_new_norm_st_dev <- numeric(length(m_list))

slope_upper    <- numeric(length(m_list))
length_upper   <- numeric(length(m_list))
curv_old_upper <- numeric(length(m_list))
curv_new_upper <- numeric(length(m_list))
curv_new_norm_upper <- numeric(length(m_list))

slope_lower    <- numeric(length(m_list))
length_lower   <- numeric(length(m_list))
curv_old_lower <- numeric(length(m_list))
curv_new_lower <- numeric(length(m_list))
curv_new_norm_lower <- numeric(length(m_list))

slope_2pc    <- numeric(length(m_list))
length_2pc   <- numeric(length(m_list))
curv_old_2pc <- numeric(length(m_list))
curv_new_2pc <- numeric(length(m_list))
curv_new_norm_2pc <- numeric(length(m_list))

slope_16pc    <- numeric(length(m_list))
length_16pc   <- numeric(length(m_list))
curv_old_16pc <- numeric(length(m_list))
curv_new_16pc <- numeric(length(m_list))
curv_new_norm_16pc <- numeric(length(m_list))

slope_50pc    <- numeric(length(m_list))
length_50pc   <- numeric(length(m_list))
curv_old_50pc <- numeric(length(m_list))
curv_new_50pc <- numeric(length(m_list))
curv_new_norm_50pc <- numeric(length(m_list))

slope_84pc    <- numeric(length(m_list))
length_84pc   <- numeric(length(m_list))
curv_old_84pc <- numeric(length(m_list))
curv_new_84pc <- numeric(length(m_list))
curv_new_norm_84pc <- numeric(length(m_list))

slope_98pc    <- numeric(length(m_list))
length_98pc   <- numeric(length(m_list))
curv_old_98pc <- numeric(length(m_list))
curv_new_98pc <- numeric(length(m_list))
curv_new_norm_98pc <- numeric(length(m_list))

tvec <- numeric(length(m_list))

# for each timestep:
#for(i in 1) {
for(i in 1:length(m_list)) {
  
  m <- m_list[i] # mask file
  mask_polygon <- raster(paste0("./DoD_binary/", m)) # mask as raster
  # mask_polygon <- aggregate(mask_polygon, fact = 10) # coarsen mask
  # mask_polygon <- rasterToPolygons(mask_polygon,
  #                                fun=NULL, 
  #                                n=4, 
  #                                na.rm=TRUE, 
  #                                digits=12, 
  #                                dissolve=TRUE) # vectorize
  # mask_polygon2<- SpatialPolygons(mask_polygon@polygons)
  
  # # import LAS as lasCatalog object (quicker than loading whole LAS in memory)
  # f <- f_list[i]
  # las_F <- readLAScatalog(paste0("./pointclouds/", f), progress = TRUE)
  # opt_filter(las_F) <- paste0("-drop_z_above ",  z_max) # elevation filter
  
  # import DEM
  d <- d_list[i]
  DEM <- raster(paste0("./DEM/", d)) %>% # mask as raster
    `*` (mask_polygon)
  #plot(DEM)
  
  t <- paste0(strsplit(d,'')[[1]][4:7], collapse = '') # get time, as character
  tvec[i] <- t
  
  # vector of 1 mm spaced x-vals for profiles:
  x.vec <- seq(0.01, rad, by = 0.001)
  df <- data.frame(x_coord = x.vec) # data frame to save profiles
  
  # empty t-step-level vectors to store outputs:
  slope <- numeric(nprof)
  prof_length <- numeric(nprof)
  curv_old <- numeric(nprof)
  curv_new <- numeric(nprof)
  curv_new_norm <- numeric(nprof)
  
  # for each profile
  #for (k in 30){ # turns out first and last profiles have no data, due to cropped DEM borders.
    
  for (k in 2:(nprof - 1)){ # turns out first and last profiles have no data, due to cropped DEM borders.
    
    # get profile
    p <- l_list[[k]]
    # points(p)
    
    # extract 1 mm spaced points
    p_3d <- as.vector(raster::extract(DEM, p))
    
    p_df <- data.frame(x = x.vec, y = p_3d) %>% # write to df with x coord
      filter(complete.cases(.)) # get rid of rows with NA.
    # plot(p_df$x, p_df$y, type = "l")
    
    # save profile to data frame, padding with NA
    df$new_y <- c(p_df$y, rep(NA, (length(x.vec) - length(p_df$y))))
    names(df)[names(df) == "new_y"] <- prof_names[k]
    
    # code for profiles from las, not DEM
    #------------------ 
    
    # # buffer around line (0.5 mm)
    # rad_prof_buff <- buffer(p, width=0.0005)
    
    # # intersect with fan mask to get just the area on the fan
    # rad_prof_buff <- raster::intersect(mask_polygon2, rad_prof_buff)
    # 
    # # extract points in buffer zone
    # rad_prof_buff_3d <- lasclip(las_F, rad_prof_buff)
    # rad_prof_buff_3d <- as.spatial(rad_prof_buff_3d)
    # #plot(rad_prof_buff_3d)
    # 
    # # add column with dist from (0.01, 0.01) for each pt
    # rad_prof_buff_3d$x.dist <- pointDistance(c(0.01, 0.01), 
    #                                          rad_prof_buff_3d, 
    #                                          lonlat = FALSE)
    # rad_prof_buff_3d <- rad_prof_buff_3d[order(rad_prof_buff_3d$x.dist), ] # order by x coord
    # 
    # # plot
    # plot(rad_prof_buff_3d$x.dist, rad_prof_buff_3d$Z, type = "l", 
    #      main = paste0("k = ", k))
    # 
    # # convert to 1 pt per mm
    # # first truncate x.vec:
    # x.out <- x.vec[x.vec <= max(rad_prof_buff_3d$x.dist)]
    # y.vec <- approx(rad_prof_buff_3d$x.dist, rad_prof_buff_3d$Z, 
    #                 xout = x.out, method = "linear")
    # plot(y.vec, type = "l")
    # 
    # # save profile to data frame
    # df$new_y <- c(y.vec$y, rep(NA, (length(x.vec) - length(y.vec$y))))
    # names(df)[names(df) == "new_y"] <- prof_names[k]
    
    #--------------
    # get profile parameters; save to empty t-step-level vectors
    # length
    prof_length[k] <- max(p_df$x)
    
    # old slope
    #-----------------
    # old slope: highest point - lowest point, div by length along x
    # smooth first: average of current sample, 6 future samples, 
    #   and 6 past samples
    # f11 <- rep(1/11,11)
    # smoo_pro <- stats::filter(p_df$y, filter = f11, 
    #                           method = "convolution",
    #                           sides = 2,
    #                           circular = FALSE) # smoothed profile, from moving average with window 25 (must be odd)
    # 
    # smoo_pro_df <- data.frame(x = p_df$x, y = smoo_pro) %>%
    #   filter(complete.cases(.))
    #lines(smoo_pro_df$x, smoo_pro_df$y, col = "green")
    
    # rise <- smoo_pro_df$y[1] - smoo_pro_df$y[length(smoo_pro_df$y)] # first y - last y
    # run <- max(smoo_pro_df$x) - min(smoo_pro_df$x) # don't count from zero - there are some NA vals
    # slope[k] <- rise/run
    #----------------
    
    # new slope: lm across profile; get slope from model object.
    mod <- lm(y ~ x, data = p_df)
    #abline(mod)
    slope[k] <- mod$coefficients[2]
    
    # curv_old: 
    #--------------
    # Use "sinuosity" i.e. length of profile line divided by length of
    #     straight line w/ same start and end 
    w_size <- 201 # window size
    filt <- rep(1/w_size, w_size) # create smoothing filter
    smoo_pro_2 <- stats::filter(p_df$y, filter = filt, # apply filter
                                method = "convolution",
                                sides = 2,
                                circular = FALSE)
    
    # create straight line with same start and end points
    smoo_pro_2 <- as.vector(smoo_pro_2) # convert to vector;
    smoo_pro_mat <- cbind(p_df$x, smoo_pro_2) %>% # as matrix of xy coords
      na.omit(.)
    smoo_pro_spat <- spLines(smoo_pro_mat) # as spatial object, to get length
    
    line_y <- c(max(smoo_pro_mat[, 2]), min(smoo_pro_mat[, 2])) # y coords of straight line from apex to toe
    line_x <- c(min(smoo_pro_mat[, 1]), max(smoo_pro_mat[, 1])) # x coords of line
    cds2 <- cbind(line_x, line_y) # as matrix
    l2 <- spLines(cds2) # as spatial object to get length
    
    # plot(p_df$x, smoo_pro_2, main = k, type = "l") # plot to check all in order
    # lines(p_df$x, p_df$y)
    # lines(cds2, col = "blue")
    
    leng_prof <- lineLength(smoo_pro_spat, byid = FALSE)
    leng_l2 <- lineLength(l2, byid = FALSE)
    curv_old[k] <- leng_prof / leng_l2
    
    # curv new
    #-----------
    # interpolate l2 at 1 mm intervals
    x.out <- smoo_pro_mat[, 1]
    l2_interp <- approx(line_x, line_y, 
                        xout = x.out, method = "linear")
    
    # plot
    # plot(l2_interp, type = "l", col = "grey",
    #      xlab = "Dist. downfan (m)",
    #      ylab = "Rel. Elevation (m)")
    # lines(smoo_pro_mat)
    
    # subtract from smoothed profile
    curv <- smoo_pro_mat[, 2] - l2_interp$y

    # save curvature to vector
    if (abs(min(curv)) < max(curv)) {
      curv_new[k] <- max(curv)
    } else if (abs(min(curv)) > max(curv)) {
      curv_new[k] <- min(curv)
    } else {
      curv_new[k] <- paste0(min(curv), " but min = max!")
    } 
    
    # curv new, normalized by profile length
    #--------------------
    curv_new_norm[k] <- curv_new[k] / prof_length[k]
  }
  
  # ammend profile database so that rows with only NA are excluded
  df_crop <- df[rowSums(is.na(df[2:(nprof-1)])) != ncol(df[2:(nprof-1)]),] # tricky indexing so extra-long x_coord is ignored
  # export profiles
  write.csv(df_crop, file = paste0("./t_series/Recalculating_slope/profiles/Profiles_t", t, ".csv"), row.names = F)
  
  # save vectors to t-step-level data frame
  # export data frame with nprof-2 rows, 5 cols for each t-step
  prof_names <- prof_names[2:(nprof-1)]
  slope <- slope[2:(nprof-1)]
  prof_length <- prof_length[2:(nprof-1)]
  curv_new <- curv_new[2:(nprof-1)]
  curv_new_norm <- curv_new_norm[2:(nprof-1)]
  curv_old <- curv_old[2:(nprof-1)]
  
  df_summary <- data.frame(prof_names, prof_length, slope, curv_old, curv_new, curv_new_norm)
  write.csv(df_summary, 
            file = paste0("./t_series/Recalculating_slope/profile_summaries/Summary_t", t, ".csv"), row.names = F)
  
  # work out average of characteristics for profiles 
  slope_mean[i] <- mean(slope)
  slope_st_dev[i] <- sd(slope)
  slope_upper[i] <- max(slope)
  slope_lower[i] <- min(slope)
  slope_2pc[i] <- quantile(slope, probs = 0.02, na.rm = T, names = F, type = 8)
  slope_16pc[i] <- quantile(slope, probs = 0.16, na.rm = T, names = F, type = 8)
  slope_50pc[i] <- quantile(slope, probs = 0.5, na.rm = T, names = F, type = 8)
  slope_84pc[i] <- quantile(slope, probs = 0.84, na.rm = T, names = F, type = 8)
  slope_98pc[i] <- quantile(slope, probs = 0.98, na.rm = T, names = F, type = 8)
  
  length_mean[i]   <- mean(prof_length)
  length_st_dev[i] <- sd(prof_length)
  length_upper[i] <- max(prof_length)
  length_lower[i] <- min(prof_length)
  length_2pc[i] <- quantile(prof_length, probs = 0.02, na.rm = T, names = F, type = 8)
  length_16pc[i] <- quantile(prof_length, probs = 0.16, na.rm = T, names = F, type = 8)
  length_50pc[i] <- quantile(prof_length, probs = 0.5, na.rm = T, names = F, type = 8)
  length_84pc[i] <- quantile(prof_length, probs = 0.84, na.rm = T, names = F, type = 8)
  length_98pc[i] <- quantile(prof_length, probs = 0.98, na.rm = T, names = F, type = 8)
  
  curv_old_mean[i] <- mean(curv_old)
  curv_old_st_dev[i] <- sd(curv_old)
  curv_old_upper[i] <- max(curv_old)
  curv_old_lower[i] <- min(curv_old)
  curv_old_2pc[i] <- quantile(curv_old, probs = 0.02, na.rm = T, names = F, type = 8)
  curv_old_16pc[i] <- quantile(curv_old, probs = 0.16, na.rm = T, names = F, type = 8)
  curv_old_50pc[i] <- quantile(curv_old, probs = 0.5, na.rm = T, names = F, type = 8)
  curv_old_84pc[i] <- quantile(curv_old, probs = 0.84, na.rm = T, names = F, type = 8)
  curv_old_98pc[i] <- quantile(curv_old, probs = 0.98, na.rm = T, names = F, type = 8)
  
  curv_new_mean[i] <- mean(curv_new)
  curv_new_st_dev[i] <- sd(curv_new)
  curv_new_upper[i] <- max(curv_new)
  curv_new_lower[i] <- min(curv_new)
  curv_new_2pc[i] <- quantile(curv_new, probs = 0.02, na.rm = T, names = F, type = 8)
  curv_new_16pc[i] <- quantile(curv_new, probs = 0.16, na.rm = T, names = F, type = 8)
  curv_new_50pc[i] <- quantile(curv_new, probs = 0.5, na.rm = T, names = F, type = 8)
  curv_new_84pc[i] <- quantile(curv_new, probs = 0.84, na.rm = T, names = F, type = 8)
  curv_new_98pc[i] <- quantile(curv_new, probs = 0.98, na.rm = T, names = F, type = 8)
  
  curv_new_norm_mean[i] <- mean(curv_new_norm)
  curv_new_norm_st_dev[i] <- sd(curv_new_norm)
  curv_new_norm_upper[i] <- max(curv_new_norm)
  curv_new_norm_lower[i] <- min(curv_new_norm)
  curv_new_norm_2pc[i] <- quantile(curv_new_norm, probs = 0.02, na.rm = T, names = F, type = 8)
  curv_new_norm_16pc[i] <- quantile(curv_new_norm, probs = 0.16, na.rm = T, names = F, type = 8)
  curv_new_norm_50pc[i] <- quantile(curv_new_norm, probs = 0.5, na.rm = T, names = F, type = 8)
  curv_new_norm_84pc[i] <- quantile(curv_new_norm, probs = 0.84, na.rm = T, names = F, type = 8)
  curv_new_norm_98pc[i] <- quantile(curv_new_norm, probs = 0.98, na.rm = T, names = F, type = 8)
  
  print(paste0(t, " is done!"))
}

# convert global-level vectors to data frame
df_t_series <- data.frame(tvec, slope_mean, slope_st_dev, slope_upper, slope_lower, 
                          slope_2pc, slope_16pc, slope_50pc, slope_84pc, slope_98pc, 
                          #
                          length_mean, length_st_dev, length_upper, length_lower, 
                          length_2pc, length_16pc, length_50pc, length_84pc, length_98pc,
                          #
                          curv_old_mean, curv_old_st_dev, curv_old_upper, curv_old_lower, 
                          curv_old_2pc, curv_old_16pc, curv_old_50pc, curv_old_84pc, curv_old_98pc,
                          #
                          curv_new_mean, curv_new_st_dev, curv_new_upper, curv_new_lower, 
                          curv_new_2pc, curv_new_16pc, curv_new_50pc, curv_new_84pc, curv_new_98pc,
                          #
                          curv_new_norm_mean, curv_new_norm_st_dev, curv_new_norm_upper, curv_new_norm_lower, 
                          curv_new_norm_2pc, curv_new_norm_16pc, curv_new_norm_50pc, 
                          curv_new_norm_84pc, curv_new_norm_98pc)
# save data frame as csv
write.csv(df_t_series, 
          file = paste0("./t_series/Recalculating_slope/Slope_length_curv.csv"), row.names = F)
