subset_brick <- function(brick, dates){
  brick[[match(dates, terra::time(brick))]]
}


# sample raster at those points at that time
sample_brick <- function(pts, brick){
  pts_unproj <- pts %>% arrange(date_adj) %>% st_transform(crs(brick))
  date_var <- time(brick)
  col_idx <- match(pts_unproj$date_adj, date_var)
  
  # extract all data. rows=locations, cols = dates
  s <- terra::extract(brick, pts_unproj, ID = F, layer = col_idx)
  
  return(s$value)
}
