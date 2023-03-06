
# create concentric rings around location i
f_rings <- function(loc, radii, crs){
  # test inputs
  stopifnot(length(radii) > 1)
  stopifnot(is.matrix(loc) | is.data.frame(loc))
  stopifnot(c('x', 'y') %in% colnames(loc))
  if(!grepl('UTM', crs$input)) 
  warnings("crs isn't in UTM coords. not tested on unprojected data.")
  
  
  # create concentic rings around c(0,0)
  geometry <- st_point(c(0,0)) %>% 
    st_sfc()
  sample_point <- st_sf(geometry)
  circles <- map(radii, ~ st_buffer(sample_point, .x)) 
  rings <- bind_rows(
    circles[[1]],
    map_df(2:length(circles), function(i){
      suppressWarnings(
        st_difference(circles[[i]], circles[[i-1]]$geometry)
      )
    })
  ) %>% 
    mutate(radius = radii, .before = geometry) 
  
  # copy and move those rings to all the locations
  sfpoints <- st_as_sf(as.data.frame(loc), coords = c('x', 'y'))
  sfpoints <- mutate(sfpoints, id = 1:nrow(sfpoints), .before = geometry)
  
  res <- map_df(split(sfpoints, 1:nrow(sfpoints)), function(x) {
    rings_copy <- rings
    rings_copy$geometry <- rings$geometry + x$geometry
    rings_copy <- mutate(rings_copy, id = x$id, .before = radius)
    st_crs(rings_copy) <- crs
    return(rings_copy)
  })
  return(res)
}
