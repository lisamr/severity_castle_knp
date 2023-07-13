
# for plotting posteriors of w as a function of distance
f_plot_w <- function(wmatrix, ...){
  plot(NA, ...)
  for(i in 1:100){
    lines(dist_vec_sim, wmatrix[i,], col= scales::alpha(1, .2))
  }
}
f_wmax <- function(wmatrix){
  wmax <- dist_vec_sim[apply(w_rep, 1, which.max)]
  return(list(dens(wmax), HPDI(wmax, .9)))
}
f_w90 <- function(wmatrix, ...){
  w90 <- map_vec(1:nrow(wmatrix), function(i){
    idx_90th <- which(cumsum(wmatrix[i,]) >= .90)[1]
    dist_vec_sim[idx_90th]
  })
  cumsum_m <- apply(wmatrix, 1, cumsum) %>% t
  
  cumsum_m[1,] %>% 
    plot(dist_vec_sim, ., type = 'l', col= scales::alpha(1, .2), ...)
  for(i in 1:50){
    cumsum_m[i,] %>% lines(dist_vec_sim, ., col= scales::alpha(1, .2)) 
  }
  lines(dist_vec_sim, colMeans(cumsum_m), lwd =3)
  abline(h = .9, lty = 2, lwd = 2, col = 'red')
  segments(x0 = HPDI(w90, .9)[1], x1 = HPDI(w90, .9)[2], y0 = .9, y1 = .9,
           lwd = 4, col = 'red4')
  points(x = mean(w90), y = .9, pch = 16, cex = 2, col = 'brown3')
  
  HPDI(w90, .9)
}
f_plot_cov <- function(delta, dist_vec_sim, ...){
  cov_fx <- sapply(1:nrow(delta), function(i){
    exp(-dist_vec_sim^2 / (2*(delta[i]^2)) )
  }) %>% t
  plot(NA, ylim = c(0, 1), ...)
  for(i in 1:50){
    lines(dist_vec_sim, cov_fx[i,], col= scales::alpha(1, .2))
  }
  colMeans(cov_fx)
  lines(dist_vec_sim, colMeans(cov_fx), lwd =3)
}


f_ddvar_wide <- function(df, variable){
  # converts ring df with distance dependent to data wide. N plots (rows) x M rings (cols)
  variable <- as.name(variable) # unquotes it
  dd_var <- df %>% # dat$rings
    mutate(var_z = zscore({{variable}})) %>% 
    pivot_wider(id_cols = id, names_from = radius, values_from = var_z, names_prefix = 'X_') %>% 
    arrange(id) %>% 
    select(-id) %>% 
    as.matrix()
  return(dd_var)
}