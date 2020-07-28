#function to get delay and embedding
library(tseriesChaos)
estimate_del_emb <- function(Data, max_delay, delay_inc = 1, max_dim, dim_inc = 1, fnn_thresh = .1, theiler_window = 1, max_lost = NA) {
  #Delay
  #get local minimum of average mutual information
  ami <- mutual(as.matrix(Data), lag.max = max_delay, plot = FALSE) 
  d <- as.numeric(unname(which(diff(sign(diff(ami)))==2)))
  if (length(d) == 0) {
    del <- max_delay
  } else {
    del <- d[1]
  }
  #get first value under false-nn threshold
  if(!is.na(max_lost)) {
    m <- floor((max_lost+del) / del) #embedding cap to not lose more than max_lost rows of data
    if (m < max_dim) {
      max_dim <- m
    }
  }
  e_thresh <- fnn_thresh
  e <- unname(false.nearest(as.matrix(Data), m = max_dim, d = del, t = theiler_window)[1,]) #false nearest neighbors
  e_inc <- seq(0, length(e), dim_inc)
  e_val <- e[e_inc]
  e_val <- e_val[!is.na(e_val)]
  emb <- as.numeric(which(e_val < e_thresh)[1]) #first value under threshold
  
  if (length(emb) == 0 | is.na(emb)) { #if none under threshold, get emb value corresponding to smallest FNN
    emb <- which(e_val == min(e_val, na.rm = TRUE))
  }
  emb <- emb[1]
  
  return(list(del = del, emb = emb))
}