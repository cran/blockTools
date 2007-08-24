optGreedy <- function(dist.mat){

  gr <- (which.min(dist.mat)%%nrow(dist.mat))

  ## random tie-breaking
  if(length(gr)>1){  
    gr <- sample(gr,1)
  }
  
  if(gr == 0){
    gr <- nrow(dist.mat)
  }
  
  gc <- unname(which.min(dist.mat[gr,]))

  ## complete random tie-breaking
  if(length(gc)>1){  
    gc <- sample(gc,1)
  }
  
  addr <- list(gr=gr, gc=gc)
  return(addr)
}
