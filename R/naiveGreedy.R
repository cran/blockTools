naiveGreedy <- function(dist.mat){

  ## if first column is all Inf, take first nonInf col:
  if(infmean(dist.mat[,1])==1){
    gc <- unname(which.min(apply(dist.mat,2,infmean)))[1]
  }else{
    gc <- 1
  }

  gr <- unname(which.min(dist.mat[,gc]))

  ## random tie-breaking within column
  if(length(gr)>1){  
    gr <- sample(gr,1)
  }
  
  addr <- list(gr=gr, gc=gc)
  return(addr)
}
