assignment <- function(block.obj, seed = NULL, namesCol = NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }

  if(is.matrix(block.obj) || is.data.frame(block.obj)){
    tmp <- list()
    tmp$blocks$"1" <- block.obj
    block.obj <- tmp
  }

  if(is.null(block.obj$level.two)){
    block.obj$level.two <- FALSE
  }

  out <- list()
  gp.names <- array(NA)

  ## perform assignment w/in groups
  for(i in 1:length(block.obj$blocks)){ 

    gp.obj <- as.matrix(block.obj$blocks[[i]])

    ncol.tab <- ncol(gp.obj)
    
    if(is.null(namesCol)){
      namesCol <- c(rep(NA, ncol.tab-1), "Distance")
      if((ncol.tab>5) || ((ncol.tab)>3 && (block.obj$level.two ==
                                            FALSE))){
        namesCol[length(namesCol)] <- "Max Distance"
      }
      
      if(block.obj$level.two == FALSE){
        for(j in 1:(ncol.tab-1)){
          namesCol[j] <- paste("Treatment ", j, sep = "")
        }
      }else{
        for(j in 1:((ncol.tab-1)/2)){
          namesCol[(2*j-1):(2*j)] <- rep(paste("Treatment ", j, sep = ""),2)
        }
      }
    }
              
    ## Put units into treatment groups with pr(u_i in g_j) = 1/|g|    
    for(j in 1:(nrow(gp.obj))){
      tmp <- gp.obj[j,]
      if(block.obj$level.two == FALSE){
        tmp[1:(ncol.tab-1)] <- tmp[sample(ncol.tab-1)]
      }else{
        s <- sample((1:(ncol.tab-1))[odd(1:(ncol.tab-1))])
        tmp[1:(ncol.tab-1)] <- tmp[c(rbind(s,s+1))]
      }
      gp.obj[j,] <- tmp
    }
    gp.obj <- as.data.frame(gp.obj)
    gp.obj[,ncol(gp.obj)] <- as.numeric(as.character(gp.obj[,ncol(gp.obj)]))
    names(gp.obj) <- namesCol
    out[[i]] <- gp.obj
    gp.names[i] <- names(block.obj$blocks)[i]
  }

  names(out) <- gp.names

  output <- list(assg = out)
  output$call <- match.call()
  class(output) <- "assg"
  return(output)
}
