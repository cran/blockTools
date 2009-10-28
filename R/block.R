block <- function(data, vcov.data = NULL, groups = NULL, n.tr = 2, id.vars,
                  block.vars = NULL, algorithm = "optGreedy", distance =
                  "mahalanobis", weight = NULL, optfactor = 10^8, row.sort = NULL, 
                  level.two = FALSE, valid.var = NULL, valid.range = NULL, seed, 
                  verbose = FALSE, ...){ 

  if(is.null(algorithm)){
    stop("Blocking algorithm is unspecified.  See documentation for
options.")
  }
  
  if(algorithm == "randGreedy"){    
    row.sort <- sample(seq(1:nrow(data)))
    data <- data[row.sort,]
  }
  
  if(algorithm == "sortGreedy"){
    if(is.null(row.sort)){
      stop("Blocking algorithm is 'sortGreedy', but vector of new row
positions unspecified.  See documentation for details.")      
    }
    if(length(row.sort)!= nrow(data)){
    stop("Length of vector 'row.sort' does not equal number of rows in
the data.  Respecify 'row.sort'.")
    }
    data <- data[row.sort,]
  }
     
  if(is.matrix(data)){
    data <- as.data.frame(data)
  }

  if(is.null(block.vars)){
    block.vars <- names(data)[!(names(data) %in% id.vars)]
  }
  
  ## subset appropriate columns for vcov calculation
  if(is.null(vcov.data)){  
    vcov.data <- data[, block.vars]
  }else{
    vcov.data <- vcov.data[, block.vars]
  }

  ## calculate variance using all groups' units
  if(is.character(distance)){
    if(distance == "mahalanobis"){
      vc.all <- var(vcov.data)
    }
    if(distance == "mcd"){
      vc.all <- cov.rob(vcov.data, method="mcd", seed = seed, ...)$cov
    }
    if(distance == "mve"){
      vc.all <- cov.rob(vcov.data, method="mve", seed = seed, ...)$cov
    }
  }

  if(!is.null(weight)){
    if(is.vector(weight)){
      if(length(weight)!=ncol(vc.all)){
  		   stop("Weight vector length must be equal to number of blocking variables.  Respecify 'weight'.")
  		 } 		  		
  		weight <- diag(weight)
    }
  	 if(is.matrix(weight)){
      if(sum(dim(weight)==dim(vc.all)) != 2){
      	stop("Weight matrix dimensions must equal number of blocking variables.  Respecify 'weight'.")
  		 }
    }
  		
  	 vc.all <- solve(t(solve(t(chol(vc.all))))%*%weight%*%solve(t(chol(vc.all))))
  }

  ## counter for groups
  gp <- 0
  
  ## list of output for every group
  out <- list()  

  if(is.null(groups)){
    data[,"groups"] <- 1 
    groups <- "groups"
  }

  if(is.factor(data[,groups])){
    data[,groups] <- as.character(data[,groups])
  }
  gp.names <- rep(NA, length(unique(data[,groups])))

  ## perform blocking w/in groups
  for(i in unique(data[,groups])){ 

    gp <- gp + 1  

    if(verbose == TRUE){
      cat("Blocking group ", i, "\n")
    }

    data.gp <- data[data[,groups]==i, c(id.vars, block.vars)]

    level.one.names <- data.gp[, id.vars[1]] 

    if(level.two == TRUE){
      if(length(id.vars) < 2){
        stop("Blocking requested at second level, but second level not
identified.  Specify a second ID variable and re-block.")
      }
      row.names(data.gp) <- data.gp[, id.vars[2]]
    }else{
      if(length(unique(data.gp[,id.vars[1]])) !=
         length(data.gp[,id.vars[1]])){
        stop("Blocking requested at top level, but some units have
identical values of the identification variable.  Respecify first
identification variable and re-block.")
      }      
      row.names(data.gp) <- data.gp[, id.vars[1]]  
    }

    data.block <- data.frame(data.gp[, !(names(data.gp) %in% id.vars)])

    ## create distance matrix if not user specified
    if(is.character(distance)){
      dist.mat <- mahal(data.block, vc.all)
    }else{      
      dist.mat <- distance[data[,groups]==i, data[,groups]==i]
    }

    row.names(dist.mat) <- row.names(data.block)
    colnames(dist.mat) <- 1:(nrow(data.block)) 

    if(!is.null(valid.var)){
      d.mat <- expand.grid(data.block[, valid.var], data.block[, valid.var])
      diffs <- abs(d.mat[,1]-d.mat[,2])
      valid.vec <- valid.range[1] <= diffs & diffs <= valid.range[2]
      if(algorithm != "optimal"){
      	dist.mat[!valid.vec] <- Inf
      }
      if(algorithm == "optimal"){
    		dist.mat[!valid.vec] <- 99999*max(dist.mat)
    		}
    }

    ## use only half the distance matrix (non-redundant)
    if(algorithm != "optimal"){
    	 dist.mat[row(dist.mat) <= col(dist.mat)] <- Inf
	 }
	      
    if(level.two == TRUE && (length(level.one.names)>1)){
      ## prohibit level 2 units from matching w/in own level 1 unit
      for(qq in 1:(length(level.one.names)-1)){
        for(rr in (qq+1):(length(level.one.names))){
          if(level.one.names[qq]==level.one.names[rr]){
            dist.mat[rr,qq] <- Inf
          }
        }
      }
    }

    ## store level 1, level 2 names, max block distance
    storage <- as.data.frame(as.matrix(t(rep(NA, 2*n.tr+1))))
    
    odd.col <- seq(1:ncol(storage))[seq(1:ncol(storage))%%2 == 1]
    even.col <- (odd.col+1)[1:(length(odd.col)-1)]

    ## BEGIN optimal BYPASS HERE
    if(algorithm != "optimal"){

      ## block counter
      j <- 0  

      ## begin block loop
      while(sum(dist.mat==Inf)< prod(dim(dist.mat))){
        j <- j+1

        if(verbose == TRUE){
          cat("Block number is", j, "\n")
        }

        if(j>1){storage <- rbind(storage, rep(NA, ncol(storage)))}

        ## initialize use.dist
        if(nrow(storage)==1){
          use.dist <- dist.mat
        }
      
        while(sum(use.dist==Inf) < prod(dim(use.dist))){

          if(algorithm=="optGreedy"){
            addr <- optGreedy(use.dist)
          } 
          if(algorithm %in% c("naiveGreedy", "randGreedy", "sortGreedy")){
          	 addr <- naiveGreedy(use.dist)
          }        

          gr <- addr$gr
          gc <- addr$gc

          both <- c(row.names(use.dist)[gr],row.names(use.dist)[gc])
          both <- both[!(both %in% storage[j,])]
          both <- both[!(both %in% as.matrix(storage[,even.col]))]      

          if(sum(is.na(storage[j,]))== length(storage[j,])){
            store.cols.l1 <- c(1,3)
            store.cols.l2 <- c(2,4)
          }else{
            store.cols.l1 <- sum(!is.na(storage[j,]))
            store.cols.l2 <- sum(!is.na(storage[j,]))+1
          }

          store.cols.l1 <- unique(store.cols.l1)
          store.cols.l2 <- unique(store.cols.l2)

          ## store names
          storage[j, store.cols.l2] <- both     

          ## create level one names
          both2 <- c(as.character(level.one.names[gr]),
                   as.character(level.one.names[gc])) 

          if(sum(both2 %in% as.matrix(storage[,odd.col]))!=0){
            both2 <- both2[!(both2 %in% as.matrix(storage[,odd.col]))]
          }

          ## store level one names
          storage[j, store.cols.l1] <- both2

          ## make submatrix of distances for this block
          m1 <- matrix(!(row(dist.mat) %in% which(rownames(dist.mat) %in%
                                               as.matrix(storage))), 
                     nrow(dist.mat),ncol(dist.mat))
          m2 <- matrix(!col(dist.mat) %in% which(rownames(dist.mat) %in%
                                               as.matrix(storage)),
                     nrow(dist.mat),ncol(dist.mat))
          submat <- dist.mat
          submat[m1 & m2] <- Inf

          ## create matrix of all chosen units (level one or two)
          already.mat <- permutations(length(which(rownames(dist.mat)
                                                   %in%
                                                   as.matrix(storage))), 2,
                                      which(rownames(dist.mat) %in%
                                            as.matrix(storage)))
        
          ## store max distance in block
          storage[j, ncol(storage)] <-
            max(dist.mat[already.mat][dist.mat[already.mat]<Inf])

          if(level.two == TRUE){
            ## all level two units from chosen level one units
            new.l2 <- which(rownames(dist.mat) %in%
                          as.character(data[,id.vars[2]][data[,id.vars[1]] %in% as.matrix(storage[,odd.col[1:(length(odd.col)-1)]])]))   
          
            new.l2 <- new.l2[!(new.l2 %in% already.mat[,])]
            submat[new.l2,] <- submat[,new.l2] <- Inf
          }

          for(alr.row in 1:nrow(already.mat)){
            submat[already.mat[alr.row,1], already.mat[alr.row,2]] <-Inf
          }
        
          use.dist <- submat

          if(sum(is.na(storage[j,]))==0){
            use.dist <- Inf
          }
        }

        if(nrow(as.matrix(dist.mat))>1){
          dist.mat <- dist.mat[-c(which(level.one.names %in%
                                      storage[j,odd.col])),
                             -c(which(level.one.names %in%
                                      storage[j,odd.col]))]
          level.one.names <- level.one.names[-c(which(level.one.names
                                                    %in%
                                                    storage[j,odd.col]))]
          use.dist <- dist.mat
        }
      }
    } ## END optimal BYPASS HERE.

    if(algorithm == "optimal"){
	   dist.mat <- optfactor*dist.mat
      optimalDistanceMatrix <- distancematrix(dist.mat)      
      optimalOutput <- nonbimatch(optimalDistanceMatrix, precision = 9)
      optimalOutput$halves$Distance <- as.double(optimalOutput$halves$Distance)/optfactor
#      optimalOutput$halves[, 1] <- level.one.names[optimalOutput$halves[, 1]] 
#      optimalOutput$halves[, 2] <- level.one.names[optimalOutput$halves[, 2]] 
	   storage <- optimalOutput$halves
	 }

    if(nrow(storage) > 1){
      ## sort by max distance
      o <- order(as.numeric(as.character(storage[,ncol(storage)])),
                 storage[,1])
      storage <- data.frame(storage[o,])
    }

    storage <- as.data.frame(storage)

    ## name columns
    names(storage)[odd.col] <- paste("Unit", 1:length(odd.col))
    names(storage)[even.col] <- paste("Subunit", 1:length(even.col))

    names(storage)[ncol(storage)] <- "Distance"
    if(n.tr>2){
      names(storage)[ncol(storage)] <- "Max Distance"
    }

    ## cut duplicate level one names when level.two==F
    if(level.two == FALSE){
      storage <- storage[, odd.col]
    }

    storage[,ncol(storage)] <-
      as.numeric(as.character(storage[,ncol(storage)]))

    rownames(storage) <- 1:(nrow(storage))
    out[[gp]] <- storage
    gp.names[gp] <- i
  }

  names(out) <- gp.names
  ## sort out by group names
  o <- order(names(out))
  out <- out[o]

  output <- list(blocks = out, level.two = level.two)
  output$call <- match.call()  
  class(output) <- "block"
  return(output)
}