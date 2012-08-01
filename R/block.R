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
  gp.names <- unique(data[,groups])
  
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
    if(is.null(valid.var)){
      valid <- 0
      validvar <- numeric(1)
      validlb <- numeric(1)
      validub <- numeric(1)
    }
    else{
      valid <- 1
      validvar <- data.gp[,valid.var]
      validlb <- valid.range[1]
      validub <- valid.range[2]
    }
    if(algorithm != "optimal"){
    if(is.character(distance)){
      if(algorithm == "optGreedy"){
        out1 <- mahaloptgreed(data.gp,
                              block.vars,
                              vcov=vc.all,
                              n.tr=n.tr,
                              l2=level.two,
                              l1names=level.one.names,
                              valid=as.integer(valid),
                              validvar = as.double(validvar),
                              validlb = as.double(validlb),
                              validub = as.double(validub),
                              verbose=as.integer(verbose))
      }
      else if(algorithm  %in%   c("naiveGreedy", "randGreedy", "sortGreedy")){
        out1 <- mahalnaive(x= data.gp,
                           block.vars=block.vars,
                           vcov=vc.all,
                           n.tr=n.tr,
                           l2=level.two,
                           l1names=level.one.names,
                           valid=as.integer(valid),
                           validvar = as.double(validvar),
                           validlb = as.double(validlb),
                           validub = as.double(validub),
                           verbose=as.integer(verbose))
      }
    }
    else{
      if(algorithm == "optGreedy"){
      dist.mat <- distance[data[,groups]==i, data[,groups]==i]
      out1 <- optgreed(dist=dist.mat,
                       n.tr = n.tr,
                       l2=level.two,
                       l1names=level.one.names,
                       valid=as.integer(valid),
                       validvar = as.double(validvar),
                       validlb = as.double(validlb),
                       validub = as.double(validub),
                       verbose=as.integer(verbose))
    }
      else if(algorithm  %in%   c("naiveGreedy", "randGreedy", "sortGreedy")){
        dist.mat <- distance[data[,groups]==i, data[,groups]==i]
      out1 <- naive(dist=dist.mat,
                    n.tr = n.tr,
                    l2=level.two,
                    l1names=level.one.names,
                    valid=as.integer(valid),
                    validvar = as.double(validvar),
                    validlb = as.double(validlb),
                    validub = as.double(validub),
                    verbose=as.integer(verbose))
    }
    }
  }
    
    
    if(algorithm == "optimal"){
    	 if(n.tr > 2){
    	  warning("You specified algorithm = optimal and n.tr > 2.  However, optimal blocking only implemented for exactly two treatment conditions.  If no other error is encountered, optimal blocks for n.tr = 2 are returned here.")
    	 }
      if(is.character(distance)){
        dist.mat <- mahal(data.block, vc.all)
      }else{      
        dist.mat <- distance[data[,groups]==i, data[,groups]==i]
      }

      if(!is.null(valid.var)){
        d.mat <- expand.grid(data.block[, valid.var], data.block[, valid.var])
        diffs <- abs(d.mat[,1]-d.mat[,2])
        valid.vec <- (valid.range[1] <= diffs) & (diffs <= valid.range[2])
      }

      dist.mat <- matrix(as.integer(optfactor*dist.mat),
                         nrow=nrow(dist.mat),
                         ncol=ncol(dist.mat))

      if(!is.null(valid.var)){
      	warning("You specified algorithm = optimal and valid.var.  However, valid.var and valid.range are only implemented for other algorithms.  If no other error is encountered, optimal blocks ignoring the restriction are returned here.")
        ##dist.mat[!valid.vec] <- 2147483647 #maximum 32 bit integer
      }

      optimalDistanceMatrix <- distancematrix(dist.mat)
      optimalOutput <- nonbimatch(optimalDistanceMatrix, precision = 9)
      optimalOutput$halves$Distance <- as.double(optimalOutput$halves$Distance)/optfactor
                                        #      optimalOutput$halves[, 1] <- level.one.names[optimalOutput$halves[, 1]] 
                                        #      optimalOutput$halves[, 2] <- level.one.names[optimalOutput$halves[, 2]]
      out1 <- optimalOutput$halves
    }
    
    storage1 <- out1
    storage1[storage1==0 & col(storage1) < ncol(storage1)] <- NA
#    storage1[storage1[,1:(ncol(storage1)-1)]==0, 1:(ncol(storage1)-1)] <- NA
    count <- 1

    if(algorithm != "optimal"){
      for(i in 1:(ncol(out1) -1)){
        storage1$temp <- as.character(data.gp[storage1[,i], id.vars[1]])
        storage1$temp2 <- as.character(data.gp[storage1[,i], id.vars[length(id.vars)]])
        names(storage1)[ncol(out1) + count] <- paste("Unit", i)
        count <- count + 1
        names(storage1)[ncol(out1) + count] <- paste("Subunit", i)
        count <- count + 1
      }
    }

    else if(algorithm == "optimal"){
      for(i in 1:(ncol(out1) -1)){
        names(storage1)[i] <- paste("Unit", i)
        storage1[,i] <-  as.character(data.gp[storage1[,i], id.vars[1]])
      }
    }
    if(algorithm != "optimal"){
      storage1$Distance <- storage1[,ncol(out1)]
      storage <- storage1[,(ncol(out1)+1):ncol(storage1)]
    }
    else if(algorithm == "optimal"){
      storage <- storage1
    }
    if(n.tr > 2){
      names(storage)[ncol(storage)] <- "Max Distance"
    }
    odd.col <- seq(1:ncol(storage))[seq(1:ncol(storage))%%2 == 1]
    even.col <- (odd.col+1)[1:(length(odd.col)-1)]
    if(level.two == FALSE && algorithm != "optimal"){
      storage <- storage[, odd.col]
    }
    if(nrow(storage) > 1){
      ## sort by max distance
      o <- order(as.numeric(as.character(storage[,ncol(storage)])),
                 storage[,1])
      storage <- data.frame(storage[o,], check.names = FALSE)
    }
    
    ## function to count NA, to remove empty rows (for valid.var)
    sum.na <- function(sum.na.vector){return(sum(is.na(sum.na.vector)))}
    ## remove empty rows
    storage <- storage[apply(storage, 1, sum.na) != (ncol(storage)-1),]
    
    rownames(storage) <- 1:(nrow(storage))
    out[[gp]] <- storage
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
