outCSV <- function(block.obj, namesCol = NULL, digits = 2, ...){

  ## takes block, assignment, or diagnose object
  if(!is.null(block.obj$blocks)){ 
    block.obj <- block.obj$blocks
  }
  
  for(i in 1:length(block.obj)){
    tab <- block.obj[[i]]
    nm <- names(block.obj)[i]
    tab[,ncol(tab)] <- 
    tab[,ncol(tab)] <- round(tab[,ncol(tab)], digits)

    if(!is.null(namesCol)){
      names(tab) <- namesCol
    }
    write.csv(tab, file=paste("Group", nm, ".csv", sep=""), ...)
  }
}
