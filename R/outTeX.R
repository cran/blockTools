outTeX <- function(block.obj, namesCol = NULL, digits = 2, ...){

  ## takes block, assignment, or diagnose object
  if(!is.null(block.obj$blocks)){ 
    block.obj <- block.obj$blocks
  }
  if(!is.null(block.obj$assg)){ 
    block.obj <- block.obj$assg
  }
  
  for(i in 1:length(block.obj)){
    tab <- block.obj[[i]]
    nm <- names(block.obj)[i]
    lab <- paste("group.", nm, sep = "")
    cap <- paste("Group ", nm, ".", sep="")

    ncol.tab <- ncol(tab)
    
    ## user-specified column names
    if(!is.null(namesCol)){
      names(tab) <- namesCol
    }
    
    tab.tex <- xtable(tab, label = lab, caption = cap, align =
                      c(rep("c",ncol(tab)+1)), 
                      digits = rep(digits, ncol(tab)+1), ...)
    print(tab.tex, file = paste("Group", nm, ".tex", sep=""),
          caption.placement = "bottom")
  }
}
