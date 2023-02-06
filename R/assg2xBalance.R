assg2xBalance <- function(assg.obj, data, id.var, bal.vars, to.report = "all"){
  
  data[["Tr"]] <- extract_conditions(assg.obj, data, id.var)
    
  fff <- formula(paste("Tr ~ ", paste(bal.vars, collapse = "+")))
  xbal.list <- list()
  n.groups <- length(assg.obj$assg)
  
  for(i in 1:n.groups){
    
    # This group's assignments:
    assg.gp <- assg.obj$assg[[i]]
    
    # The data from this group:
    data.tr.gp <- data[data[[id.var]] %in% unlist(assg.gp[, 1:(ncol(assg.gp) - 1)]), ]
    
    xbal.out <- RItools::xBalance(fff, data = data.tr.gp, report = c(to.report))

    xbal.list[[paste("Group", i, sep = "")]] <- xbal.out
  }
  
  data <- data[!(is.na(data$Tr)), ]
  xbal.list[["Overall"]] <- RItools::xBalance(fff, data = data, report = c(to.report))

  return(xbal.list)	
}